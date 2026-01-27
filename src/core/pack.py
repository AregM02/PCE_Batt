import numpy as np
import re
from typing import List, Dict, Union, Any
from src.core import Cell
from src.utils import create_logger
from types import SimpleNamespace

np.random.seed(1234)

class Pack:
    def __init__(self, topology_str: str, cell_params: Dict[str, Any], n_virtual_packs: int = 200):
        """
        :param topology_str: e.g., "1-p(5)" (1 series + block of 5 parallel)
        :param cell_params: Dictionary containing reference cell parameters
        :param n_virtual_packs: Number of full pack realizations to simulate
        """

        self.logger = create_logger(self.__class__.__name__)
        
        self.topology_str = topology_str
        self.n_virtual_packs = n_virtual_packs
        self.cell_params = cell_params
        self.components = self._parse_topology(topology_str)

        # Statistical cell model
        self.ref_cell = Cell(**self.cell_params)
        self.time = None
        
        # Result containers
        self.voltage = SimpleNamespace(mean=None, std=None, samples=None)
        
        # module-level current and soc distributions
        self.currents = SimpleNamespace(
            samples=[], # for every parallel block: array of shape (time, n_virtual_packs, N_parallel_cells)
            per_block_mean=[],
            per_block_std=[]
        )
        self.soc = SimpleNamespace(divergence=None)
        
        self.logger.info(f"Initialized Pack with topology: {topology_str} and {n_virtual_packs} virtual realizations.")


    def _parse_topology(self, topology: str) -> List[Dict[str, Union[str, int]]]:
        """
        Parses the topology string into structured component data.

        """

        parts = topology.split('-')
        components = []
        for part in parts:
            part = part.strip()
            p_match = re.match(r'p\((\d+)\)', part)
            if p_match:
                components.append({'type': 'parallel', 'n': int(p_match.group(1))})
            else:
                try:
                    components.append({'type': 'series', 'n': int(part)})
                except ValueError:
                    self.logger.error(f"Invalid topology segment: '{part}'")
        return components


    def solve(self, current: np.ndarray, time: np.ndarray, temperature: np.ndarray) -> None:
        """
        Solves the stochastic pack model over the given time series.
        """

        self.logger.info(f"Starting Pack solver for {len(time)} steps.")
        n_steps = len(time)
        pack_v_history = np.zeros((self.n_virtual_packs, n_steps)) # each row represents one random instance of the pack
        self.time = time

        # Reset current lists for new solve
        self.currents.samples = []
        self.currents.per_block_mean = []
        self.currents.per_block_std = []

        for i, comp in enumerate(self.components):
            n = comp['n']
            total_instances = self.n_virtual_packs * n
            self.logger.debug(f"Solving component {i}: {comp['type']} with n={n}")
            
            # (1) Solve base cell for the nominal current share
            i_nominal = current / n if comp['type'] == 'parallel' else current
            self.ref_cell.solve(current=i_nominal, time=time, temperature=temperature)
            
            models = self.ref_cell.get_model_dict()
            ecm = models.get('ecm')
            ocv = models.get('ocv')
            hyst = models.get('hysteresis')

            if not ecm or not ocv:
                self.logger.warning(f"Required models (ECM/OCV) missing for component {i}.")

            # (2) Extract Stochastic Physical Components
            # Here sampled PCE for the transient part only
            ecm_samples = ecm.joint.sample(total_instances)
            v_transient = ecm.vb_transient(*ecm_samples) # (n_steps, total_instances)
            
            # Sample resistance
            z_r0 = np.random.normal(0, 1, total_instances)
            r0_samples = ecm.mu_r0[:, np.newaxis] + (z_r0[np.newaxis, :] * ecm.sigma_r0[:, np.newaxis])

            # Sample OCV & Hysteresis
            z_ocv = np.random.normal(0, 1, total_instances)
            v_ocv = ocv.voltage.mean[:, np.newaxis] + (z_ocv[np.newaxis, :] * ocv.voltage.std[:, np.newaxis])
            
            v_hyst = 0
            if hyst:
                z_hyst = np.random.normal(0, 1, total_instances)
                v_hyst = hyst.voltage.mean[:, np.newaxis] + (z_hyst[np.newaxis, :] * hyst.voltage.std[:, np.newaxis])

            # Combine into the Source Potential (zero-load voltage)
            v_source = v_ocv + v_hyst + v_transient

            # (3) Apply Topology Physics
            if comp['type'] == 'series':
                # If in series, sum all
                v_cell = v_source + (current[:, np.newaxis] * r0_samples)
                v_reshaped = v_cell.reshape(n_steps, self.n_virtual_packs, n)
                pack_v_history += np.sum(v_reshaped, axis=2).T
                
            else:
                # Parallel Admittance Logic: V_block = (I_total + sum(V_src/R0)) / sum(1/R0)
                inv_r = 1.0 / r0_samples
                v_src_reshaped = v_source.reshape(n_steps, self.n_virtual_packs, n)
                inv_r_reshaped = inv_r.reshape(n_steps, self.n_virtual_packs, n)
                
                sum_v_over_r = np.sum(v_src_reshaped * inv_r_reshaped, axis=2)
                sum_inv_r = np.sum(inv_r_reshaped, axis=2)
                
                v_block = (current[:, np.newaxis] + sum_v_over_r) / sum_inv_r
                pack_v_history += v_block.T
                
                # I_i = (V_block - V_src_i) / R0_i
                i_cells = (v_block[:, :, np.newaxis] - v_src_reshaped) * inv_r_reshaped
                
                # Store per-block distribution data
                self.currents.samples.append(i_cells)
                self.currents.per_block_mean.append(np.mean(i_cells, axis=(1, 2)))
                self.currents.per_block_std.append(np.std(i_cells, axis=(1, 2)))

        # (4) Final Aggregation
        self.voltage.samples = pack_v_history
        self.voltage.mean = np.mean(pack_v_history, axis=0)
        self.voltage.std = np.std(pack_v_history, axis=0)
        self.logger.info("Pack simulation completed successfully.")


    def get_weakest_link_capacity(self) -> Dict[str, float]:
        """Calculates pack capacity based on the minimum series cell."""
        c_mean = self.cell_params.get('capacity', 0.0)
        c_std = self.cell_params.get('capacity_unc', 0.0)
        n_series = sum(c['n'] for c in self.components if c['type'] == 'series') or 1
        
        self.logger.debug(f"Calculating weakest link with n_series={n_series}")
        
        samples = np.random.normal(c_mean, c_std, (self.n_virtual_packs, n_series))
        pack_capacities = np.min(samples, axis=1)
        
        return {
            'theoretical_mean': float(c_mean),
            'actual_pack_mean': float(np.mean(pack_capacities)),
            'actual_pack_std': float(np.std(pack_capacities))
        }
