import re
import numpy as np
from types import SimpleNamespace
from typing import List, Dict, Any
from src.core import Cell
from src.utils import create_logger, get_project_root
from src.gaussian_process.interpolators import BatteryParameterInterpolator

np.random.seed(1234)


class BasePack:
    """
    Base class providing common utilities for Pack Solvers, including 
    topology parsing, result containers, and logging.
    """
    def __init__(self, topology_str: str, cell_params: Dict[str, Any], n_virtual_packs: int):
        self.logger = create_logger(self.__class__.__name__)
        self.topology_str = topology_str
        self.n_virtual_packs = n_virtual_packs
        self.cell_params = cell_params
        
        # Parse "1-p(3)" into structured data
        self.components = self._parse_topology(topology_str)

        # Common result containers
        self.time = None
        self.voltage = SimpleNamespace(mean=None, std=None, samples=None)
        self.currents = SimpleNamespace(
            samples=[],          # Raw current splits: (time, n_virtual, n_parallel)
            per_block_mean=[],   # Average current per block over time
            per_block_std=[]     # Current distribution (sigma) per block
        )

    def _parse_topology(self, topology: str) -> List[Dict[str, Any]]:
        """Parses the topology string into structured component data."""
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

        self.logger.info(f"Parsed pack topology with following modules in series:")
        for i, comp in enumerate(components):
            if comp['n'] == 1:
                self.logger.info(f"({i}) - {comp['n']} cell in series")
            else:
                self.logger.info(f"({i}) - {comp['n']} cells in {comp['type']}")

        return components

    def _finalize_results(self, pack_v_history: np.ndarray):
        """Aggregates sample history into mean and standard deviation."""
        self.voltage.samples = pack_v_history
        self.voltage.mean = np.mean(pack_v_history, axis=0)
        self.voltage.std = np.std(pack_v_history, axis=0)
        self.logger.info("Simulation completed and results finalized.")


class PCPack(BasePack):
    """
    Polynomial Chaos Expansion based pack solver. 
    Uses a reference cell solve to get statistical profiles.
    """
    def __init__(self, topology_str: str, cell_params: Dict[str, Any], n_virtual_packs: int = 200):
        super().__init__(topology_str, cell_params, n_virtual_packs)
        self.ref_cell = Cell(**self.cell_params)
        self.logger.info(f"Initialized PCPack with {n_virtual_packs} realizations.")

    def solve(self, current: np.ndarray, time: np.ndarray, temperature: np.ndarray) -> None:
        self.logger.info(f"Starting PCPack solver...")
        n_steps = len(time)
        dt = np.diff(time, prepend=time[0])
        self.time = time
        pack_v_history = np.zeros((self.n_virtual_packs, n_steps))

        # Reset current distribution logs
        self.currents.samples, self.currents.per_block_mean, self.currents.per_block_std = [], [], []

        for block_n, comp in enumerate(self.components):
            self.logger.info(f"Solving module {block_n}: {comp['type']} with n={comp['n']}")
            n = comp['n']
            total_instances = self.n_virtual_packs * n
            
            # Solve base cell for nominal stats to feed the interpolator
            i_nominal = current / n if comp['type'] == 'parallel' else current
            self.ref_cell.solve(current=i_nominal, time=time, temperature=temperature)
            
            models = self.ref_cell.get_model_dict()
            ecm, ocv, hyst = models.get('ecm'), models.get('ocv'), models.get('hysteresis')
            
            # Parameter Profile Extraction
            interpolator = BatteryParameterInterpolator()
            dist_params = interpolator.get_interpolated_params(temperature, self.ref_cell.soc.mean)

            # Assign fixed identities (z-scores) for every cell in this block
            z = {k: np.random.normal(0, 1, total_instances) for k in ['r0', 'ocv', 'hyst', 't1', 'c1', 't2', 'c2']}
            
            # Time-varying baseline matrices
            r0_samples = ecm.mu_r0[:, np.newaxis] + (z['r0'][np.newaxis, :] * ecm.sigma_r0[:, np.newaxis])
            v_ocv = ocv.voltage.mean[:, np.newaxis] + (z['ocv'][np.newaxis, :] * ocv.voltage.std[:, np.newaxis])
            v_hyst = hyst.voltage.mean[:, np.newaxis] + (z['hyst'][np.newaxis, :] * hyst.voltage.std[:, np.newaxis]) if hyst else 0
            
            v_rc1, v_rc2 = np.zeros(total_instances), np.zeros(total_instances)
            block_v_history = np.zeros((n_steps, self.n_virtual_packs))
            if comp['type'] == 'parallel':
                comp_current_samples = np.zeros((n_steps, self.n_virtual_packs, n))

            eps = 1e-12
            for k in range(n_steps):
                # Update RC parameters for current SOC/Temp state
                t1_inv = np.maximum(np.abs(dist_params['mu_tau1_inv'][k] + z['t1'] * dist_params['sigma_tau1_inv'][k]), eps)
                c1_inv = dist_params['mu_c1_inv'][k] + z['c1'] * dist_params['sigma_c1_inv'][k]
                t2_inv = np.maximum(np.abs(dist_params['mu_tau2_inv'][k] + z['t2'] * dist_params['sigma_tau2_inv'][k]), eps)
                c2_inv = dist_params['mu_c2_inv'][k] + z['c2'] * dist_params['sigma_c2_inv'][k]

                v_cell_pot = v_ocv[k, :] + v_hyst[k, :] + v_rc1 + v_rc2
                r0_k = r0_samples[k, :]

                if comp['type'] == 'series':
                    i_cells = current[k]
                    v_terminal = v_cell_pot + (i_cells * r0_k)
                    block_v_history[k, :] = np.sum(v_terminal.reshape(self.n_virtual_packs, n), axis=1)
                else:
                    # Parallel Admittance solver: I_total = sum((V_block - V_pot_i) / R0_i)
                    inv_r0 = 1.0 / r0_k
                    v_pot_r = (v_cell_pot * inv_r0).reshape(self.n_virtual_packs, n).sum(axis=1)
                    sum_inv_r = inv_r0.reshape(self.n_virtual_packs, n).sum(axis=1)
                    
                    v_block = (current[k] + v_pot_r) / sum_inv_r
                    block_v_history[k, :] = v_block
                    
                    # Back-calculate current split for state feedback
                    i_cells = (np.repeat(v_block, n) - v_cell_pot) * inv_r0
                    comp_current_samples[k, :, :] = i_cells.reshape(self.n_virtual_packs, n)

                # Discrete-time state update
                e1, e2 = np.exp(-t1_inv * dt[k]), np.exp(-t2_inv * dt[k])
                v_rc1 = v_rc1 * e1 + (i_cells * c1_inv / t1_inv) * (1 - e1)
                v_rc2 = v_rc2 * e2 + (i_cells * c2_inv / t2_inv) * (1 - e2)

            pack_v_history += block_v_history.T
            if comp['type'] == 'parallel':
                self._log_current_dist(comp_current_samples)

        self._finalize_results(pack_v_history)

    def _log_current_dist(self, samples):
        """Helper to store parallel current statistics."""
        self.currents.samples.append(samples)
        self.currents.per_block_mean.append(np.mean(samples, axis=(1, 2)))
        self.currents.per_block_std.append(np.std(samples, axis=(1, 2)))


class MCPack(BasePack):
    """
    Monte Carlo based pack solver. 
    Performs full time-stepping for every cell realization.
    """
    def __init__(self, topology_str: str, cell_params: Dict[str, Any], n_virtual_packs: int = 200):
        super().__init__(topology_str, cell_params, n_virtual_packs)
        self.cell = Cell(**self.cell_params)
        self.interpolator = BatteryParameterInterpolator()
        self.logger.info(f"Initialized MCPack with {n_virtual_packs} realizations.")

    def solve(self, current: np.ndarray, time: np.ndarray, temperature: np.ndarray):
        self.logger.info(f"Starting MCPack solver...")
        n_steps = len(time)
        dt = np.diff(time, prepend=time[0])
        self.time = time
        pack_v_history = np.zeros((self.n_virtual_packs, n_steps))

        # Reset logs
        self.currents.samples, self.currents.per_block_mean, self.currents.per_block_std = [], [], []

        for block_n, comp in enumerate(self.components):
            self.logger.info(f"Solving module {block_n}: {comp['type']} with n={comp['n']}")
            n = comp['n']
            n_total_cells = n * self.n_virtual_packs

            # Nominal solve for SOC
            i_nom = current / n if comp['type'] == 'parallel' else current
            self.cell.solve(current=i_nom, time=time, temperature=temperature)
            
            models = self.cell.get_model_dict()
            ocv, hyst = models.get('ocv'), models.get('hysteresis')
            dist_ecm = self.interpolator.get_interpolated_params(temperature, self.cell.soc.mean)

            # Extract OCV/Hysteresis mean and sigma profiles
            v_ocv_mu, v_ocv_sigma = ocv.solve(soc=[self.cell.soc.mean, self.cell.soc.std], T=temperature)
            v_hyst_mu, v_hyst_sigma = (hyst.solve(current=i_nom, time=time, soc=[self.cell.soc.mean, self.cell.soc.std], T=temperature) 
                                       if hyst else (np.zeros(n_steps), np.zeros(n_steps)))

            # Sample z scores
            z = {k: np.random.normal(0, 1, n_total_cells) for k in ['R0', 'tau1', 'c1', 'tau2', 'c2', 'ocv', 'hyst']}

            # (4) Physics Loop
            v_rc1, v_rc2 = np.zeros(n_total_cells), np.zeros(n_total_cells)
            block_v_history = np.zeros((n_steps, self.n_virtual_packs))
            if comp['type'] == 'parallel':
                comp_current_samples = np.zeros((n_steps, self.n_virtual_packs, n))

            eps = 1e-12
            for k in range(n_steps):
                # Sample parameters at time k
                r0 = dist_ecm['mu_R0'][k] + z['R0'] * dist_ecm['sigma_R0'][k]
                t1_inv = np.maximum(np.abs(dist_ecm['mu_tau1_inv'][k] + z['tau1'] * dist_ecm['sigma_tau1_inv'][k]), eps)
                c1_inv = dist_ecm['mu_c1_inv'][k] + z['c1'] * dist_ecm['sigma_c1_inv'][k]
                t2_inv = np.maximum(np.abs(dist_ecm['mu_tau2_inv'][k] + z['tau2'] * dist_ecm['sigma_tau2_inv'][k]), eps)
                c2_inv = dist_ecm['mu_c2_inv'][k] + z['c2'] * dist_ecm['sigma_c2_inv'][k]

                # Reconstruct source potential
                v_src = (v_ocv_mu[k] + z['ocv'] * v_ocv_sigma[k]) + (v_hyst_mu[k] + z['hyst'] * v_hyst_sigma[k])
                v_cell_pot = v_src + v_rc1 + v_rc2

                if comp['type'] == 'series':
                    i_cells = current[k]
                    v_term = v_cell_pot + i_cells * r0
                    block_v_history[k, :] = np.sum(v_term.reshape(self.n_virtual_packs, n), axis=1)
                else:
                    inv_r0 = 1.0 / r0
                    v_pot_r = (v_cell_pot * inv_r0).reshape(self.n_virtual_packs, n).sum(axis=1)
                    sum_inv_r = inv_r0.reshape(self.n_virtual_packs, n).sum(axis=1)
                    
                    v_block = (current[k] + v_pot_r) / sum_inv_r
                    block_v_history[k, :] = v_block
                    
                    i_cells = (np.repeat(v_block, n) - v_cell_pot) * inv_r0
                    comp_current_samples[k, :, :] = i_cells.reshape(self.n_virtual_packs, n)

                # State feedback
                e1, e2 = np.exp(-t1_inv * dt[k]), np.exp(-t2_inv * dt[k])
                v_rc1 = v_rc1 * e1 + (i_cells * c1_inv / t1_inv) * (1 - e1)
                v_rc2 = v_rc2 * e2 + (i_cells * c2_inv / t2_inv) * (1 - e2)

            pack_v_history += block_v_history.T
            if comp['type'] == 'parallel':
                self.currents.samples.append(comp_current_samples)
                self.currents.per_block_mean.append(np.mean(comp_current_samples, axis=(1, 2)))
                self.currents.per_block_std.append(np.std(comp_current_samples, axis=(1, 2)))

        self._finalize_results(pack_v_history)
