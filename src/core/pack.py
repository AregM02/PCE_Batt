import numpy as np
import re
from src.core import Cell
from types import SimpleNamespace

class Pack:
    def __init__(self, topology_str, cell_params, n_virtual_packs=200):
        """
        :param topology_str: e.g., "1-p(5)" (1 series + block of 5 parallel)
        :param n_virtual_packs: Number of full pack realizations to simulate
        """
        self.topology_str = topology_str
        self.n_virtual_packs = n_virtual_packs
        self.components = self._parse_topology(topology_str)
        self.cell = Cell(**cell_params)
        
        # Result container
        self.voltage = SimpleNamespace()

    def _parse_topology(self, topology):
        parts = topology.split('-')
        components = []
        for part in parts:
            part = part.strip()
            p_match = re.match(r'p\((\d+)\)', part)
            if p_match:
                components.append({'type': 'parallel', 'n': int(p_match.group(1))})
            else:
                components.append({'type': 'series', 'n': int(part)})
        return components

    def solve(self, current, time, temperature):
        # Accumulator for pack voltage traces
        # Shape: (n_virtual_packs, time_steps)
        pack_v_history = np.zeros((self.n_virtual_packs, len(time)))

        for comp in self.components:
            n = comp['n']
            # Divide current for parallel blocks, full current for series
            i_cell = current / n if comp['type'] == 'parallel' else current
            
            # 1. Solve the representative base cell
            self.cell.solve(current=i_cell, time=time, temperature=temperature)
            
            # 2. Extract models using the standardized keys
            models = self.cell.get_model_dict()
            ecm = models.get('ecm')
            ocv = models.get('ocv')
            hyst = models.get('hysteresis') 

            # Total number of unique cells needed for this component
            total_instances = self.n_virtual_packs * n
            
            # --- 3. STOCHASTIC REALIZATIONS ---
            
            # A. ECM: Evaluate PCE Polynomial
            # Draw samples from the joint distribution stored in ecm
            ecm_samples = ecm.joint.sample(total_instances)
            v_ecm = ecm.solution_pce(*ecm_samples) # (time, total_instances)
            
            # B. OCV: Gaussian Sampling
            v_ocv = np.random.normal(
                loc=ocv.voltage.mean[:, np.newaxis], 
                scale=ocv.voltage.std[:, np.newaxis], 
                size=(len(time), total_instances)
            )
            
            # C. HYSTERESIS: Gaussian Sampling
            # Assume hysteresis variation is also normally distributed around its mean
            v_hyst = 0
            if hyst:
                v_hyst = np.random.normal(
                    loc=hyst.voltage.mean[:, np.newaxis],
                    scale=hyst.voltage.std[:, np.newaxis],
                    size=(len(time), total_instances)
                )

            # Sum all stochastic contributions for every individual cell
            v_cells_all = v_ecm + v_ocv + v_hyst

            # --- 4. TOPOLOGY AGGREGATION ---
            # Reshape to (time, n_virtual_packs, n_cells_per_pack)
            v_reshaped = v_cells_all.reshape(len(time), self.n_virtual_packs, n)

            if comp['type'] == 'series':
                # Series: Voltages add up
                pack_v_history += np.sum(v_reshaped, axis=2).T
            else:
                # Parallel: Voltages average (Current sharing assumption)
                pack_v_history += np.mean(v_reshaped, axis=2).T

        # --- 5. FINALIZE STATS ---
        self.voltage.samples = pack_v_history
        self.voltage.mean = np.mean(pack_v_history, axis=0)
        self.voltage.std = np.std(pack_v_history, axis=0)