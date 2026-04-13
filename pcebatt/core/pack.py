import numpy as np
import chaospy as cp
from types import SimpleNamespace
from typing import Iterable, List, Dict, Any, Literal
from ..gaussian_process.interpolators import BatteryParameterInterpolator
from ..models.submodels.ocv import qOCV_Differential
import scipy.linalg as la
from ..utils.logger import create_logger


class BasePack:
    def __init__(self, topology: List[Iterable[Any]],
                 cell_params: Dict[str, Any],
                 n_virtual_packs: int = None):
        
        self.n_virtual_packs = n_virtual_packs
        self.cell_params = cell_params
        self.N = 0
        
        self.components = self._parse_topology(topology)

        self.time = None
        self.voltage = SimpleNamespace(mean=None, std=None, samples=None)
        self.currents = SimpleNamespace(samples=None, mean=None, std=None)
        self.soc = SimpleNamespace(samples=None, mean=None, std=None)

    def _parse_topology(self, topology: List[Iterable[Any]]) -> List[int]:
        """ Series/parallel topology is interpreted as a list of parallel blocks"""
        components = []
        for branch in topology:
            self.N += branch[1]
            if branch[0] == "series":
                for _ in range(branch[1]):
                    components.append(1)
            elif branch[0] == "parallel":
                components.append(branch[1])
            else:
                raise Exception("Wrong topology entered!")
                
        return components

    def _finalize_results(self, pack_v_history: np.ndarray, soc_history: np.ndarray, i_history: np.ndarray):
        self.voltage.samples = pack_v_history
        self.voltage.mean = np.mean(pack_v_history, axis=0)
        self.voltage.std = np.std(pack_v_history, axis=0)
        self.soc.samples = soc_history
        self.currents.samples = i_history
        print("Simulation completed and results finalized.")


class MCPack(BasePack):
    def __init__(self,
                 topology: List[Iterable[Any]],
                 cell_params: Dict[str, Any],
                 n_virtual_packs = 200,
                 mode: Literal["thevenin", "const_ecm", "2rc"] = "2rc"):
        
        super().__init__(topology, cell_params, n_virtual_packs)
        self.mode = mode
        self.interpolator = BatteryParameterInterpolator()

    def _generate_default_z(self, size: int) -> Dict[str, np.ndarray]:
        keys = ['R0', 'tau1', 'c1', 'tau2', 'c2', 'cap']
        return {k: np.random.normal(0, 1, size) for k in keys}
    
    def solve(self, current: np.ndarray, time: np.ndarray, temperature: np.ndarray,     
              z_samples: Dict[str, np.ndarray] = None):
        n_steps = len(time)
        dt_vec = np.diff(time, append=time[-1])
        self.time = time
        
        # Pre-instantiate models to avoid overhead inside the loop
        ocv_model = qOCV_Differential(logger_level="ERROR")

        n_v = self.n_virtual_packs
        n_total_samples = self.N * n_v
        
        # Generate or use provided realizations ONCE for all cells
        if z_samples is not None:
            z_master = z_samples
        else:
            np.random.seed(1234) # Seed here to ensure identical configs get identical random arrays
            z_master = self._generate_default_z(n_total_samples)

        total_v_history = np.zeros((n_v, n_steps))
        all_soc_history = []
        all_i_history = []

        cursor = 0
        eps = 1e-12

        # Iterate through each "parallel block" (n=1 is just a series block)
        for n in self.components:
            n_comp_total = n * n_v
            
            # Slice the chunk for this specific component
            z_comp = {k: v[cursor : cursor + n_comp_total] for k, v in z_master.items()}
            cursor += n_comp_total

            # Initialize states for this specific block
            soc = np.full(n_comp_total, self.cell_params['initial_soc'])
            v_rc1, v_rc2 = np.zeros(n_comp_total), np.zeros(n_comp_total)
            q_cells = self.cell_params['capacity'] + z_comp['cap'] * self.cell_params['capacity_unc']

            block_soc_h = np.zeros((n_steps, n_comp_total))
            block_v_h = np.zeros((n_steps, n_v))
            block_i_h = np.zeros((n_steps, n_comp_total))

            for k in range(n_steps):
                block_soc_h[k] = soc
                
                # Parameter Interpolation
                p = self.interpolator.get_interpolated_params(T=temperature[k], SoC=np.mean(soc))
                r0 = np.maximum(p['mu_R0'] + z_comp['R0'] * p['sigma_R0'], eps)
                
                v_ocv, _ = ocv_model.solve(soc=[soc, 0], T=temperature[k])
                v_pot = v_ocv if self.mode == 'thevenin' else (v_ocv + v_rc1 + v_rc2)
                
                # Kirchhoff (handles both n=1 and n>1)
                inv_r = 1.0 / r0
                sum_v_inv_r = np.sum((v_pot * inv_r).reshape(n_v, n), axis=1)
                sum_inv_r = np.sum(inv_r.reshape(n_v, n), axis=1)
                
                v_block = (current[k] + sum_v_inv_r) / sum_inv_r
                block_v_h[k, :] = v_block
                
                # Current flows into cell if v_block > v_pot
                i_cells = (np.repeat(v_block, n) - v_pot) * inv_r
                block_i_h[k] = i_cells

                # State integration (step k -> k+1)
                if k < n_steps - 1:
                    dt_next = dt_vec[k+1]
                    soc = soc + (i_cells * dt_next) / (q_cells * 3600.0)
                    
                    if self.mode != 'thevenin':
                        t1_inv = np.maximum(p['mu_tau1_inv'] + z_comp['tau1'] * p['sigma_tau1_inv'], eps)
                        c1_inv = p['mu_c1_inv'] + z_comp['c1'] * p['sigma_c1_inv']
                        e1 = np.exp(-t1_inv * dt_next)
                        v_rc1 = v_rc1 * e1 + (i_cells * (c1_inv / t1_inv)) * (1.0 - e1)
                        
                        t2_inv = np.maximum(p['mu_tau2_inv'] + z_comp['tau2'] * p['sigma_tau2_inv'], eps)
                        c2_inv = p['mu_c2_inv'] + z_comp['c2'] * p['sigma_c2_inv']
                        e2 = np.exp(-t2_inv * dt_next)
                        v_rc2 = v_rc2 * e2 + (i_cells * (c2_inv / t2_inv)) * (1.0 - e2)

            # Accumulate block voltage onto total pack voltage
            total_v_history += block_v_h.T
            all_soc_history.append(block_soc_h)
            all_i_history.append(block_i_h)

        # Finalize once at the very end
        self._finalize_results(total_v_history,
                               np.concatenate(all_soc_history, axis=1), 
                               np.concatenate(all_i_history, axis=1))
        

class Galerkin_Pack(BasePack):
    def __init__(self,
                 topology: List[Iterable[Any]],
                 cell_params: Dict[str, Any],
                 poly_order: int = 2,
                 mode: Literal["thevenin", "const_ecm", "2rc"] = "2rc",
                 logger_level: str = "INFO"):
        
        super().__init__(topology, cell_params, n_virtual_packs=None)
        self.mode = mode
        self.poly_order = poly_order
        self.interpolator = BatteryParameterInterpolator()
        self.logger = create_logger("GalerkinPack", logger_level)
        
        # 1. Stochastic Setup
        self.logger.info(f"Setting up stochastic space for {self.N} cells...")
        self.dist = cp.J(*[cp.Normal(0, 1) for _ in range(self.N)])
        self.basis = cp.generate_expansion(self.poly_order, self.dist, normed=False)
        self.P = len(self.basis)
        self.logger.info(f"Stochastic space initialized. PCE Basis size (P): {self.P}")

        # 2. Pre-compute
        self._precompute_galerkin_terms()

    def _precompute_galerkin_terms(self):
        self.logger.info("Starting precomputation of Galerkin tensors...")
        nodes, weights = cp.generate_quadrature(self.poly_order * 2, self.dist, rule="gaussian")
        poly_evals = self.basis(*nodes) 
        
        self.norms = np.sum((poly_evals**2) * weights, axis=1)
        
        self.logger.info(f"Computing M-Tensor (size {self.P}^3). This may take a moment...")
        E_mnk = np.einsum('mr, nr, kr, r -> mnk', poly_evals, poly_evals, poly_evals, weights)
        self.M_tensor = E_mnk / self.norms
        
        self.logger.info("Projecting inverse capacity into PCE space...")
        mu_q = self.cell_params['capacity']
        sigma_q = self.cell_params['capacity_unc']
        self.q_inv_coeffs = np.zeros((self.N, self.P))
        
        for i in range(self.N):
            q_realizations = 1.0 / (mu_q + sigma_q * nodes[i, :])
            self.q_inv_coeffs[i, :] = np.sum(q_realizations * poly_evals * weights, axis=1) / self.norms
        self.logger.info("Precomputation finished.")


    def solve(self, current: np.ndarray, time: np.ndarray, temperature: np.ndarray):
        n_steps = len(time)
        dt_vec = np.diff(time, append=time[-1])
        ocv_model = qOCV_Differential(logger_level="ERROR")

        soc = np.zeros((self.N, self.P))
        soc[:, 0] = self.cell_params['initial_soc']
        v_rc1, v_rc2 = np.zeros((self.N, self.P)), np.zeros((self.N, self.P))

        total_v_mean, total_v_var = np.zeros(n_steps), np.zeros(n_steps)
        all_soc_history, all_i_history = np.zeros((n_steps, self.N)), np.zeros((n_steps, self.N))
        eps = 1e-12

        self.logger.info(f"Starting time integration ({n_steps} steps)...")
        log_interval = max(1, n_steps // 10) # Log every 10%

        for k in range(n_steps):
            if k % log_interval == 0:
                self.logger.info(f"Progress: { (k/n_steps)*100:.1f}% | Step {k}/{n_steps}")

            dt = dt_vec[k] if k < n_steps - 1 else 0
            all_soc_history[k] = soc[:, 0]
            
            p = self.interpolator.get_interpolated_params(T=temperature[k], SoC=np.mean(soc[:, 0]))
            r0 = np.maximum(p['mu_R0'], eps)
            
            # Linearized OCV: V = V_mean + (dV/dSoC) * sigma_SoC
            v_mean, dvdz = ocv_model.solve(soc=[soc[:, 0], np.ones(self.N)], T=temperature[k])
            
            v_ocv = np.zeros((self.N, self.P))
            v_ocv[:, 0] = v_mean
            v_ocv[:, 1:] = dvdz[:, None] * soc[:, 1:] 
            
            v_pot = v_ocv if self.mode == 'thevenin' else (v_ocv + v_rc1 + v_rc2)

            cursor, i_cells = 0, np.zeros((self.N, self.P))
            v_pack_block_mean, v_pack_block_var = 0, 0
            
            for n_cells in self.components:
                idx = slice(cursor, cursor + n_cells)
                v_pot_block = v_pot[idx, :] 
                inv_r = 1.0 / r0
                
                current_k = np.zeros(self.P)
                current_k[0] = current[k] 
                
                v_block = (current_k + np.sum(v_pot_block * inv_r, axis=0)) / (n_cells * inv_r)
                
                v_pack_block_mean += v_block[0]
                v_pack_block_var += np.sum((v_block[1:] ** 2) * self.norms[1:])
                i_cells[idx, :] = (v_block[None, :] - v_pot_block) * inv_r
                cursor += n_cells
            
            total_v_mean[k], total_v_var[k] = v_pack_block_mean, v_pack_block_var
            all_i_history[k] = i_cells[:, 0]

            if dt > 0:
                # The Intrusive Galerkin Product for dot_SoC = I * (1/Q)
                dot_soc = np.einsum('mnk, im, in -> ik', self.M_tensor, i_cells, self.q_inv_coeffs) / 3600.0
                soc += dot_soc * dt
                
                if self.mode != 'thevenin':
                    t1_inv, c1_inv = np.maximum(p['mu_tau1_inv'], eps), p['mu_c1_inv']
                    e1 = np.exp(-t1_inv * dt)
                    v_rc1 = v_rc1 * e1 + (i_cells * (c1_inv / t1_inv)) * (1.0 - e1)
                    
                    t2_inv, c2_inv = np.maximum(p['mu_tau2_inv'], eps), p['mu_c2_inv']
                    e2 = np.exp(-t2_inv * dt)
                    v_rc2 = v_rc2 * e2 + (i_cells * (c2_inv / t2_inv)) * (1.0 - e2)

        self.voltage.mean, self.voltage.std = total_v_mean, np.sqrt(total_v_var)
        self.soc.mean, self.currents.mean = all_soc_history, all_i_history
        self.logger.info("Simulation complete.")


class Reduced_Galerkin_Pack(BasePack):
    def __init__(self,
                 topology: List[Iterable[Any]],
                 cell_params: Dict[str, Any],
                 poly_order: int = 2,
                 m_modes: int = 3,
                 mode: Literal["thevenin", "const_ecm", "2rc"] = "2rc",
                 logger_level: str = "INFO"):
        
        super().__init__(topology, cell_params, n_virtual_packs=None)
        self.mode = mode
        self.poly_order = poly_order
        self.m_modes = min(m_modes, self.N) 
        self.interpolator = BatteryParameterInterpolator()
        self.ocv_model = qOCV_Differential(logger_level="ERROR")
        self.logger = create_logger("ReducedGalerkinPack", logger_level)
        
        # base stochastic setup (m modes)
        self.dist = cp.J(*[cp.Normal(0, 1) for _ in range(self.m_modes)])
        self.basis = cp.generate_expansion(self.poly_order, self.dist, normed=False)
        self.P = len(self.basis)
        self.nodes, self.weights = cp.generate_quadrature(self.poly_order * 2, self.dist, rule="gaussian")
        
        # galerkin tensors
        poly_evals = self.basis(*self.nodes) 
        self.norms = np.sum((poly_evals**2) * self.weights, axis=1)
        E_mnk = np.einsum('mr, nr, kr, r -> mnk', poly_evals, poly_evals, poly_evals, self.weights)
        self.M_tensor = E_mnk / self.norms

        # spectral decomp
        init_soc = self.cell_params['initial_soc']
        self.Phi_r = self._compute_spectral_matrix(init_soc)
        
        # capacity uncertainty projection
        self.q_inv_coeffs = self._compute_q_inv_projections(self.Phi_r)
        
        self.logger.info(f"Full Projection RPCE Setup: m={self.m_modes} modes, P={self.P} terms")


    def _compute_spectral_matrix(self, target_soc: float) -> np.ndarray:
        """ Constructs the Block-Diagonal topological mapping matrix. """

        # get baseline slope for the Jacobian
        _, h_vec = self.ocv_model.solve(soc=[np.full(self.N, target_soc), np.ones(self.N)], T=25.0)
        
        Phi_out = np.zeros((self.N, self.m_modes))
        cursor_row, cursor_col = 0, 0
        
        for n_cells in self.components:
            modes_for_block = min(n_cells, self.m_modes - cursor_col)
            if modes_for_block <= 0: break

            L = np.eye(n_cells) - (np.ones((n_cells, n_cells)) / n_cells)
            vals, vecs = la.eigh(L)
            
            idx_mean = np.argmin(np.abs(vals))
            idx_diff = np.argsort(np.abs(vals))[-(modes_for_block-1):] if modes_for_block > 1 else []
            idx_final = np.unique(np.append(idx_mean, idx_diff)).astype(int)
            
            Phi_block = vecs[:, idx_final]
            q, _ = la.qr(Phi_block, mode='economic')
            
            Phi_out[cursor_row : cursor_row + n_cells, cursor_col : cursor_col + modes_for_block] = q
            cursor_row += n_cells
            cursor_col += modes_for_block

        return Phi_out


    def _compute_q_inv_projections(self, Phi_matrix: np.ndarray) -> np.ndarray:
        """ Projects inverse capacity into the reduced stochastic space. """

        poly_evals = self.basis(*self.nodes)
        mapped_nodes = Phi_matrix @ self.nodes 
        
        q_inv_out = np.zeros((self.N, self.P))
        mu_q, sigma_q = self.cell_params['capacity'], self.cell_params['capacity_unc']
        for i in range(self.N):
            q_inv_real = 1.0 / (mu_q + sigma_q * mapped_nodes[i, :])
            q_inv_out[i, :] = np.sum(q_inv_real * poly_evals * self.weights, axis=1) / self.norms
        return q_inv_out


    def _get_ocv_at_nodes(self, z_nodes: np.ndarray, temperature: float) -> np.ndarray:
        """
        Wrapper to evaluate the non-linear OCV model at every quadrature node.
        """
        original_shape = z_nodes.shape
        z_flat = z_nodes.flatten()
        
        # Solve with 0 sigma to get raw OCV values without Taylor expansion
        v_flat, _ = self.ocv_model.solve(soc=[z_flat, np.zeros_like(z_flat)], T=temperature)
        
        return v_flat.reshape(original_shape)


    def solve(self, current: np.ndarray, time: np.ndarray, temperature: np.ndarray):
        n_steps = len(time)
        dt_vec = np.diff(time, append=time[-1])
        eps = 1e-12

        soc = np.zeros((self.N, self.P))
        soc[:, 0] = self.cell_params['initial_soc']
        v_rc1, v_rc2 = np.zeros((self.N, self.P)), np.zeros((self.N, self.P))

        total_v_mean, total_v_var = np.zeros(n_steps), np.zeros(n_steps)
        all_soc_history, all_i_history = np.zeros((n_steps, self.N)), np.zeros((n_steps, self.N))
        
        # pre-evaluate basis for fast projection
        poly_evals = self.basis(*self.nodes) # shape (P, num_nodes)
        
        self.logger.info(f"Solving Full-Projection Reduced System... (Steps: {n_steps})")
        
        for k in range(n_steps):
            dt = dt_vec[k]
            mean_soc = np.mean(soc[:, 0])
            all_soc_history[k] = soc[:, 0]
            
            # reconstruct exact SoC for all N cells at all quadrature nodes
            z_nodes = soc @ poly_evals # shape: (N, num_nodes)
            
            # evaluate exact non-linear OCV
            v_ocv_nodes = self._get_ocv_at_nodes(z_nodes, temperature[k])
            
            # Project the non-linear OCV curve
            v_ocv = (v_ocv_nodes @ (poly_evals.T * self.weights[:, np.newaxis])) / self.norms # shape: (N, P)
            
            # CIRCUIT POTENTIAL
            v_pot = v_ocv if self.mode == 'thevenin' else (v_ocv + v_rc1 + v_rc2)

            # CURRENT DISTRIBUTION
            cursor, i_cells = 0, np.zeros((self.N, self.P))
            v_pack_block_mean, v_pack_block_var = 0, 0
            
            p = self.interpolator.get_interpolated_params(T=temperature[k], SoC=mean_soc)
            current_k = np.zeros(self.P)
            current_k[0] = current[k] 
            inv_r = 1.0 / np.maximum(p['mu_R0'], eps)

            for n_cells in self.components:
                idx = slice(cursor, cursor + n_cells)
                v_pot_block = v_pot[idx, :] 
                
                v_block = (current_k + np.sum(v_pot_block * inv_r, axis=0)) / (n_cells * inv_r)
                
                v_pack_block_mean += v_block[0]
                v_pack_block_var += np.sum((v_block[1:] ** 2) * self.norms[1:])
                i_cells[idx, :] = (v_block[None, :] - v_pot_block) * inv_r
                cursor += n_cells
            
            total_v_mean[k], total_v_var[k] = v_pack_block_mean, v_pack_block_var
            all_i_history[k] = i_cells[:, 0]

            # INTEGRATION
            if dt > 0:
                dot_soc = np.einsum('mnk, im, in -> ik', self.M_tensor, i_cells, self.q_inv_coeffs) / 3600.0
                soc += dot_soc * dt
                
                if self.mode != 'thevenin':
                    for v_rc, t_inv, c_inv in zip([v_rc1, v_rc2], 
                                                 [p['mu_tau1_inv'], p['mu_tau2_inv']], 
                                                 [p['mu_c1_inv'], p['mu_c2_inv']]):
                        e = np.exp(-np.maximum(t_inv, eps) * dt)
                        v_rc[:] = v_rc * e + (i_cells * (c_inv / np.maximum(t_inv, eps))) * (1.0 - e)

        self.voltage.mean, self.voltage.std = total_v_mean, np.sqrt(total_v_var)
        self.soc.mean, self.currents.mean = all_soc_history, all_i_history
        self.logger.info("Simulation Complete.")