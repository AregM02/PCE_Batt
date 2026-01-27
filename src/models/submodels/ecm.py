import sys
import chaospy
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from scipy.linalg import block_diag
from abc import abstractmethod
from src.utils import create_logger, print_progress
from ..base import Model
from src.gaussian_process.interpolators import BatteryParameterInterpolator
from collections.abc import Iterable
from typing import Tuple
from numpy.typing import NDArray
from types import SimpleNamespace

MAX_SOLVER_STEP = 20.


class ECM(Model):
    """
    Template for an ECM object. Must have a solve method, which must return the simulated voltage array and the standard deviation of its uncertainty. If the model is deterministic, return zeros for the std.
    
    """

    def __init__(self):
        self.interpolator = None  # assigned by subclasses
        self.max_solver_step = MAX_SOLVER_STEP
        self.model_type = 'ecm'
        self.voltage = SimpleNamespace()


    def validate_interpolator(self, required_vars):
        """Check that interpolator provides all required variables."""

        # get one sample of params (e.g. mid-SOC, mid-temp)
        try:
            params = self.interpolator.test_sample()
        except:
            self.logger.warning(f"[{self.__class__.__name__}] Could not sample interpolator during validation.")
            return

        missing = [v for v in required_vars if v not in params.keys()]

        if missing:
            from utils.parameter_loader import load_vars
            import inspect

            raise ValueError(
                f"[{self.__class__.__name__}] Interpolator missing required parameters: {missing}. Adjust {inspect.getfile(load_vars)} to initiate all necessary variables for {self.__class__.__name__}."
            )


    @abstractmethod
    def solve(self, *args, **kwargs) -> Iterable[NDArray, NDArray]:
        pass


class nRC(ECM):
    """Standard implementation of an nRC ECM."""

    def __init__(self, N: int = 2, logger_level: str = 'DEBUG'):
        
        # define required variables
        required_vars = ['R0'] + [var for n in range(1, N+1) for var in (f'tau{n}_inv', f'c{n}_inv')]

        super().__init__()
        self.logger = create_logger(__class__.__name__, logger_level)
        self.interpolator = BatteryParameterInterpolator()
        self.validate_interpolator(required_vars)

        self.logger.info(f'[{__class__.__name__}] Initialized!')


    def solve(self, **kwargs) -> Iterable[NDArray, NDArray]:

        current = kwargs['current']
        time = kwargs['time'] 
        soc = kwargs['soc'][0]
        T = kwargs['T']

        self.logger.info(f'[{__class__.__name__}] Starting solver...')

        # interpolators for input arrays
        current_ip = interp1d(time, current, kind='linear', bounds_error=False, fill_value=(current[0], current[-1]))
        soc_ip = interp1d(time, soc, kind='linear', bounds_error=False, fill_value=(soc[0], soc[-1]))
        T_ip = interp1d(time, T, kind='linear', bounds_error=False, fill_value=(T[0], T[-1]))

        total_duration = time[-1]-time[0]
        # Right-hand-side for the system of differential equations
        def rhs(t:float, x:NDArray) -> NDArray:
            print_progress(t, time[0], total_duration) # show progress

            # Get interpolated soc, current and temperature
            soc_t = soc_ip(t)
            current_t = current_ip(t)
            T_t = T_ip(t)

            # Load and parse interpolated parameters
            params_t = self.interpolator.get_interpolated_params(soc_t, T_t)
            params_t = [params_t[name] for name in params_t.keys()]
            
            # add up contributions from all branches
            x_dot = params_t[0]*current_t
            for tau_inv, c_inv in zip(params_t[1::2], params_t[2::2]):
                x_dot += -tau_inv*x + c_inv*current_t

            return x_dot
    
        voltage = solve_ivp(
                        fun=rhs, max_step=self.max_solver_step, y0=np.array([0.]),
                        t_span=(time[0], time[-1]), t_eval=time
                        ).y.T
        voltage = voltage.flatten() 
        
        self.voltage.mean, self.voltage.std = voltage, np.zeros_like(voltage)

        return voltage, np.zeros_like(voltage)
    

class GalerkinPCE(ECM):
    _instance = None

    def __new__(cls, *args, **kwargs):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self, logger_level: str = 'DEBUG'):
        if not hasattr(self, '_initialized'):

            required_vars = [
                'mu_R0', 'sigma_R0',
                'mu_tau1_inv', 'sigma_tau1_inv', 
                'mu_c1_inv', 'sigma_c1_inv', 
                'mu_tau2_inv', 'sigma_tau2_inv',
                'mu_c2_inv', 'sigma_c2_inv', 
                'CORR',
            ]

            super().__init__()
            self._initialized = True  # Flag to block re-init
            self.logger = create_logger(__class__.__name__, logger_level)
            self.interpolator = BatteryParameterInterpolator() 
            self.validate_interpolator(required_vars)

            # ~~~~ GALERKIN PCE ~~~~~
            # create the joint distribution, basis and parameter variables
            joint = chaospy.J(chaospy.Normal(), chaospy.Normal(), chaospy.Normal(),
                              chaospy.Normal(), chaospy.Normal())
            self.joint = joint

            Phi = chaospy.generate_expansion(3, joint, normed=True)
            self.Phi = Phi

            nu = chaospy.variable(5)
            self.E_Phi = np.array([chaospy.E(Phi, joint)]+[chaospy.E(nu[i]*Phi, joint) for i in range(5)]) # shape = (6, basis_size)
            self.T_Phi = np.array([chaospy.E(chaospy.outer(Phi, Phi), joint)]+
                         [chaospy.E(nu[i]*chaospy.outer(Phi, Phi), joint) for i in range(5)]) # shape = (6, basis_size, basis_size)
            self.nu_r = nu[0] # for R0

            self.output_pce = None

            self.logger.info(f'[{__class__.__name__}] Initialized!')


    def compile_statistical_matrix(self, T: float, soc: float)-> NDArray:
        params_t = self.interpolator.get_interpolated_params(T, soc)
        (
        _, sigma_r0,
        mu_tau1_inv, sigma_tau1_inv, 
        mu_c1_inv, sigma_c1_inv, 
        mu_tau2_inv, sigma_tau2_inv, 
        mu_c2_inv, sigma_c2_inv,
        corr_mat
        ) = (
        params_t['mu_R0'], params_t['sigma_R0'],
        params_t['mu_tau1_inv'], params_t['sigma_tau1_inv'],
        params_t['mu_c1_inv'], params_t['sigma_c1_inv'],
        params_t['mu_tau2_inv'], params_t['sigma_tau2_inv'],
        params_t['mu_c2_inv'], params_t['sigma_c2_inv'],
        params_t['CORR']
        )

        Mu = np.asarray([0, mu_tau1_inv, mu_c1_inv, mu_tau2_inv, mu_c2_inv])
        sigma_diag = np.diag(np.asarray([sigma_r0, sigma_tau1_inv, sigma_c1_inv, sigma_tau2_inv, sigma_c2_inv])) # no correlations for now
        corr_mat = (corr_mat + corr_mat.T)/2  # ensure symmetry
        SIGMA = sigma_diag @ corr_mat @ sigma_diag
        try:
            L = np.linalg.cholesky(SIGMA)  # shape = (5, 5)
        except:
            L = sigma_diag # fallback to no correlations if not positive definite
            
        return np.hstack([Mu[:, None], L]) # shape = (5, 6) [Mu, SIGMA]


    def solve(self, **kwargs) -> Tuple[NDArray, NDArray]:
        """
        Parameters
        ----------
            **kwargs : dict of NDArray
                Input arrays. Keys: 'current', 'time', 'soc', 'T'.
        """
        
        current = kwargs['current']
        time = kwargs['time'] 
        soc = kwargs['soc'][0]
        T = kwargs['T']

        self.logger.info(f'[{__class__.__name__}] Starting solver...')

        # Set up interpolators for input arrays; required for the solver
        current_ip = interp1d(time, current, kind='linear', bounds_error=False, fill_value=(current[0], current[-1]))
        soc_ip = interp1d(time, soc, kind='linear', bounds_error=False, fill_value=(soc[0], soc[-1]))
        T_ip = interp1d(time, T, kind='linear', bounds_error=False, fill_value=(T[0], T[-1]))

        total_duration = time[-1]-time[0]

        # Right-hand-side for the system of differential equations
        def rhs(t:float, x:NDArray) -> NDArray:
            print_progress(t, time[0], total_duration) # show progress

            # Get interpolated soc, current and temperature
            soc_t = soc_ip(t)
            current_t = current_ip(t)
            T_t = T_ip(t)

            # Load the statistical matrix
            stat_matrix = self.compile_statistical_matrix(T_t, soc_t) # shape = (5, 6)

            T_all = np.tensordot(stat_matrix, self.T_Phi, axes=([1],[0])) # shape = (5, basis_size, basis_size)
            C_all = np.tensordot(stat_matrix, self.E_Phi, axes=([1],[0])) # shape = (5, basis_size)

            H = block_diag(T_all[1], T_all[3]) # shape = (2*basis_size, 2*basis_size)
            b = np.concatenate([C_all[2], C_all[4]]) # shape = (2*basis_size,)

            return -np.sum(x * H, -1) + current_t * b
        

        coefficients = solve_ivp(
                                fun=rhs, max_step=self.max_solver_step, y0=np.zeros(2 * len(self.Phi)),
                                t_span=(time[0], time[-1]), t_eval=time
                                ).y.T
        

        coefficients1, coefficients2 = (
                                        coefficients[:, :len(self.Phi)], 
                                        coefficients[:, len(self.Phi):]
                                        )

        # Get the interpolator for R0 separately
        params = self.interpolator.get_interpolated_params(T, soc)
        mu_r0 = params['mu_R0']
        sigma_r0 = params['sigma_R0']

        vb_approx = (
                    (mu_r0 + self.nu_r * sigma_r0) * current  
                    + chaospy.sum(self.Phi * coefficients1, -1) 
                    + chaospy.sum(self.Phi * coefficients2, -1)
                    )
        
        # --- (FOR PACK SIMULATION) --- #
        self.output_pce = vb_approx  
        self.mu_r0 = mu_r0
        self.sigma_r0 = sigma_r0
        self.vb_transient = vb_approx - (mu_r0 + self.nu_r * sigma_r0) * current
        # ----------------------------- #

        mean = chaospy.E(vb_approx, self.joint)
        variance = np.abs(chaospy.Var(vb_approx, self.joint))
        sigma = np.sqrt(variance)

        sys.stdout.write('\n')
        sys.stdout.flush()
        self.logger.info(f'[{__class__.__name__}] Finished solving!')

        self.voltage.mean, self.voltage.std = mean, sigma # save mean and std
        return mean, sigma
    
    def get_output_expansion(self) -> Tuple[NDArray[chaospy.polynomial], chaospy.joint]:
        assert self.output_pce is not None, "No output available at this point: run solve() to get one!"

        return self.output_pce, self.joint
