import sys
import chaospy
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from scipy.linalg import block_diag
from collections.abc import Iterable
from abc import abstractmethod
from utils import create_logger, print_progress
from ..base import Model
from gaussian_process.interpolators import BatteryParameterInterpolator


class ECM(Model):
    """
    Template for an ECM object. Must have a solve method, which must return the simulated voltage array and the standard deviation of its uncertainty. If the model is deterministic, return zeros for the std.
    
    """

    def __init__(self):
        self.interpolator = None  # assigned by subclasses

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
    def solve(self, *args, **kwargs) -> Iterable[np.ndarray, np.ndarray]:
        pass


class GalerkinPCE(ECM):
    """
    Implements system matrices and a solver method for the 2RC Galerkin PCE.
    Immutable singleton model object for the cell class. Initializes only one GalerkinPCModel object to save up on memory; any subsequent call to the constructor returns a reference to this one object.

    """

    _instance = None

    def __new__(cls, *args, **kwargs):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self, with_correlations: bool = False, logger_level: str = 'DEBUG'):
        if not hasattr(self, '_initialized'):  # Prevent re-initialization

            required_vars = [
                'mu_r0', 'sigma_r0', # mean and std. R0
                'mu_tau1_inv', 'sigma_tau1_inv', # mean and std. 1/tau1
                'mu_c1_inv', 'sigma_c1_inv', # mean and std. 1/c1
                'mu_tau2_inv', 'sigma_tau2_inv', # mean and std. 1/tau2
                'mu_c2_inv', 'sigma_c2_inv', # mean and std. 1/c2
                'rho1', 'rho2', # correlations (1/c1, 1/tau1) and (1/c2, 1/tau2)
            ]

            super().__init__()
            self._initialized = True  # Flag to block re-init
            self.with_correlations = with_correlations
            self.logger = create_logger(__class__.__name__, logger_level)
            self.max_solver_step = 20. # set appropriate step size (20 is fine usually)
            self.interpolator = BatteryParameterInterpolator()  # interpolator for ECM parameters
            self.validate_interpolator(required_vars) # are all necessary variables available?

            # ~~~~ GALERKIN PCE ~~~~~
            # create the joint distribution, basis and parameter variables
            joint = chaospy.J(chaospy.Normal(), chaospy.Normal(), chaospy.Normal(),
                              chaospy.Normal(), chaospy.Normal())
            self.joint = joint

            self.expansion, self.norms = chaospy.generate_expansion(3, joint, retall=True)
            self.nu0, nu1, nu2, nu3, nu4 = chaospy.variable(5)

            # predefine necessary tensors
            phi_phi = chaospy.outer(self.expansion, self.expansion)
            phi = self.expansion

            # mean/std tensors for 1/RC
            self.e_phi_phi = chaospy.E(phi_phi, joint)  # H_mean
            self.e_nu_phi_phi1 = chaospy.E(nu1 * phi_phi, joint)  # H_sigma1
            self.e_nu_phi_phi3 = chaospy.E(nu3 * phi_phi, joint)  # H_sigma2

            # mean/std tensors for 1/C
            self.e_phi = chaospy.E(phi, joint)  # b_mean

            self.e_nu_phi1 = chaospy.E(nu1 * phi, joint)  # for rho1
            self.e_nu_phi2 = chaospy.E(nu2 * phi, joint)  # b_sigma1
            self.e_nu_phi3 = chaospy.E(nu3 * phi, joint)  # for rho2
            self.e_nu_phi4 = chaospy.E(nu4 * phi, joint)  # b_sigma2

            self.logger.info('[GalerkinPCE] Initialized!')

    def solve(self, **kwargs) -> Iterable[np.ndarray, np.ndarray]:
        """
        Parameters
        ----------
            **kwargs : dict of np.ndarray
                Input arrays. Keys: 'current', 'time', 'soc', 'T'.
        """
        
        current = kwargs['current']
        time = kwargs['time'] 
        soc = kwargs['soc']
        T = kwargs['T']

        self.logger.info('[GalerkinPCE] Starting solver...')

        # Set up interpolators for input arrays; required for the solver
        current_ip = interp1d(time, current, kind='linear', bounds_error=False, fill_value=(current[0], current[-1]))
        soc_ip = interp1d(time, soc, kind='linear', bounds_error=False, fill_value=(soc[0], soc[-1]))
        T_ip = interp1d(time, T, kind='linear', bounds_error=False, fill_value=(T[0], T[-1]))

        total_duration = time[-1]-time[0]
        # Right-hand-side for the system of differential equations
        def rhs(t:float, x:np.ndarray) -> np.ndarray:
            print_progress(t, time[0], total_duration) # show progress

            # Get interpolated soc, current and temperature
            soc_t = soc_ip(t)
            current_t = current_ip(t)
            T_t = T_ip(t)

            # Load and parse interpolated parameters
            params_t = self.interpolator.get_interpolated_params(soc_t, T_t)
            (
             mu_tau1_inv, sigma_tau1_inv, 
             mu_c1_inv, sigma_c1_inv, 
             mu_tau2_inv, sigma_tau2_inv, 
             mu_c2_inv, sigma_c2_inv,
             rho1, rho2,
            ) = (
             params_t['mu_tau1_inv'], params_t['sigma_tau1_inv'],
             params_t['mu_c1_inv'], params_t['sigma_c1_inv'],
             params_t['mu_tau2_inv'], params_t['sigma_tau2_inv'],
             params_t['mu_c2_inv'], params_t['sigma_c2_inv'],
             params_t['rho1'], params_t['rho2'],
            )

            if not self.with_correlations:
                rho1, rho2 = 0, 0

            #Calculate interpolated tensors:
            # H_mean*param_mean + H_sigma*param_sigma or b_mean*param_mean + b_sigma*param_sigma
            # apply additional Cholesky transformation to inject correlations
            e_tau1_phi_phi = mu_tau1_inv * self.e_phi_phi + sigma_tau1_inv * self.e_nu_phi_phi1
            e_c1_phi_phi = mu_c1_inv * self.e_phi + sigma_c1_inv*rho1*self.e_nu_phi1 + sigma_c1_inv*np.sqrt(1-rho1**2)*self.e_nu_phi2
            e_tau2_phi_phi = mu_tau2_inv * self.e_phi_phi + sigma_tau2_inv * self.e_nu_phi_phi3
            e_c2_phi_phi = mu_c2_inv * self.e_phi + sigma_c2_inv*rho2*self.e_nu_phi3 + sigma_c2_inv*np.sqrt(1-rho2**2)*self.e_nu_phi4

            # Combine separate RC equations into larger tensors (not necessary as of now, but can be useful for stabilizing transforms later)
            H = block_diag(e_tau1_phi_phi, e_tau2_phi_phi)
            b = np.concatenate([e_c1_phi_phi, e_c2_phi_phi])

            return (-np.sum(x * H, -1) + current_t * b) / np.concatenate([self.norms, self.norms])
        

        coefficients = solve_ivp(
                                fun=rhs, max_step=self.max_solver_step, y0=np.zeros(2 * len(self.expansion)),
                                t_span=(time[0], time[-1]), t_eval=time
                                ).y.T
        
        
        coefficients1, coefficients2 = (
                                        coefficients[:, :len(self.expansion)], 
                                        coefficients[:, len(self.expansion):]
                                        )

        # Get the interpolator for R0 separately
        mu_r0 =self.interpolator.interpolators[0]
        sigma_r0 = self.interpolator.interpolators[1]
        vb_approx = (
                    (mu_r0(T, soc) + self.nu0 * sigma_r0(T, soc)) * current + 
                    chaospy.sum(self.expansion * coefficients1, -1) +
                    chaospy.sum(self.expansion * coefficients2, -1)
                    )

        mean = chaospy.E(vb_approx, self.joint)
        variance = np.abs(chaospy.Var(vb_approx, self.joint))
        sigma = np.sqrt(variance)

        sys.stdout.write('\n')
        sys.stdout.flush()
        self.logger.info('[GalerkinPCE] Finished solving!')

        return mean, sigma


class nRC(ECM):
    """Standard implementation of an nRC ECM."""

    def __init__(self, N: int = 2, logger_level: str = 'DEBUG'):
        
        # define required variables
        required_vars = ['R0'] + [var for n in range(1, N+1) for var in (f'tau{n}_inv', f'c{n}_inv')]

        super().__init__()
        self.logger = create_logger(__class__.__name__, logger_level)
        self.max_solver_step = 20.
        self.interpolator = BatteryParameterInterpolator()
        self.validate_interpolator(required_vars)

        self.logger.info(f'[{__class__.__name__}] Initialized!')


    def solve(self, **kwargs) -> Iterable[np.ndarray, np.ndarray]:

        current = kwargs['current']
        time = kwargs['time'] 
        soc = kwargs['soc']
        T = kwargs['T']

        self.logger.info(f'[{__class__.__name__}] Starting solver...')

        # interpolators for input arrays
        current_ip = interp1d(time, current, kind='linear', bounds_error=False, fill_value=(current[0], current[-1]))
        soc_ip = interp1d(time, soc, kind='linear', bounds_error=False, fill_value=(soc[0], soc[-1]))
        T_ip = interp1d(time, T, kind='linear', bounds_error=False, fill_value=(T[0], T[-1]))

        total_duration = time[-1]-time[0]
        # Right-hand-side for the system of differential equations
        def rhs(t:float, x:np.ndarray) -> np.ndarray:
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
        
        return voltage, np.zeros_like(voltage)