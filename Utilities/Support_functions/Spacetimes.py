import numpy as np
from scipy.stats import beta
from scipy.optimize import root_scalar

class Schwarzschild:

    def __init__(self):
        self.MASS = 1
        self.HAS_PHOTON_SPHERE = True
        self.HAS_R_OF_B = False

    def metric(self, r, theta) -> np.array:

        g_tt = -(1 - 2 * self.MASS / r)
        g_rr = -1 / g_tt
        g_thth   = r ** 2
        g_phiphi = g_thth * np.sin(theta)**2

        if type(r) != np.ndarray:
            return np.array([g_tt, g_rr, g_thth, g_phiphi])
        else:
            return np.column_stack([g_tt, g_rr, g_thth, g_phiphi])

    def photon_sphere(self):

        return 3 * self.MASS

    def ISCO(self):

        return 6 * self.MASS

class Regular_Black_Hole:

    def __init__(self, parameter):
        self.MASS = 1
        self.PARAMETER = parameter
        self.HAS_PHOTON_SPHERE = True
        self.HAS_R_OF_B = False

    def metric(self, r, theta) -> np.array:

        r_eff = np.sqrt(r**2 + self.PARAMETER**2)

        g_tt = -(1 - 2 * self.MASS / r_eff)
        g_rr = -1 / g_rr
        g_thth   = r_eff ** 2
        g_phiphi = g_thth * np.sin(theta)**2

        if type(r) != np.ndarray:
            return np.array([g_tt, g_rr, g_thth, g_phiphi])
        else:
            return np.column_stack([g_tt, g_rr, g_thth, g_phiphi])
    def photon_sphere(self):

        return np.sqrt( (3 * self.MASS)**2 - self.PARAMETER**2 )

    def ISCO(self):

        return np.sqrt( (6 * self.MASS)**2 - (self.PARAMETER)**2 )

class Wormhole:

    def __init__(self, r_throat, parameter):
        self.MASS = 1
        self.R_THROAT = r_throat
        self.PARAMETER = parameter
        self.HAS_PHOTON_SPHERE = True
        self.HAS_R_OF_B = False

    def metric(self, r, theta, use_global_coords: bool = True) -> np.array:


        if use_global_coords:
                
            """
            The metric uses global coordinates - the range of ell is [-inf, + inf] 
            
            """

            ell = r

            r = np.sqrt(ell**2 + self.R_THROAT**2)

            g_tt = -np.exp(-2 * self.MASS/ r - 2 * self.PARAMETER * (self.MASS / r)**2)
            g_rr = 1 + self.R_THROAT / r
            g_thth   = r**2
            g_phiphi = g_thth * np.sin(theta)**2

        else:

            g_tt = -np.exp(-2 * self.MASS/ r - 2 * self.PARAMETER * (self.MASS / r)**2)
            g_rr = 1 / (1 - self.R_THROAT / r)
            g_thth   = r**2
            g_phiphi = g_thth * np.sin(theta)**2

        if type(r) != np.ndarray:
            return np.array([g_tt, g_rr, g_thth, g_phiphi])
        else:
            return np.column_stack([g_tt, g_rr, g_thth, g_phiphi])

    def photon_sphere(self, use_global_coords: bool = True):

        r_ph = self.MASS / 2 * (1 + np.sqrt(1 + 8 * self.PARAMETER))

        if use_global_coords:
            return np.sqrt(r_ph**2 - self.R_THROAT**2)
        
        else:
            return r_ph

    def ISCO(self):

        return 2 * self.MASS * (np.sqrt(4 / 9 * (6 * self.PARAMETER + 1) ) * np.cosh(1 / 3 * np.arccosh( (1 + 9*self.PARAMETER + 27 / 2 * self.PARAMETER**2)/pow(6 * self.PARAMETER + 1, 3 / 2) )) + 1 / 3)

class JNW_Naked_Singularity:

    def __init__(self, parameter):
        self.MASS = 1
        self.PARAMETER = parameter

        if self.PARAMETER < 0.5:
            self.HAS_PHOTON_SPHERE = False

        else:
            self.HAS_PHOTON_SPHERE = True

        self.HAS_R_OF_B = True

    def metric(self, r, theta) -> np.array:

        r_singularity = 2 * self.MASS / self.PARAMETER

        g_tt = -pow(1 - r_singularity / r, self.PARAMETER)

        if g_tt != 0:
            g_rr = - 1 / g_tt
        else:
            g_rr = np.inf
     
        g_thth   = r**2 * pow(1 - r_singularity / r, 1 - self.PARAMETER) 
        g_phiphi = g_thth * np.sin(theta)**2

        if type(r) != np.ndarray:
            return np.array([g_tt, g_rr, g_thth, g_phiphi])
        else:
            return np.column_stack([g_tt, g_rr, g_thth, g_phiphi])

    def photon_sphere(self):

        if self.PARAMETER > 0.5:

            return self.MASS / self.PARAMETER * (2 * self.PARAMETER + 1) 

        else:

            print("Photon sphere does not exist for values of gamma < 0.5!")

            return None

    def ISCO(self):

        if (self.PARAMETER > 1 / np.sqrt(5)):

            r_ISCO_outer = 1 / self.PARAMETER * (3 * self.PARAMETER + 1 + np.sqrt(5 * self.PARAMETER**2 - 1))
            r_ISCO_inner = 1 / self.PARAMETER * (3 * self.PARAMETER + 1 - np.sqrt(5 * self.PARAMETER**2 - 1))

            return np.array([r_ISCO_outer, r_ISCO_inner])
        
        else:

            return 2 * self.MASS / self.PARAMETER
        
    def get_turning_point_equation(self, r_turning, impact_param):
                
        alpha = 1 / (1 / 2 - self.PARAMETER)

        return pow(r_turning, alpha) - 2 / self.PARAMETER * pow(r_turning, alpha - 1) - pow(impact_param, alpha)

    def get_impact_params_and_turning_points(self, GRANULARITY, r_source):
                
        density_parameter  = 5.5
        distribution_range = np.linspace(0, 1, GRANULARITY)

        if self.HAS_PHOTON_SPHERE:

            # In the presence of a photon sphere, the turning points will be between it and the source

            higher_order_turning_points = r_source + (beta.cdf(distribution_range, density_parameter, density_parameter)) * (self.photon_sphere() - r_source)

            metrics_at_turning_points   = self.metric(higher_order_turning_points, np.pi / 2)

            # Calculate the impact parameters for the correspoinding turning points

            test = -metrics_at_turning_points.T[3] / metrics_at_turning_points.T[0]

            higher_order_impact_params  = np.sqrt(-metrics_at_turning_points.T[3] / metrics_at_turning_points.T[0])

        else:

            # In the absence of a photon sphere the turning points are located between the singularity and the source

            r_singularity = 2 * self.MASS / self.PARAMETER
            b_source      = r_source * pow(1 - r_singularity / r_source, 1 / 2 - self.PARAMETER)

            # NOTE: To properly construct the images we need a very uneven distribution of photon impact parameters
            # The majority of the points must be distributed at the ends of the interval - this is neatly done with a beta distribution

            higher_order_impact_params = b_source * (1 - beta.cdf(distribution_range, density_parameter, density_parameter))

                    #-------- Calculate the turning points for the correspoinding impact parameters --------#

            higher_order_turning_points = []

            for impact_param in higher_order_impact_params:

                roots = root_scalar(self.get_turning_point_equation, 
                                            x0      = r_source,
                                            args    = (impact_param), 
                                            maxiter = 100000, 
                                            xtol    = 1e-28, 
                                            method  = 'toms748', 
                                            bracket = [r_singularity - 1, r_source + 1])

                if np.absolute(roots.root - r_singularity) > 1e-10:
                    higher_order_turning_points.append(roots.root)

                else:
                    higher_order_turning_points.append(r_singularity)

        return higher_order_impact_params, higher_order_turning_points
        
class Gaus_Bonnet_Naked_Singularity:

    def __init__(self, parameter):

        self.MASS = 1
        self.PARAMETER = parameter

        if self.PARAMETER < 1.35:

            self.HAS_PHOTON_SPHERE = True

        else:

            self.HAS_PHOTON_SPHERE = False

        self.HAS_R_OF_B = True

    def metric(self, r, theta) -> np.array:

        if self.PARAMETER != 0:
            f = 1 + r**2 / 2 / self.PARAMETER * (1 - np.sqrt(1 + 8 * self.PARAMETER * self.MASS / r**3))
        else:
            f = 1 - 2 * self.MASS / r

        g_tt = -f

        if g_tt != 0:
            g_rr = - 1 / g_tt
        else:
            g_rr = np.inf

        g_thth = r**2
        g_phiphi = g_thth * np.sin(theta)**2

        if type(r) != np.ndarray:
            return np.array([g_tt, g_rr, g_thth, g_phiphi])
        else:
            return np.column_stack([g_tt, g_rr, g_thth, g_phiphi])
    
    def photon_sphere(self):

        polynomial_coefs = [1, 0, -9 * self.MASS**2, 8 * self.MASS * self.PARAMETER]

        roots = np.roots(polynomial_coefs)
        roots = roots[roots > 0]

        return np.max(roots)
    
    def get_impact_params_and_turning_points(self, GRANULARITY: int, r_source: float) -> tuple:

        density_parameter  = 5.5

        metric_source       = self.metric(r_source, np.pi / 2)
        impact_param_source = np.sqrt(-metric_source[3] / metric_source[0])

        if self.HAS_PHOTON_SPHERE:
            distribution_range_1 = np.linspace(0, 1, int(GRANULARITY / 2))
            distribution_range_2 = np.linspace(0, 1, GRANULARITY - int(GRANULARITY / 2))

            metric_photon_sphere       = self.metric(self.photon_sphere(), np.pi / 2) 
            impact_param_photon_sphere = np.sqrt(-metric_photon_sphere[3] / metric_photon_sphere[0])

            higher_order_impact_params1 = impact_param_source        + beta.cdf(distribution_range_1, density_parameter, density_parameter) * (impact_param_photon_sphere - impact_param_source)
            higher_order_impact_params2 = impact_param_photon_sphere + beta.cdf(distribution_range_2, density_parameter, density_parameter) * (1e-2 - impact_param_photon_sphere)
                                                                                                                                                
            higher_order_impact_params = np.append(higher_order_impact_params1, higher_order_impact_params2)

        else:
            distribution_range = np.linspace(0, 1, GRANULARITY)
            higher_order_impact_params = impact_param_source + beta.cdf(distribution_range, density_parameter, density_parameter) * (1e-2 - impact_param_source)
        
        #-------- Calculate the impact parameters for the correspoinding turning points --------#
        
        higher_order_turning_points = []
        higher_order_impact_params_to_return = []

        for index, impact_param in enumerate(higher_order_impact_params):

            a = (impact_param**2 / 2 / self.PARAMETER - 1)**2 - impact_param**4 / 4 / self.PARAMETER**2
            b = 0
            c = 2 * impact_param**2 * (impact_param**2 / 2 / self.PARAMETER - 1)
            d = -2 * impact_param**4 / self.PARAMETER
            e = impact_param**4

            polynom_coeffs = [a, b, c, d, e]
                              
            roots = np.roots(polynom_coeffs)
            roots = roots[np.abs(np.imag(roots)) < 1e-14]
            roots = np.real(roots)
            roots = roots[roots > 0]
            # roots = roots[roots < self.photon_sphere()]

            if self.PARAMETER <= 1 and max(roots) < self.photon_sphere():
                continue

            if len(roots) == 0:
                higher_order_turning_points.append(1e-10)      
            elif len(roots) != 2:
                higher_order_turning_points.append(np.max(roots))
            else:
                higher_order_turning_points.append(np.min(roots))

            higher_order_impact_params_to_return.append(impact_param)

        return higher_order_impact_params_to_return, higher_order_turning_points
