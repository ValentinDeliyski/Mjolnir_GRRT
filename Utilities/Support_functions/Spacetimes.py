import numpy as np

class Schwarzschild:

    def __init__(self):
        self.MASS = 1
        self.HAS_PHOTON_SPHERE = True

    def metric(self, r, theta) -> np.array:

        g_tt = -(1 - 2 * self.MASS / r)
        g_rr = -1 / g_tt
        g_thth   = r ** 2
        g_phiphi = g_thth * np.sin(theta)**2

        return np.array([g_tt, g_rr, g_thth, g_phiphi])

    def photon_sphere(self):

        return 3 * self.MASS

    def ISCO(self):

        return 6 * self.MASS

class Regular_Black_Hole:

    def __init__(self, parameter):
        self.MASS = 1
        self.PARAMETER = parameter
        self.HAS_PHOTON_SPHERE = True

    def metric(self, r, theta) -> np.array:

        r_eff = np.sqrt(r**2 + self.PARAMETER**2)

        g_tt = -(1 - 2 * self.MASS / r_eff)
        g_rr = -1 / g_rr
        g_thth   = r_eff ** 2
        g_phiphi = g_thth * np.sin(theta)**2

        return np.array([g_tt, g_rr, g_thth, g_phiphi])

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

    def metric(self, r_global, theta) -> np.array:

        #------------------------------------------------------------------------#
        #                                                                        #
        #   Metric uses global coordinates - the range of ell is [-inf, + inf]   #
        #                                                                        #
        #------------------------------------------------------------------------#

        r = np.sqrt(r_global**2 + self.R_THROAT**2)

        g_tt = -np.exp(-2 * self.MASS/ r - 2 * self.PARAMETER * (self.MASS / r)**2)
        g_rr = 1 + self.R_THROAT / r
        g_thth   = r_global**2 + self.R_THROAT**2
        g_phiphi = g_thth * np.sin(theta)**2

        return np.array([g_tt, g_rr, g_thth, g_phiphi])

    def photon_sphere(self):

        r_ph = self.MASS / 2 * (1 + np.sqrt(1 + 8 * self.PARAMETER))

        return np.sqrt(r_ph**2 - self.R_THROAT**2)

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

    def metric(self, r, theta) -> np.array:

        r_singularity = 2 * self.MASS / self.PARAMETER

        g_tt = -pow(1 - r_singularity / r, self.PARAMETER)
        g_rr = - 1 / g_tt
        g_thth   = r**2 * pow(1 - r_singularity / r, 1 - self.PARAMETER) 
        g_phiphi = g_thth * np.sin(theta)**2

        return np.array([g_tt, g_rr, g_thth, g_phiphi])

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