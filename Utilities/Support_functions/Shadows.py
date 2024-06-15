import numpy as np
try:
    import Spacetimes
except:
    from Support_functions import Spacetimes

import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.stats import beta

class Wormhole_Shadow:

    def __init__(self, spin: float, gamma: float, obs_distance: float, obs_inclanation: float, Granularity: int) -> None:

        if spin == 0:
            print("Zero spin parameter breaks some of the expressions - setting it to 1e-10...")
            spin = 1e-10

        self.SPIN = spin
        self.GAMMA = gamma
        self.R_OBS = obs_distance
        
        if obs_inclanation == 0:
            print("Zero observer inclination breaks some of the expressions - setting it to 0.01...")
            obs_inclanation = 0.01

        self.THETA_OBS = obs_inclanation
        self.GRANULARITY = Granularity
        self.R_THROAT = 1

    def get_integrals_of_motion(self, r_ph: np.ndarray) -> tuple:

        omega   = 2 * self.SPIN / r_ph**3
        dr_omega = -3 * omega / r_ph

        N    = np.exp(-1 / r_ph - self.GAMMA / r_ph**2)
        dr_N = (1 / r_ph**2 + 2 * self.GAMMA / r_ph**3) * N
        
        Sigma = dr_N / N - 1 / r_ph

        ksi = Sigma / (Sigma * omega - dr_omega)
        eta = r_ph**2 / N**2 * (1 - omega * ksi)**2 

        return eta, ksi
    
    def get_branch_intersection_condition(self, r: float) -> float:

        eta, ksi = self.get_integrals_of_motion(r)

        _, _, _, _, N_throat, omega_throat, _ = self.get_metric(self.R_THROAT, self.THETA_OBS)

        #--- The above calculates the integrals of motion, using the parametric equations for photon orbits outside the throat, then imposes 
        #--- the condition they must satisfy ON the throat. This finds the intersection points between the two shadow branches (parametrized via r_ph)

        return eta * N_throat**2 / self.R_THROAT**2 - (1 - omega_throat * ksi)**2
    
    def get_branch_intersection_point(self):

        return fsolve(self.get_branch_intersection_condition, x0 = 10, xtol = 1e-6)

    def get_orbital_condition(self, r_ph) -> float:

        eta, ksi = self.get_integrals_of_motion(r_ph)
    
        return eta - ksi**2

    def get_retrograde_photon_orbit(self) -> float:

        return fsolve(self.get_orbital_condition, x0 = 10, xtol = 1e-3)
    
    def get_metric(self, r: float, theta: float):

        N = np.exp(-1 / r - self.GAMMA / r**2)
        omega = 2 * self.SPIN / r**3
        sin_theta = np.sin(theta)

        g_tt = -N**2 + omega**2 * r**2 * sin_theta**2

        try:
            g_rr = 1 / (1 - self.R_THROAT / r)

        except:
            g_rr = 0

        g_th_th = r**2
        g_ph_ph = r**2 * sin_theta**2
        g_t_ph  = -omega * r**2 * sin_theta**2

        return g_tt, g_th_th, g_ph_ph, g_t_ph, N, omega, g_rr
    
    def generate_shadow_rim(self, Eta, Ksi):

        # ========================== Metric at the observer ========================== #

        g_tt, g_th_th, g_ph_ph, g_t_ph, N, omega, g_rr = self.get_metric(self.R_OBS, self.THETA_OBS)

        g2 = g_t_ph**2 - g_tt*g_ph_ph

        inv_lapse  = np.sqrt(g_ph_ph / g2)
        gamma_coef = -g_t_ph / g_ph_ph * inv_lapse

        Theta_potential =  Eta - Ksi**2 / np.sin(self.THETA_OBS)**2
        Rad_potential   = -Eta * N**2 / self.R_OBS**2 + (1 - omega * Ksi)**2

        #-- Pick out only the array components that have a non-negative Theta potential.
        #-- The radial potential should never be able to go negative

        Rad_potential   = Rad_potential[Theta_potential >= 0]
        Ksi             = Ksi[Theta_potential >= 0]   
        Eta             = Eta[Theta_potential >= 0]    
        Theta_potential = Theta_potential[Theta_potential >= 0]    

        # ================ Local Contravariant Momenta at the observer ================ #

        p_t     = inv_lapse - gamma_coef * Ksi
        p_r     = np.sqrt(Rad_potential) / N / np.sqrt(g_rr)
        p_theta = np.sqrt(Theta_potential) / np.sqrt(g_th_th)
        p_phi   = Ksi / np.sqrt(g_ph_ph)

        # =========================== Celestial coordinates =========================== #

        alpha_coords = np.arctan( p_phi / p_r )
        beta_coords  = np.arcsin( p_theta / p_t )
        
        #-- Get rid of that weird arc that goes "backwards" and does not contribute to the shadow...

        if alpha_coords.size != 0:

            if self.SPIN > 0:
                index = np.argmax(alpha_coords)
            else:
                index = np.argmin(alpha_coords)

            alpha_coords = alpha_coords[index:]
            beta_coords = beta_coords[index:]

        return alpha_coords, beta_coords

    def generate_shadow(self, Subplot, plot_center_cross: bool = False, linestyle: str = "-", Plot_title: bool = False) -> tuple:

        # =========================================================================================================== #
        # ============================== Branch due to photon orbits outside the throat ============================= #
        # =========================================================================================================== #

        r_ph_retrograde = self.get_retrograde_photon_orbit()

        #-- The beta CDF maps more points in the uniform interval [0, 1] towards 1. Scaling this by the range of possible photon orbits
        #-- results in a distribution that is more "dense" towards higher radii - this "closes" the shadow on the "fat" side a lot better

        r_ph_scan = beta.cdf(x = np.linspace(0, 1, self.GRANULARITY), a = 1, b = 10) * (r_ph_retrograde - self.R_THROAT) + self.R_THROAT

        Eta, Ksi = self.get_integrals_of_motion(r_ph = r_ph_scan)

        alpha_first_coords, beta_first_coords = self.generate_shadow_rim(Eta, Ksi)

        # =========================================================================================================== #
        # ================================ Branch due to photon orbits on the throat ================================ #
        # =========================================================================================================== #

        #-- Pick out the value for the azimuthal impact parameter at which the two shadow branches meet (a.e the first part of allowed range for Ksi)

        _, Ksi_meet = self.get_integrals_of_motion(self.R_THROAT)

        # ============================== Metric at the throat ============================== #

        _, _, _, _, N_throat, omega_throat, _ = self.get_metric(self.R_THROAT, self.THETA_OBS)

        # ================================================================================== #

        Ksi_intersect = 1 / (N_throat / self.R_THROAT / np.sin(self.THETA_OBS) + omega_throat)

        Ksi = np.linspace(Ksi_intersect, Ksi_meet, self.GRANULARITY)
        Eta = self.R_THROAT**2 / N_throat**2 * (1 - omega_throat * Ksi)**2

        alpha_second_coords, beta_second_coords = self.generate_shadow_rim(Eta, Ksi)

        # =========================================================================================================== #
        # ======================================== Proper branch intersecting ======================================= #
        # =========================================================================================================== #

        eta_int, ksi_int = self.get_integrals_of_motion(self.get_branch_intersection_point())

        alpha_int, beta_int = self.generate_shadow_rim(eta_int, ksi_int)
    
        if beta_int.size > 0:

            if self.SPIN > 0:

                throat_bitmask  = beta_second_coords**2  < beta_int**2
                outside_bitmask = alpha_first_coords < alpha_int

            else:

                throat_bitmask  = beta_second_coords**2  < beta_int**2
                outside_bitmask = alpha_first_coords > alpha_int

            alpha_first_coords  = alpha_first_coords[outside_bitmask]
            beta_first_coords   = beta_first_coords[outside_bitmask]
            alpha_second_coords = alpha_second_coords[throat_bitmask]
            beta_second_coords  = beta_second_coords[throat_bitmask]

        if plot_center_cross:

            Subplot.plot([-100, 100], [0, 0], color = "k", linestyle = "--", linewidth = 1)
            Subplot.plot([0, 0], [-100, 100], color = "k", linestyle = "--", linewidth = 1)

        if len(beta_second_coords) > 50:

            X_plot_first  = np.concatenate([[min(alpha_second_coords)], alpha_first_coords, np.flip(alpha_first_coords), [min(alpha_second_coords)]]) * self.R_OBS
            X_plot_second = np.concatenate([np.flip(alpha_second_coords), alpha_second_coords]) * self.R_OBS
            
            Y_plot_first  = np.concatenate([[-max(beta_second_coords)], -beta_first_coords, np.flip(beta_first_coords), [max(beta_second_coords)]]) * self.R_OBS
            Y_plot_second = np.concatenate([-np.flip(beta_second_coords), beta_second_coords]) * self.R_OBS

            Subplot.plot(X_plot_first, Y_plot_first, color = "black", linestyle = linestyle)
            Subplot.plot(X_plot_second, Y_plot_second, color = "red", linestyle = linestyle)

        else:
            
            X_plot_first  = np.concatenate([alpha_first_coords, np.flip(alpha_first_coords), [alpha_first_coords[0]]]) * self.R_OBS
            Y_plot_first  = np.concatenate([-beta_first_coords, np.flip(beta_first_coords), [-beta_first_coords[0]]]) * self.R_OBS

            Subplot.plot(X_plot_first, Y_plot_first, color = "black", linestyle = linestyle)

        Subplot.set_xlabel(r"$r_{\text{obs}}\alpha$,  [M]", fontsize = 15)
        Subplot.set_ylabel(r"$r_{\text{obs}}\beta$,  [M]", fontsize = 15)

        if Plot_title:
            Subplot.set_title(r"$\gamma = {},\,\, a = {}$".format(self.GAMMA, round(self.SPIN,1)), fontsize = 15)

        # Subplot.set_xlim([-8, 8])
        # Subplot.set_ylim([-8, 8])
        # return alpha_coords, beta_coords


def add_Kerr_Shadow(spin: float, obs_distance: float, obs_inclanation: float, figure: plt) -> None:

    #--- Equatorial Photon Orbits ---#

    r_ph_max = 2*(1 + np.cos(2/3 * np.arccos(-spin)))
    r_ph_min = 2*(1 + np.cos(2/3 * np.arccos( spin)))

    r_ph = np.linspace(r_ph_min, r_ph_max, 10000)

    #--- Critical Integrals of Motion ---#

    J = (spin**2 * (r_ph + 1) + r_ph**3 - 3*r_ph**2)/(spin*(r_ph - 1))
    K = (-9*r_ph**4 + 4*spin**2*r_ph**3 + 6*r_ph**5 - r_ph**6)/(r_ph - 1)**2/spin**2

    cos_0 = np.cos(obs_inclanation)
    sin_0 = np.sin(obs_inclanation)
    tan_0 = np.tan(obs_inclanation)

    #--- Metric at the Observer ---#
    g_phph = (obs_distance**2 + spin**2 + 2*spin**2*sin_0**2/obs_distance)*sin_0**2
    g_tph  = -2*spin*obs_distance*sin_0**2/(obs_distance**2 + spin**2*cos_0**2)
    g_tt   = -(1 - 2*obs_distance/(obs_distance**2 + spin**2*cos_0**2))

    g2    = g_tph**2 - g_tt*g_phph
    gamma = -g_tph/g_phph*np.sqrt(g_phph/g2)
    ksi   = np.sqrt(g_phph/g2)

    rho2  = obs_distance**2 + spin**2*cos_0**2
    delta = obs_distance**2 - 2*obs_distance + spin**2 

    Rad_potential = obs_distance**4 + (spin**2 - K - J**2)*obs_distance**2 + 2*obs_distance*(K + (spin - J)**2) - K*spin**2
    
    #--- Celestial Coordinates of the Shadow Rim ---#

    alpha = np.arctan(J * np.sqrt(rho2 * delta/Rad_potential/g_phph))

    #--- This should be non-negative on the shadow rim ---#
    problem_term = (K + spin**2*cos_0**2 - J**2/tan_0**2) * (K + spin**2*cos_0**2 - J**2/tan_0**2 > 0)
    beta = np.arcsin(np.sqrt(problem_term)/np.sqrt(obs_distance**2 + spin**2*cos_0**2)/(ksi - J*gamma)) 

    alpha = obs_distance * alpha[(beta != 0)]
    beta  = obs_distance * beta [(beta != 0)]

    alpha = np.concatenate([alpha, np.flip(alpha), [alpha[0]]])
    beta  = np.concatenate([beta ,-np.flip(beta),  [beta[0]] ])

    figure.plot(alpha,beta,"k--")

if __name__ == "__main__":

    Fig = plt.figure(figsize = (18, 18))

    obs_inclination = np.pi / 2
    Granularity = 100000

    Subplot = Fig.add_subplot(331)
    Subplot.set_aspect(1)
    Shadow_1 = Wormhole_Shadow(0, 0.0, 1e5, obs_inclination, Granularity)
    Shadow_1.generate_shadow(Subplot)

    Subplot = Fig.add_subplot(332)
    Subplot.set_aspect(1)
    Shadow_2 = Wormhole_Shadow(0, 1.0, 1e5, obs_inclination, Granularity)
    Shadow_2.generate_shadow(Subplot)

    Subplot = Fig.add_subplot(333)
    Subplot.set_aspect(1)
    Shadow_3 = Wormhole_Shadow(0, 2.0, 1e5, obs_inclination, Granularity)
    Shadow_3.generate_shadow(Subplot)

    Subplot = Fig.add_subplot(334)
    Subplot.set_aspect(1)
    Shadow_4 = Wormhole_Shadow(0.5, 0.0, 1e5, obs_inclination, Granularity)
    Shadow_4.generate_shadow(Subplot)

    Subplot = Fig.add_subplot(335)
    Subplot.set_aspect(1)
    Shadow_5 = Wormhole_Shadow(0.5, 1.0, 1e5, obs_inclination, Granularity)
    Shadow_5.generate_shadow(Subplot)

    Subplot = Fig.add_subplot(336)
    Subplot.set_aspect(1)
    Shadow_6 = Wormhole_Shadow(0.5, 2.0, 1e5, obs_inclination, Granularity)
    Shadow_6.generate_shadow(Subplot)

    Subplot = Fig.add_subplot(337)
    Subplot.set_aspect(1)
    Shadow_7 = Wormhole_Shadow(1, 0.0, 1e5, obs_inclination, Granularity)
    Shadow_7.generate_shadow(Subplot)

    Subplot = Fig.add_subplot(338)
    Subplot.set_aspect(1)
    Shadow_8 = Wormhole_Shadow(1, 1.0, 1e5, obs_inclination, Granularity)
    Shadow_8.generate_shadow(Subplot)

    Subplot = Fig.add_subplot(339)
    Subplot.set_aspect(1)
    Shadow_9 = Wormhole_Shadow(1, 2.0, 1e5, obs_inclination, Granularity)
    Shadow_9.generate_shadow(Subplot)

    Fig.savefig("WH_Shadows", bbox_inches = 'tight')

    # plt.xlim([-10 / test_class.R_OBS, 10 / test_class.R_OBS])
    # plt.ylim([-10 / test_class.R_OBS, 10 / test_class.R_OBS])
    plt.show()
