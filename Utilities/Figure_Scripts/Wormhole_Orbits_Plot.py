import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq, fsolve

class Wormhole_Orbits:

    def __init__(self, spin: float, gamma: float) -> None:

        if spin == 0:
            print("Zero spin parameter breaks some of the expressions - setting it to 1e-10...")
            spin = 1e-10
        
        self.SPIN     = spin
        self.GAMMA    = gamma
        self.R_THROAT = 1

    def get_metric(self, r: float, theta: float):

        N = np.exp(-1 / r - self.GAMMA / r**2)
        omega = 2 * self.SPIN / r**3
        sin_theta = np.sin(theta)

        g_tt = -N**2 + omega**2 * r**2 * sin_theta**2
        g_rr = 1 / (1 - self.R_THROAT / r)
        g_th_th = r**2
        g_ph_ph = r**2 * sin_theta**2
        g_t_ph  = -omega * r**2 * sin_theta**2

        return g_tt, g_rr, g_th_th, g_ph_ph, g_t_ph, N, omega
    
    def get_dr_metric(self, r: float, theta: float):

        _, _, _, _, _, N, omega = self.get_metric(r = r, theta = theta)

        dr_N = N * (1 / r**2 + 2 * self.GAMMA / r**3)
        dr_omega = -3 * omega / r

        dr_g_tt = -2 * N * dr_N + 2 * r * omega * (omega + r * dr_omega) * np.sin(theta)**2
        dr_g_t_phi = -r * (2 * omega + r * dr_omega) * np.sin(theta)**2
        dr_g_rr = - 1 / (1 - self.R_THROAT / r)**2 * (self.R_THROAT / r**2)
        dr_g_th_th = 2 * r
        dr_g_ph_ph = dr_g_th_th * np.sin(theta)**2

        return dr_g_tt, dr_g_rr, dr_g_th_th, dr_g_ph_ph, dr_g_t_phi, dr_N, dr_omega
    
    def get_d2r_metric(self, r: float, theta: float):

        _, _, _, _, _, N, omega = self.get_metric(r = r, theta = theta)
        _, _, _, _, _, dr_N, dr_omega = self.get_dr_metric(r = r, theta = theta)

        d2r_N = dr_N * (1 / r**2 + 2 * self.GAMMA / r**3) - N * (2 / r**3 + 6 * self.GAMMA / r**4)
        d2r_omega = -3 * dr_omega / r + 3 * omega / r**2 

        d2r_g_tt = -2 * dr_N**2 - 2 * N * d2r_N + 2 * ((omega + r * dr_omega)**2 + r * omega * (2 * dr_omega + r * d2r_omega)) * np.sin(theta)**2
        d2r_g_t_phi = -(2 * omega + r * dr_omega + r * (3 * dr_omega + r * d2r_omega)) * np.sin(theta)**2
        d2r_g_rr = 2 / (1 - self.R_THROAT / r)**2 * ((self.R_THROAT / r**2)**2 / (1 - self.R_THROAT / r) + self.R_THROAT / r**3)
        d2r_g_th_th = 2
        d2r_g_ph_ph = d2r_g_th_th * np.sin(theta)**2

        return d2r_g_tt, d2r_g_rr, d2r_g_th_th, d2r_g_ph_ph, d2r_g_t_phi, d2r_N, d2r_omega

    def get_photon_orbit_condition(self, r: float):

        N = np.exp(-1 / r - self.GAMMA / r**2)
        omega = 2 * self.SPIN / r**3
        dr_omega = -3 * omega / r
        Sigma = - 1 / r + 1 / r**2 + 2 * self.GAMMA / r**3
        Ksi = Sigma / (Sigma * omega - dr_omega)
        Eta = r**2 / N**2 * (1 - omega * Ksi)**2

        return Eta - Ksi**2
    
    def get_Kepler_velocity(self, r: float, direction: int) -> float:

        dr_g_tt,  _, _, dr_g_ph_ph,  dr_g_t_phi,  _, _ = self.get_dr_metric(r, theta = np.pi / 2)

        return -dr_g_t_phi / dr_g_ph_ph + direction * np.sqrt(dr_g_t_phi**2 - dr_g_tt * dr_g_ph_ph) / dr_g_ph_ph
    
    def get_Energy(self, r: float, direction: int) -> float:

        g_tt, _, _, g_ph_ph, g_t_ph, _, _ = self.get_metric(r, theta = np.pi / 2)

        Kepler_vel = self.get_Kepler_velocity(r, direction)

        return -(g_tt + g_t_ph * Kepler_vel) / np.sqrt(-g_tt - 2 * g_t_ph * Kepler_vel - g_ph_ph * Kepler_vel**2)

    def get_mb_orbit_condition(self, r: float , direction: int) -> float:

        if np.isnan(self.get_Energy(r, direction)):

            return  500
        
        return self.get_Energy(r, direction) - 1
    
    def get_ISCO_condition(self, r: float, direction: int):

        g_tt,     _, _, g_ph_ph,     g_t_ph,      _, _ = self.get_metric(r, theta = np.pi / 2)
        dr_g_tt,  _, _, dr_g_ph_ph,  dr_g_t_phi,  _, _ = self.get_dr_metric(r, theta = np.pi / 2)
        d2r_g_tt, _, _, d2r_g_ph_ph, d2r_g_t_phi, _, _ = self.get_d2r_metric(r, theta = np.pi / 2)

        Kepler_vel = self.get_Kepler_velocity(r, direction)

        E_norm = -g_tt - 2 * g_t_ph * Kepler_vel - g_ph_ph * Kepler_vel**2

        if E_norm > 0:

            E   = self.get_Energy(r, direction)
            L_z =  (g_t_ph + g_ph_ph * Kepler_vel) / np.sqrt(E_norm)
            ell = L_z / E

            d2r_g2 = 2 * dr_g_t_phi**2 + 2 * g_t_ph * d2r_g_t_phi - d2r_g_ph_ph * g_tt - 2 * dr_g_ph_ph * dr_g_tt - g_ph_ph * d2r_g_tt; 

            d2r_V_eff = -d2r_g2 + E**2 * (d2r_g_ph_ph + 2 * ell * d2r_g_t_phi + ell**2 * d2r_g_tt)

            return d2r_V_eff
    
        else:

            return 500

    def analytic_ISCO(self):
        return 2 * (np.sqrt(4 / 9 * (6 * self.GAMMA + 1)) * np.cosh( 1 / 3 * np.arccosh((1 + 9 * self.GAMMA + 27/2 * self.GAMMA**2) / (6 * self.GAMMA + 1 )**(3/2)) )  + 1 / 3 )

    def analytic_r_ph(self):
        return (1 + np.sqrt(1 + 8 * self.GAMMA)) / 2


cmap = matplotlib.colormaps['plasma']

spin_range = np.linspace(1e-10, 2, 1000)
gamma_range = np.linspace(1e-10, 3, 11)

MAGIC_NUMBER = 1.17 / 2

Photon_orbit_figure = plt.figure().add_subplot(111)

MB_orbit_figure = plt.figure()
MB_prograde_figure = MB_orbit_figure.add_subplot(121)
MB_retrograde_figure = MB_orbit_figure.add_subplot(122)

ISCO_figure = plt.figure()
ISCO_prograde_figure = ISCO_figure.add_subplot(121)
ISCO_retrograde_figure = ISCO_figure.add_subplot(122)

for gamma in gamma_range:

    r_ISCO_prograde = []
    r_ISCO_retrograde = []
    r_ph = []
    r_mb_retrograde = []
    r_mb_prograde = []

    for spin in spin_range:

        WH_Orbits_Class = Wormhole_Orbits(spin, gamma)

        r_ph.append(fsolve(WH_Orbits_Class.get_photon_orbit_condition, x0 = 10))

        r_ISCO_retrograde.append(brentq(WH_Orbits_Class.get_ISCO_condition, args = -1, a = 1.5 * WH_Orbits_Class.R_THROAT , b = 15, maxiter = 100000))
        r_mb_retrograde.append(brentq(WH_Orbits_Class.get_mb_orbit_condition, args = -1, a = 1.2 , b = 15, maxiter = 100000))

        try:
            r_ISCO_prograde.append(brentq(WH_Orbits_Class.get_ISCO_condition, args = 1, a = MAGIC_NUMBER * WH_Orbits_Class.analytic_ISCO(), b = 8, maxiter = 100000))
        except:
            True  

        try:
            r_mb_prograde.append(brentq(WH_Orbits_Class.get_mb_orbit_condition, args = 1, a = r_ISCO_prograde[-1] / 2 , b = 8, maxiter = 100000))
        except:
            True 

    dx    = spin_range[int(len(spin_range) / 2)] - spin_range[int(len(spin_range) / 2 - 1)]
    dy    = r_ph[int(len(spin_range) / 2)] - r_ph[int(len(spin_range) / 2 - 1)]
    angle = np.rad2deg(np.arctan2(dy, dx))[0]

    Photon_orbit_figure.plot(spin_range, r_ph, color = cmap(gamma / (gamma_range[-1] + 0.2)))
    Photon_orbit_figure.text(spin_range[int(len(spin_range) / 2)], r_ph[int(len(spin_range) / 2)], r'$\gamma=${}'.format(round(gamma, 2)), ha='left', va='bottom',
                             transform_rotates_text=True, rotation=angle, rotation_mode='anchor')

    dx    = spin_range[int(len(spin_range) / 2)] - spin_range[int(len(spin_range) / 2 - 1)]
    dy    = r_ISCO_retrograde[int(len(spin_range) / 2)] - r_ISCO_retrograde[int(len(spin_range) / 2 - 1)]
    angle = np.rad2deg(np.arctan2(dy, dx))

    ISCO_retrograde_figure.plot(spin_range, r_ISCO_retrograde, color = cmap(gamma / (gamma_range[-1] + 0.2)))
    ISCO_retrograde_figure.text(spin_range[int(len(spin_range) / 2)], r_ISCO_retrograde[int(len(spin_range) / 2)], r'$\gamma=${}'.format(round(gamma, 2)), ha='left', va='bottom',
                                transform_rotates_text=True, rotation=angle, rotation_mode='anchor')

    dx    = spin_range[int(len(spin_range) / 2)] - spin_range[int(len(spin_range) / 2 - 1)]
    dy    = r_mb_retrograde[int(len(spin_range) / 2)] - r_mb_retrograde[int(len(spin_range) / 2 - 1)]
    angle = np.rad2deg(np.arctan2(dy, dx))

    MB_retrograde_figure.plot(spin_range[:len(r_mb_retrograde)], r_mb_retrograde, color = cmap(gamma / (gamma_range[-1] + 0.2)))
    MB_retrograde_figure.text(spin_range[int(len(spin_range) / 2)], r_mb_retrograde[int(len(spin_range) / 2)], r'$\gamma=${}'.format(round(gamma, 2)), ha='left', va='bottom',
                              transform_rotates_text=True, rotation=angle, rotation_mode='anchor')

    try:
        ISCO_prograde_figure.plot(spin_range[:len(r_ISCO_prograde)], r_ISCO_prograde, color = cmap(gamma / (gamma_range[-1] + 0.2)))

        dx    = spin_range[int(len(r_ISCO_prograde) / 3)] - spin_range[int(len(r_ISCO_prograde) / 3 - 1)]
        dy    = r_ISCO_prograde[int(len(r_ISCO_prograde) / 3)] - r_ISCO_prograde[int(len(r_ISCO_prograde) / 3 - 1)]
        angle = np.rad2deg(np.arctan2(dy, dx))

        ISCO_prograde_figure.text(spin_range[int(len(r_ISCO_prograde) / 3)], r_ISCO_prograde[int(len(r_ISCO_prograde) / 3)], r'$\gamma=${}'.format(round(gamma, 2)), ha='left', va='bottom',
                                  transform_rotates_text=True, rotation=angle, rotation_mode='anchor')

    except:
        True

    try:
        MB_prograde_figure.plot(spin_range[:len(r_mb_prograde)], r_mb_prograde, color = cmap(gamma / (gamma_range[-1] + 0.2)))

        dx    = spin_range[int(len(r_mb_prograde) / 3)] - spin_range[int(len(r_mb_prograde) / 3 - 1)]
        dy    = r_mb_prograde[int(len(r_mb_prograde) / 3)] - r_mb_prograde[int(len(r_mb_prograde) / 3 - 1)]
        angle = np.rad2deg(np.arctan2(dy, dx))

        MB_prograde_figure.text(spin_range[int(len(r_mb_prograde) / 3)], r_mb_prograde[int(len(r_mb_prograde) / 3)], r'$\gamma=${}'.format(round(gamma, 2)), ha='left', va='bottom',
                                transform_rotates_text=True, rotation=angle, rotation_mode='anchor')

    except:
        True

ISCO_prograde_figure.set_xlabel(r"$a\,\,[-]$", fontsize = 16)
ISCO_prograde_figure.set_ylabel(r"$r_{\text{ISCO}}^+\,\,[M]$", fontsize = 16)
ISCO_prograde_figure.set_xlim([spin_range[0], 0.42])

ISCO_retrograde_figure.set_xlabel(r"$a\,\,[-]$", fontsize = 16)
ISCO_retrograde_figure.set_ylabel(r"$r_{\text{ISCO}}^-\,\,[M]$", fontsize = 16)
ISCO_retrograde_figure.set_xlim([spin_range[0], spin_range[-1]])

MB_retrograde_figure.set_xlabel(r"$a\,\,[-]$", fontsize = 16)
MB_retrograde_figure.set_ylabel(r"$r_{\text{mb}}^-\,\,[M]$", fontsize = 16)
MB_retrograde_figure.set_xlim([spin_range[0], spin_range[-1]])

MB_prograde_figure.set_xlabel(r"$a\,\,[-]$", fontsize = 16)
MB_prograde_figure.set_ylabel(r"$r_{\text{mb}}^+\,\,[M]$", fontsize = 16)
MB_prograde_figure.set_xlim([spin_range[0], 0.42])

Photon_orbit_figure.set_xlabel(r"$a\,\,[-]$", fontsize = 16)
Photon_orbit_figure.set_ylabel(r"$r_{\text{ph}}^-\,\,[M]$", fontsize = 16)
Photon_orbit_figure.set_xlim([spin_range[0], spin_range[-1]])

plt.show()