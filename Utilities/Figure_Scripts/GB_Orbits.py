import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from scipy.optimize import brentq, root_scalar
import csv

def get_GB_f(r, gamma):

    return 1 + r**2 / 2 / gamma * (1 - np.sqrt(1 + 8 * gamma / r**3))

def get_GB_dr_f(r, gamma):

    return 2 / r * (get_GB_f(r, gamma) - 1) + 6 / r**2 / np.sqrt(1 + 8 * gamma / r**3)

def get_GB_d2r_f(r, gamma):

    return 2 / r**2 * (get_GB_f(r, gamma) - 1) + 12 / r / np.sqrt(r**4 + 8 * gamma * r) * (1 - (r**3 + 2 * gamma) / (r**3 + 8 * gamma))

def get_ISCO_Condition(r, gamma):

    g_tt = -get_GB_f(r, gamma)
    dr_g_tt = -get_GB_dr_f(r, gamma)
    d2r_g_tt = -get_GB_d2r_f(r, gamma)

    g_phph = r**2
    dr_g_phph = 2 * r
    d2r_g_phph = 2

    d2r_g_2 = -d2r_g_tt * g_phph - 2 * dr_g_phph * dr_g_tt - g_tt * d2r_g_phph

    root_function = -dr_g_tt / dr_g_phph
    Omega = np.sqrt(root_function) 


    ell = - g_phph / g_tt * Omega
    root_function = -g_tt - g_phph * Omega**2
    E = - g_tt / np.sqrt(root_function)

    d2r_V_eff = -d2r_g_2 + E**2 * (d2r_g_phph + ell**2 * d2r_g_tt)

    return d2r_V_eff

def get_mb_condition(r, gamma):

    g_tt = -get_GB_f(r, gamma)
    dr_g_tt = -get_GB_dr_f(r, gamma)

    g_phph = r**2
    dr_g_phph = 2 * r
  
    root_function = -dr_g_tt / dr_g_phph
    if root_function > 0:
        Omega = np.sqrt(root_function)
    else:
        return 0.01
    
    return g_tt**2 + g_tt + g_phph * Omega**2

def parse_file(file_name: str):

    with open(file_name, 'r') as file:

        csvreader = csv.reader(file, delimiter = ",")

        r_orbit = []
        gamma_value = []

        for row in csvreader:

            gamma_value.append(float(row[0]))
            r_orbit.append(float(row[1]))

    return gamma_value, r_orbit


params = {"ytick.color" : "black",
          "xtick.color" : "black",
          "axes.labelcolor" : "black",
          "axes.edgecolor" : "black",
          "text.usetex" : True,
          "font.family" : "serif",
          "font.serif" : ["Computer Modern Serif"]}
plt.rcParams.update(params)

gamma_range = np.linspace(0.01, 1.8, 3000)

r_mb = []
r_ph_outer = []

r_h_p = []
r_h_m = []
gamma_r_h = []

r_ISCO_outer = []

for gamma in gamma_range:

    coeffs_ph_orbits = [1, 0, -9, 8 * gamma]

    ph_orbits = np.roots(coeffs_ph_orbits)

    ph_orbits = ph_orbits[abs(np.imag(ph_orbits)) < 1e-10]
    ph_orbits = ph_orbits[ph_orbits > 0]

    if gamma < 1:
        ph_orbits = ph_orbits[ph_orbits > 1 + np.sqrt(1 - gamma)]

    if len(ph_orbits) > 0:
        r_ph_outer.append(max(ph_orbits))

    if gamma <= 1:
        r_h_p.append(1 + np.sqrt(1 - gamma))
        r_h_m.append(1 - np.sqrt(1 - gamma))
        gamma_r_h.append(gamma)

    r_ISCO_outer.append(root_scalar(get_ISCO_Condition, args = gamma, x0 = 4.5, xtol= 1e-12, method = "Newton", maxiter=1000).root)
    r_mb.append(brentq(get_mb_condition, args = gamma, a = 1, b = 10, xtol= 1e-10))

r_ISCO_outer = r_ISCO_outer[:np.where(np.array(r_ISCO_outer) < 6)[0][-1]]
r_mb         = r_mb[:np.where(np.abs(r_mb - gamma_range**(1/3)) < 0.1)[0][0]]

gamma_r_ISCO_inner, r_ISCO_inner = parse_file("Figure_Scripts\\GB_ISCO_Inner_Gamma.csv")
r_ISCO_inner.append(r_ISCO_outer[-10])
gamma_r_ISCO_inner.append(gamma_r_ISCO_inner[-1])

gamma_r_ph_inner, r_ph_inner = parse_file("Figure_Scripts\\GB_r_ph_Inner_Gamma.csv")

plt.tight_layout()

# ============================================================================================================= #
dx    = gamma_r_ph_inner[int(len(gamma_r_ph_inner) / 2)] - gamma_r_ph_inner[int(len(gamma_r_ph_inner) / 2 - 1)]
dy    = r_ph_inner[int(len(r_ph_inner) / 2)] - r_ph_inner[int(len(r_ph_inner) / 2 - 1)]

plt.text(gamma_r_ph_inner[int(len(r_ph_inner) / 2)], r_ph_inner[int(len(r_ph_inner) / 2)] + 0.05, 
         r'$r^-_{ph}$', 
         ha = 'left', 
         va = 'bottom',
         transform_rotates_text = True, 
         rotation = np.rad2deg(np.arctan2(dy, dx)), 
         rotation_mode = 'anchor', 
         fontsize = 16)

plt.plot(gamma_r_ph_inner, r_ph_inner, "r--")

# ============================================================================================================= #

dx    = gamma_range[int(len(gamma_range) / 2)] - gamma_range[int(len(gamma_range) / 2 - 1)]
dy    = r_ph_outer[int(len(r_ph_outer) / 2)] - r_ph_outer[int(len(r_ph_outer) / 2 - 1)]

plt.text(gamma_range[int(len(r_ph_outer) / 2)], r_ph_outer[int(len(r_ph_outer) / 2)] + 0.05, 
         r'$r^+_{ph}$', 
         ha = 'left', 
         va = 'bottom',
         transform_rotates_text = True, 
         rotation = np.rad2deg(np.arctan2(dy, dx)), 
         rotation_mode = 'anchor', 
         fontsize = 16)

plt.plot(gamma_range[:len(r_ph_outer)], r_ph_outer, "r-")

# ============================================================================================================= #

dx    = gamma_range[int(len(gamma_range) / 2)] - gamma_range[int(len(gamma_range) / 2 - 1)]
dy    = r_mb[int(len(r_mb) / 2)] - r_mb[int(len(r_mb) / 2 - 1)]

plt.text(gamma_range[int(len(r_mb) / 2)], r_mb[int(len(r_mb) / 2)] + 0.05, 
         r'$r_{mb}$', 
         ha = 'left', 
         va = 'bottom',
         transform_rotates_text = True, 
         rotation = np.rad2deg(np.arctan2(dy, dx)), 
         rotation_mode = 'anchor', 
         fontsize = 16)

plt.plot(gamma_range[:len(r_mb)], r_mb, "b:")

# ============================================================================================================= #

dx    = gamma_range[int(len(gamma_range) / 2)] - gamma_range[int(len(gamma_range) / 2 - 1)]
dy    = r_ISCO_outer[int(len(r_ISCO_outer) / 2)] - r_ISCO_outer[int(len(r_ISCO_outer) / 2 - 1)]

plt.text(gamma_range[int(len(r_ISCO_outer) / 2)], r_ISCO_outer[int(len(r_ISCO_outer) / 2)] + 0.05, 
         r'$r_{ISCO}^+$', 
         ha = 'left', 
         va = 'bottom',
         transform_rotates_text = True, 
         rotation = np.rad2deg(np.arctan2(dy, dx)), 
         rotation_mode = 'anchor', 
         fontsize = 16)

plt.plot(gamma_range[:len(r_ISCO_outer)], r_ISCO_outer, "b-")

# ============================================================================================================= #

dx    = gamma_r_ISCO_inner[int(len(gamma_r_ISCO_inner) / 1.5)] - gamma_r_ISCO_inner[int(len(gamma_r_ISCO_inner) / 1.5 - 1)]
dy    = r_ISCO_inner[int(len(r_ISCO_inner) / 1.5)] - r_ISCO_inner[int(len(r_ISCO_inner) / 1.5 - 1)]

plt.text(gamma_r_ISCO_inner[int(len(r_ISCO_inner) / 1.5)], r_ISCO_inner[int(len(r_ISCO_inner) / 1.5)] + 0.05, 
         r'$r_{ISCO}^-$', 
         ha = 'left', 
         va = 'bottom',
         transform_rotates_text = True, 
         rotation = np.rad2deg(np.arctan2(dy, dx)) * 1.05, 
         rotation_mode = 'anchor', 
         fontsize = 16)

plt.plot(gamma_r_ISCO_inner, r_ISCO_inner, "b--")

# ============================================================================================================= #

gamma_range = np.linspace(0.01, 2, 3000)
gamma_range = gamma_range[gamma_range > 1]

r_in = gamma_range**(1/3)

dx    = gamma_range[int(len(gamma_range) / 2)] - gamma_range[int(len(gamma_range) / 2 - 1)]
dy    = r_in[int(len(r_in) / 2)] - r_in[int(len(r_in) / 2 - 1)]

plt.text(gamma_range[int(len(r_in) / 2)], r_in[int(len(r_in) / 2)] + 0.05, 
         r'$r_{in}$', 
         ha = 'left', 
         va = 'bottom',
         transform_rotates_text = True, 
         rotation = np.rad2deg(np.arctan2(dy, dx)), 
         rotation_mode = 'anchor', 
         fontsize = 16)

plt.plot(gamma_range, r_in, "orange")

# ============================================================================================================= #

dx    = gamma_r_h[int(len(gamma_r_h) / 2)] - gamma_r_h[int(len(gamma_r_h) / 2 - 1)]
dy    = r_h_p[int(len(r_h_p) / 2)] - r_h_p[int(len(r_h_p) / 2 - 1)]

plt.text(gamma_r_h[int(len(r_h_p) / 2)], r_h_p[int(len(r_h_p) / 2)] + 0.05, 
         r'$r_{h}^+$', 
         ha = 'left', 
         va = 'bottom',
         transform_rotates_text = True, 
         rotation = np.rad2deg(np.arctan2(dy, dx)), 
         rotation_mode = 'anchor', 
         fontsize = 16)

plt.plot(gamma_r_h, r_h_p, "k-")

# ============================================================================================================= #

dx    = gamma_r_h[int(len(gamma_r_h) / 2)] - gamma_r_h[int(len(gamma_r_h) / 2 - 1)]
dy    = r_h_m[int(len(r_h_m) / 2)] - r_h_m[int(len(r_h_m) / 2 - 1)]

plt.text(gamma_r_h[int(len(r_h_m) / 2)], r_h_m[int(len(r_h_m) / 2)] + 0.05, 
         r'$r_{h}^-$', 
         ha = 'left', 
         va = 'bottom',
         transform_rotates_text = True, 
         rotation = np.rad2deg(np.arctan2(dy, dx)), 
         rotation_mode = 'anchor', 
         fontsize = 16)

plt.plot(gamma_r_h, r_h_m, "k--")

plt.xlim([0, 2])
plt.ylim([0, 7])

plt.tick_params(labelsize = 26)
plt.xlabel(r"$\gamma\,\,[M^2]$", fontsize = 32)
plt.ylabel(r"$r\,\,[M]$", fontsize = 32)



plt.show()