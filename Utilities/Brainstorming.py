import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from scipy.optimize import brentq

def JNW_ms_orbit_condition(r: float, gamma: float) -> float:

    b = 2 / gamma

    if r < b:
        return 100

    return (1 - gamma) * (1 - b / r)**gamma * b + 2 * (1 - b / r)**(gamma + 1) * r + (2 * gamma + 1) * b - 2 * r

def get_JNW_ISCO(gamma: float, direction: int) -> float:

    return 1 / gamma * (3 * gamma + 1 + direction * np.sqrt(5 * gamma**2 - 1))

def get_JNW_r_ph(gamma: float) -> float:

    b = 2 / gamma

    return (2 * gamma + 1) * b / 2

gamma_range = np.linspace(1 / np.sqrt(5) + 1e-10, 1, 10000)

r_ms = []
r_ph = []
r_ISCO_outer = []
r_ISCO_inner = []
gamma_plot_r_ms = []
gamma_plot_r_ph = []
gamma_plot_r_ISCO_inner = []

for gamma in gamma_range:

    r_ISCO_outer.append(get_JNW_ISCO(gamma, 1))

    solution = get_JNW_ISCO(gamma, -1)

    if solution > 2 / gamma:
        r_ISCO_inner.append(solution)
        gamma_plot_r_ISCO_inner.append(gamma)

    solution = get_JNW_r_ph(gamma)

    if solution > 2 / gamma:
        gamma_plot_r_ph.append(gamma)
        r_ph.append(solution)

    try:
        solution = brentq(JNW_ms_orbit_condition, args = gamma, a = 4, b = 10, xtol= 1e-10)

        if solution > (2 / gamma) + 1e-10:
            r_ms.append(solution)
            gamma_plot_r_ms.append(gamma)

    except:
        continue

cmap = matplotlib.colormaps['plasma']

# ================================================================================================================================================= #

gamma_plot_r_ISCO_inner.insert(0, gamma_range[0])
r_ISCO_inner.insert(0, r_ISCO_outer[0])

dx    = gamma_plot_r_ISCO_inner[int(len(r_ISCO_inner) / 4)] - gamma_plot_r_ISCO_inner[int(len(r_ISCO_inner) / 4 - 1)]
dy    = r_ISCO_inner[int(len(r_ISCO_inner) / 4)] - r_ISCO_inner[int(len(r_ISCO_inner) / 4 - 1)]
angle = np.rad2deg(np.arctan2(dy, dx))

plt.text(gamma_plot_r_ISCO_inner[int(len(r_ISCO_inner) / 4)], r_ISCO_inner[int(len(r_ISCO_inner) / 4)] + 0.1, r'$r_-$', ha='left', va='bottom',
         transform_rotates_text=True, rotation=angle, rotation_mode='anchor')

plt.plot(gamma_plot_r_ISCO_inner, r_ISCO_inner, color = "black")

plt.plot(gamma_plot_r_ISCO_inner[-1], r_ISCO_inner[-1], "o", color = "red")

# ================================================================================================================================================= #

dx    = gamma_range[int(len(r_ISCO_outer) / 2)] - gamma_range[int(len(r_ISCO_outer) / 2 - 1)]
dy    = r_ISCO_outer[int(len(r_ISCO_outer) / 2)] - r_ISCO_outer[int(len(r_ISCO_outer) / 2 - 1)]
angle = np.rad2deg(np.arctan2(dy, dx))

plt.text(gamma_range[int(len(r_ISCO_outer) / 2)], r_ISCO_outer[int(len(r_ISCO_outer) / 2)], r'$r_+$', ha='left', va='bottom',
         transform_rotates_text=True, rotation=angle, rotation_mode='anchor')

plt.plot(gamma_range, r_ISCO_outer, color = "black")

plt.plot(gamma_range[0], r_ISCO_outer[0], "o", color = "red")

# ================================================================================================================================================= #

dx    = gamma_range[int(len(r_ISCO_outer) / 2)]
angle = np.rad2deg(np.arctan2(-2, gamma_range[int(len(r_ISCO_outer) / 2)] ** 2))

plt.text(gamma_range[int(len(r_ISCO_outer) / 2)], 2 / gamma_range[int(len(r_ISCO_outer) / 2)], r'$r_{\text{cs}}$', ha='left', va='bottom',
         transform_rotates_text=True, rotation=angle, rotation_mode='anchor')

plt.plot(gamma_range, 2 / gamma_range, "--", color = 'red')

# ================================================================================================================================================= #

dx    = gamma_plot_r_ms[int(len(r_ms) / 1.5)] - gamma_plot_r_ms[int(len(r_ms) / 1.5 - 1)]
dy    = r_ms[int(len(r_ms) / 1.5)] - r_ms[int(len(r_ms) / 1.5 - 1)]
angle = np.rad2deg(np.arctan2(dy, dx))

plt.text(gamma_range[int(len(r_ISCO_outer) / 2)], r_ms[int(len(r_ms) / 1.5)] + 0.15, r'$r_\text{ms}$', ha='left', va='bottom',
         transform_rotates_text=True, rotation=angle, rotation_mode='anchor')

plt.plot(gamma_plot_r_ms, r_ms)

plt.plot(gamma_plot_r_ms[0], r_ms[0], "o", color = "blue")

# ================================================================================================================================================= #

dx    = gamma_plot_r_ph[int(len(r_ms) / 1.5)] - gamma_plot_r_ph[int(len(r_ms) / 1.5 - 1)]
dy    = r_ph[int(len(r_ms) / 1.5)] - r_ph[int(len(r_ms) / 1.5 - 1)]
angle = np.rad2deg(np.arctan2(dy, dx))

plt.text(gamma_range[int(len(r_ISCO_outer) / 2)], r_ph[int(len(r_ph) / 1.5)] + 0.16, r'$r_\text{ph}$', ha='left', va='bottom',
         transform_rotates_text=True, rotation=angle, rotation_mode='anchor')

plt.plot(gamma_plot_r_ph, r_ph)

plt.plot(gamma_plot_r_ph[0], r_ph[0], "o", color = "red")

plt.xlabel(r"$\gamma\,[-]$", fontsize = 16)
plt.ylabel(r"$r_\text{mb},\, [M]$", fontsize = 16)

plt.xlim([1 / np.sqrt(5) - 0.01, 1])

plt.show()