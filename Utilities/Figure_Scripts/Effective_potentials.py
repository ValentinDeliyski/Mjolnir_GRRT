import matplotlib.pyplot as plt
import matplotlib
import numpy as np

def find_nearest(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or np.fabs(value - array[idx-1]) < np.fabs(value - array[idx])):
        return array[idx-1], idx - 1
    else:
        return array[idx], idx

def plot_JNW_potential(subplot: plt.subplot, gamma: float, impact_param: float, color, linestyle: str) -> None:

    r_singularity = 2 / gamma

    r_range = np.linspace(r_singularity, 11, 1000)

    potential = impact_param**2 / r_range**2 * pow(1 - r_singularity / r_range, 2 * gamma - 1)

    subplot.plot(r_range, potential, color = color, linestyle = linestyle)

    if gamma == gamma_range[-1]:
        
        x_coord = 6.3

        dx    = r_range[find_nearest(r_range, x_coord)[1]] - r_range[find_nearest(r_range, x_coord)[1] - 1]
        dy    = potential[find_nearest(r_range, x_coord)[1]] - potential[find_nearest(r_range, x_coord)[1] - 1]
        angle = np.rad2deg(np.arctan2(dy, dx))

        subplot.text(find_nearest(r_range, x_coord)[0], potential[find_nearest(r_range, x_coord)[1]] - 0.5, r'$\gamma=${}'.format(round(gamma, 2)), ha='left', va='bottom',
                     transform_rotates_text = True, rotation = angle, rotation_mode='anchor', fontsize = 16)
        
    if gamma == gamma_range[0]:
        x_coord = 6.5

        dx    = r_range[find_nearest(r_range, x_coord)[1]] - r_range[find_nearest(r_range, x_coord)[1] - 1]
        dy    = potential[find_nearest(r_range, x_coord)[1]] - potential[find_nearest(r_range, x_coord)[1] - 1]
        angle = np.rad2deg(np.arctan2(dy, dx))

        subplot.text(find_nearest(r_range, x_coord)[0], potential[find_nearest(r_range, x_coord)[1]], r'$\gamma=${}'.format(round(gamma, 2)), ha='left', va='bottom',
                     transform_rotates_text = True, rotation = 0.9*angle , rotation_mode='anchor', fontsize = 16)
        

    subplot.set_ylim([0, 10])
    subplot.set_xlim([1, 10])

def plot_GB_effective_potential(subplot: plt.subplot, gamma: float, impact_param: float, color, linestyle: str) -> None:


    r_range = np.linspace(0, 11, 1000)

    f = 1 + r_range**2 / 2 / gamma * (1 - np.sqrt(1 + 8 * gamma / r_range**3))

    potential = impact_param**2 / r_range**2 * f

    subplot.plot(r_range, potential, color = color, linestyle = linestyle)

    if gamma == gamma_range[-1]:
        
        x_coord = 4

        dx    = r_range[find_nearest(r_range, x_coord)[1]] - r_range[find_nearest(r_range, x_coord)[1] - 1]
        dy    = potential[find_nearest(r_range, x_coord)[1]] - potential[find_nearest(r_range, x_coord)[1] - 1]
        angle = np.rad2deg(np.arctan2(dy, dx))

        subplot.text(2.2, potential[find_nearest(r_range, 2.5)[1]] + 0.5, r'$\gamma=${}'.format(round(gamma, 2)), ha='left', va='bottom',
                     transform_rotates_text = True, rotation = angle, rotation_mode='anchor', fontsize = 16)
        
    if gamma == gamma_range[0]:

        x_coord = 4

        dx    = r_range[find_nearest(r_range, x_coord)[1]] - r_range[find_nearest(r_range, x_coord)[1] - 1]
        dy    = potential[find_nearest(r_range, x_coord)[1]] - potential[find_nearest(r_range, x_coord)[1] - 1]
        angle = np.rad2deg(np.arctan2(dy, dx))

        subplot.text(2.2, potential[find_nearest(r_range, 2.5)[1]]- 0.5, r'$\gamma=${}'.format(round(gamma, 2)), ha='left', va='bottom',
                     transform_rotates_text = True, rotation = angle, rotation_mode='anchor', fontsize = 16)

    subplot.set_ylim([0, 10])
    subplot.set_xlim([0, 10])

    return

Figure = plt.figure()
subplot = Figure.add_subplot(111)

cmap = matplotlib.colormaps['plasma']

gamma_range = np.linspace(0.5 + 1e-2, 1, 10)

for gamma in gamma_range:

    plot_JNW_potential(subplot = subplot, 
                       gamma = gamma, 
                       impact_param = 12,
                       color = cmap((gamma - gamma_range[0]) / ((gamma_range[-1] - gamma_range[0]) + 0.05)), 
                       linestyle = "-")
    


subplot.tick_params(axis='both', which='major', labelsize=14)
subplot.set_xlabel(r"$r\,\,[M]$", fontsize = 16)
subplot.set_ylabel(r"$V_\text{eff}(r)$",  fontsize = 16)

Figure = plt.figure()
subplot = Figure.add_subplot(111)

gamma_range = np.linspace(1.0, 3 * np.sqrt(3) / 4 - 1e-5, 10)

for gamma in gamma_range:

    plot_GB_effective_potential(subplot = subplot, 
                                gamma = gamma, 
                                impact_param = 12,
                                color = cmap((gamma - gamma_range[0]) / ((gamma_range[-1] - gamma_range[0]) + 0.05)), 
                                linestyle = "-")
    
    # subplot.plot([1, 1], [-10, 100])

subplot.tick_params(axis='both', which='major', labelsize=16)
subplot.set_xlabel(r"$r\,\,[M]$", fontsize = 16)
subplot.set_ylabel(r"$V_\text{eff}(r)$",  fontsize = 16)


plt.show()