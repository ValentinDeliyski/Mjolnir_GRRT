import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

def get_density_profile(rho: float, 
                        z: float, 
                        tan_theta: float = 0.1,
                        R_0: float = 4.5, 
                        R_cut: float = 4.5, 
                        Cutoff_Scale: float = 0.4) -> float:

    r = np.sqrt(rho**2 + z**2)

    density_profile = (R_0 / r)**2 * np.exp(-(z / (rho * tan_theta))**2 / 2)

    if (r < R_cut):
        density_profile = density_profile * np.exp(-((r - R_cut)/Cutoff_Scale)**2)

    return density_profile

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):

    new_cmap = colors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
                                                        cmap(np.linspace(minval, maxval, n)))
    
    return new_cmap

def make_disk_profile_plot(Scale: float = 10, 
                           Resolution: int = 500,  
                           tan_theta: float = 0.1,
                           R_0: float = 4.5, 
                           R_cut: float = 4.5, 
                           Cutoff_Scale: float = 0.4,
                           Save_figure: bool = False):

    r_sch = 2
    Grid = np.zeros((Resolution, Resolution))

    # ----------- Compute the density profile on the grid ----------- #

    for left_idx in range(Resolution):

        for right_idx in range(Resolution):

            rho = right_idx / (Resolution - 1) * 2 * Scale
            z   = Scale * (1 - left_idx / (Resolution / 2))

            if rho**2 + z**2 > r_sch**2 and rho > 0:

                Grid[left_idx][right_idx] = get_density_profile(rho = rho, 
                                                                z = z, 
                                                                tan_theta = tan_theta, 
                                                                R_0 = R_0, 
                                                                R_cut = R_cut, 
                                                                Cutoff_Scale = Cutoff_Scale)

    # ----------- Create the Schwawrzschild event horizon plot ----------- #

    angle = np.linspace(-np.pi / 2, np.pi / 2, 1000)            
    
    rho_horizon = r_sch * np.cos(angle)
    z_horizon   = r_sch * np.sin(angle) 

    # ----------- Plot the density profile and the horizon ----------- #

    cmap = plt.get_cmap('seismic')
    new_cmap = truncate_colormap(cmap, 0.5, 1)

    Main_Figure = plt.figure()

    Subplot = Main_Figure.add_subplot(111)

    Density_profile = Subplot.imshow(Grid, cmap = new_cmap, extent=(0, 2 * Scale, - Scale, Scale), vmin = 0, vmax = 1)
    Horizon_plot = Subplot.plot(rho_horizon, z_horizon, "k--")
    Subplot.set_xlabel(r"$\rho,\,[M]$")
    Subplot.set_ylabel(r"$z,\,[M]$")

    colorbar = Main_Figure.colorbar(Density_profile, ax = Subplot, fraction=0.046, pad=0.04)
    colorbar.set_label(r'$n(\vec{r}\,)/n_0$')

    if Save_figure:

        Main_Figure.savefig("Disc_Contour_plot_" + str(tan_theta) + ".png", bbox_inches = 'tight')

    plt.show()

make_disk_profile_plot(tan_theta = 0.5, Save_figure = True)