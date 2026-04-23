from Support_functions.Surface_Cubic_B_spline import Surface_Cubic_B_spline
from Support_functions.Parsers import Units_class, Simulation_Parser
from Metric_spline_creator import Numerical_metric_parser_class
import matplotlib.pyplot as plt
from matplotlib import colors

from numpy import pi, dot, float64, cosh, sinh, exp, array, flip, cos, sin, ones, arccos, arctan, sqrt, linspace, zeros
from numpy.linalg import norm
from numpy.typing import NDArray

import os 

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):

    new_cmap = colors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
                                                        cmap(linspace(minval, maxval, n)))
    
    return new_cmap

def split_EVPA(EVPA):

        DISCONTINUITY_TRESHOLD = pi / 2

        branch_index = []
        branch_index.append(0)

        for index, _ in enumerate(EVPA):

            if index > 1 and abs(EVPA[index] - EVPA[index - 1]) > DISCONTINUITY_TRESHOLD:

                branch_index.append(index)

        branch_index.append(len(EVPA) - 1)

        return branch_index


if __name__ == "__main__":
    
    plt.rcParams['axes.titlepad'] = 15
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = 'serif' # or 'sans-serif' or 'monospace'
    plt.rcParams['font.serif'] = 'cmr10'
    plt.rcParams['font.sans-serif'] = 'cmss10'
    plt.rcParams['font.monospace'] = 'cmtt10'
    plt.rcParams["axes.formatter.use_mathtext"] = True
    plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}' + '\n' + r'\usepackage{xcolor}'

    Numerical_metric_parser = Numerical_metric_parser_class(File_path = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\Utilities\\Numerical_metrics\\Zero_curvature\\rh=0.01_om=0.835271834271409_h0=0.0580074274088523.dat",
                                                            Grid_R_size = 120,
                                                            Grid_Theta_size = 30)
    
    x_coord, theta_coord, F_0, F_1, F_2, W, Scalar_field = Numerical_metric_parser.get_parsed_results()

    Scalar_field_spline_instance = Surface_Cubic_B_spline(x_grid = theta_coord, 
                                                          y_grid = x_coord, 
                                                          z_grid = Scalar_field, 
                                                          X_patch_number = 2 * Numerical_metric_parser.GRID_THETA_SIZE - 1, 
                                                          Y_patch_number = Numerical_metric_parser.GRID_R_SIZE)

    Max_scalar_field = max(Scalar_field.flatten())

    Scalar_field_grid = zeros((1000, 1000))
    
    for rho_idx, rho in enumerate(linspace(-1, 1, 1000)):
        
        for z_idx, z in enumerate(linspace(-1, 1, 1000)):
  
            r = sqrt(rho**2 + z**2)
            
            if r > 0.01:
            
                theta = arccos(z / r)
                
                x = sqrt(r**2 - 0.01**2)
                x = x / (x + 1)
                
                Scalar_field_grid[z_idx][rho_idx] = Scalar_field_spline_instance.evaluate_spline_single_point(x = theta, 
                                                                                                              y = x, 
                                                                                                              x_grid_len = 2 * Numerical_metric_parser.GRID_THETA_SIZE - 1,
                                                                                                              y_grid_len = Numerical_metric_parser.GRID_R_SIZE) / Max_scalar_field
                
    Sim_parser_Kerr = Simulation_Parser("C:/Users/Valur/Documents/Repos/Mjolnir_GRRT/Utilities/Numerical_metric_runs/Zero_curvature/Config_II_70_deg_B_0.0_1.0_0.0_min_order_0_max_more_order_0_zoom/Kerr_analog_single_geodesic/Kerr_photon_log")
    Sim_parser_Numerical = Simulation_Parser("C:/Users/Valur/Documents/Repos/Mjolnir_GRRT/Utilities/Numerical_metric_runs/Zero_curvature/Config_II_70_deg_B_0.0_1.0_0.0_min_order_0_max_more_order_0_zoom/Numerical_single_geodesic/Numerical_photon_log")
    
    Position_tuple_Kerr, _, _, _, Affine_param_Kerr, _, Polarization_tuple_Kerr = Sim_parser_Kerr.get_photon_log()
    Position_tuple_Numerical, _, _, _, Affine_param_Numerical, _, Polarization_tuple_Numerical = Sim_parser_Numerical.get_photon_log()
    
    # ==================================================================  Polarization plotting ===================================================================
    
    t_Kerr, r_Kerr, theta_Kerr, phi_Kerr = Position_tuple_Kerr
    t_Numerical, r_Numerical, theta_Numerical, phi_Numerical = Position_tuple_Numerical
    
    M = float(Sim_parser_Kerr.Simulation_metadata["Mass [M]"])
    a = float(Sim_parser_Kerr.Simulation_metadata["Spin Parameter [M]"])
    
    R_H = M + sqrt(M**2 - a**2)
    r_Kerr = (array(r_Kerr) - a**2 / R_H) / M
    r_Numerical = array(r_Numerical) / M
    
    Pol_x_Kerr, Pol_y_Kerr = Polarization_tuple_Kerr
    EVPA_Kerr = arctan(-array(Pol_x_Kerr) / array(Pol_y_Kerr))
    
    Pol_x_Numerical, Pol_y_Numerical = Polarization_tuple_Numerical
    EVPA_Numerical = arctan(-array(Pol_x_Numerical) / array(Pol_y_Numerical))
    
    Figure, (Trajectory_Plot, EVPA_subplot) = plt.subplots(1, 2, gridspec_kw = {'width_ratios': [1, 1]}, layout = "compressed", figsize = (8, 4))
    
    Pol_start_idx_Kerr = (array(Pol_x_Kerr) != 0).argmax(axis = 0)
    Pol_start_idx_Numerical = (array(Pol_x_Numerical) != 0).argmax(axis = 0)
    Start_r = min(r_Kerr[Pol_start_idx_Kerr], r_Numerical[Pol_start_idx_Numerical])
    
    EVPA_Kerr_idx_split = split_EVPA(EVPA = EVPA_Kerr)
    EVPA_Numerical_idx_split = split_EVPA(EVPA = EVPA_Numerical)
    
    EVPA_subplot.plot([-100, -100], [-100, -99], "C0")
    EVPA_subplot.plot([-100, -100], [-100, -99], "m")
    
    EVPA_subplot.legend(["Kerr Analog", r"Model II$^0_{0.01}$"], loc = "lower right", fontsize = 12)
                
    for idx in range(len(EVPA_Kerr_idx_split) - 1):
        EVPA_subplot.plot(r_Kerr[EVPA_Kerr_idx_split[idx] : EVPA_Kerr_idx_split[idx + 1]], EVPA_Kerr[EVPA_Kerr_idx_split[idx] : EVPA_Kerr_idx_split[idx + 1]], "C0")       
                
    for idx in range(len(EVPA_Numerical_idx_split) - 1):        
        EVPA_subplot.plot(r_Numerical[EVPA_Numerical_idx_split[idx] : EVPA_Numerical_idx_split[idx + 1]], EVPA_Numerical[EVPA_Numerical_idx_split[idx] : EVPA_Numerical_idx_split[idx + 1]], "m")
         
    for idx in range(len(EVPA_Kerr_idx_split) - 1):
        EVPA_subplot.plot(r_Kerr[EVPA_Kerr_idx_split[idx + 1] - 1 : EVPA_Kerr_idx_split[idx + 1] + 1], EVPA_Kerr[EVPA_Kerr_idx_split[idx + 1] - 1 : EVPA_Kerr_idx_split[idx + 1] + 1], "C0--")               
    
    for idx in range(len(EVPA_Numerical_idx_split) - 1):
        EVPA_subplot.plot(r_Numerical[EVPA_Numerical_idx_split[idx + 1] - 1 : EVPA_Numerical_idx_split[idx + 1] + 1], EVPA_Numerical[EVPA_Numerical_idx_split[idx + 1] - 1 : EVPA_Numerical_idx_split[idx + 1] + 1], "m--")      

    # EVPA_subplot.plot(r_Numerical, EVPA_Numerical, "m")
    # EVPA_subplot.plot(r_Kerr, EVPA_Kerr, "C0")
    
    EVPA_subplot.plot([-100, 100], [EVPA_Kerr[-1], EVPA_Kerr[-1]], "k--")
    EVPA_subplot.plot([-100, 100], [EVPA_Numerical[-1], EVPA_Numerical[-1]], "k--")
    
    EVPA_subplot.plot(r_Kerr[Pol_start_idx_Kerr], EVPA_Kerr[Pol_start_idx_Kerr], "ko")
    EVPA_subplot.plot(r_Numerical[Pol_start_idx_Numerical], EVPA_Numerical[Pol_start_idx_Numerical], "ko")
    
    EVPA_subplot.set_xlim([0, Start_r + 10])
    EVPA_subplot.set_ylim([-1.1 * pi / 2, 1.1  * pi / 2]) 
    EVPA_subplot.set_ylabel(r"EVPA [rad]", fontsize = 18)
    EVPA_subplot.set_xlabel(r"$r$ [$M_\text{ADM}$]", fontsize = 18)
    
    EVPA_subplot.tick_params(which = 'minor', length = 4, labelsize = 18)
    EVPA_subplot.tick_params(which = 'major', length = 8, labelsize = 18)
    EVPA_subplot.set_title("EVPA Evolution", fontsize = 18)
    
    y_ticks = [-pi / 2, -pi / 4, 0, pi / 4, pi / 2]
    y_tick_labels = [r"$-\frac{\pi}{2}$", r"$-\frac{\pi}{4}$", "0", r"$\frac{\pi}{4}$", r"$\frac{\pi}{2}$"]
                
    EVPA_subplot.set_yticks(y_ticks)
    EVPA_subplot.set_yticklabels(y_tick_labels) 
    
    r_Kerr = r_Kerr[Pol_start_idx_Kerr:-1]
    theta_Kerr = theta_Kerr[Pol_start_idx_Kerr:-1]
    
    r_Numerical = r_Numerical[Pol_start_idx_Numerical:-1]
    theta_Numerical = theta_Numerical[Pol_start_idx_Numerical:-1]
    
    Trajectory_Plot.plot(r_Kerr * sin(theta_Kerr), r_Kerr * cos(theta_Kerr), "C0")
    Trajectory_Plot.plot([-100, 100], [0, 0], "k--")
    Trajectory_Plot.plot(r_Numerical * sin(theta_Numerical), r_Numerical * cos(theta_Numerical), "m")
    
    Trajectory_Plot.plot(r_Kerr[0] * sin(theta_Kerr[0]), r_Kerr[0] * cos(theta_Kerr[0]), "ko")
    Trajectory_Plot.plot(r_Numerical[0] * sin(theta_Numerical[0]), r_Numerical[0] * cos(theta_Numerical[0]), "ko")
    
    Trajectory_Plot.set_xlim([-0.05, 0.4])
    Trajectory_Plot.set_ylim([-0.05, 0.4])
    
    x_ticks = [0, 0.2, 0.4]
    x_tick_labels = ["0", "0.2", "0.4"]
    
    Trajectory_Plot.set_xticks(x_ticks)
    Trajectory_Plot.set_xticklabels(x_tick_labels)
    
    Trajectory_Plot.set_xlabel(r"$r\sin\theta$ [$M_\text{ADM}$]", fontsize = 18)
    Trajectory_Plot.set_ylabel(r"$r\cos\theta$ [$M_\text{ADM}$]", fontsize = 18)
    
    Trajectory_Plot.tick_params(which = 'minor', length = 4, labelsize = 18)
    Trajectory_Plot.tick_params(which = 'major', length = 8, labelsize = 18) 
    Trajectory_Plot.set_title("Geodesic Trajectories", fontsize = 18)
    
    Horizon_circle = plt.Circle((0.0, 0.0), 0.01 / M, color = "k") # type: ignore
    Trajectory_Plot.add_patch(Horizon_circle)
    
    cmap = plt.get_cmap('Reds')
    new_cmap = truncate_colormap(cmap, 0, 0.6)
    
    Scalar_field_plot = Trajectory_Plot.imshow(Scalar_field_grid, extent = array([-1, 1, -1, 1]) / M, cmap = new_cmap)
    colorbar = Figure.colorbar(Scalar_field_plot, ax = Trajectory_Plot, fraction = 0.046, pad = 0.04)
    colorbar.ax.set_title(r"$\frac{\phi}{\max\phi}$", fontsize = 18)
    colorbar.set_ticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    
    colorbar.ax.tick_params(labelsize = 18)
    
    plt.show()