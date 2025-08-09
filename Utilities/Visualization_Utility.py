from Visualizer_Class import Sim_Visualizer
import matplotlib.pyplot as plt
from Support_functions.Parsers import*

if __name__ == "__main__":
    
    fig_title = r"Numerical Config VII$^{-5}_{0.05}$"
            
    plt.rcParams['axes.titlepad'] = 20
    plt.rcParams['text.usetex'] = True
    plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

    EHT_Array           = ["ngEHT"]
    Sim_path            = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\Utilities\\Reference_simulations\\Galin_Numerical_config_III_4_230\\Numerical"
    Sim_Frequency_Bins  = ["230"] # In units of [GHz]

    Visualizer = Sim_Visualizer(Sim_path           = Sim_path, 
                                Sim_Frequency_Bins = Sim_Frequency_Bins,
                                Array              = EHT_Array,
                                Font_size          = 32, 
                                Label_Pad          = 8, 
                                Common_file_name   = "Numerical",
                                Respect_folder_structure = False)

    Visualizer.plot_ray_tracer_results(Export_data_for_Ehtim = False, 
                                       Save_Figures          = True, 
                                       Radiation_Component   = "Stokes I",
                                       Custom_fig_title      = fig_title,
                                       Obs_effective_distance = Visualizer.Units.M87_DISTANCE_GEOMETRICAL,
                                       Colormap = "hot")

    # Visualizer.plot_EHTIM_results(Make_contour_plots = False,                                                      
    #                               Contour_specs      = [], 
    #                               Save_Figures       = True,
    #                               Plot_no_blur       = False,
    #                               Custom_fig_title = fig_title) 
    
    # Visualizer.plot_EHTIM_results(Make_contour_plots = True,                                                      
    #                               Contour_specs      = [([0.15,  0.2,  0.3], ["r", "w", "k"]),
    #                                                     ([0.04, 0.06, 0.1], ["r", "w", "k"])],
    #                             Save_Figures       = True,
    #                             Plot_no_blur       = False,
    #                             Custom_fig_title = fig_title) 
                 
    # # Visualizer.plot_EHTIM_results(Make_contour_plots = False,                                                      
    # #                             Contour_specs      = [], 
    # #                             Save_Figures       = True,
    # #                             Plot_no_blur       = True,
    # #                             Custom_fig_title = fig_title) 
 
    # Visualizer.plot_VIDA_style(Center_plot = False, 
    #                            Save_Figures = True, 
    #                            Custom_fig_title = fig_title)

    # Visualizer.create_EHTIM_superposition()

    # Visualizer.plot_superposition(Center_plot = False,
    #                               Save_Figures = True,
    #                               Custom_fig_title = r"Gauss-Bonnet ($\gamma = 1.15$)")

    # Visualizer.compare_superpos_w_single_freq(Contour_specs = [([0.12, 0.2, 0.3], ["b", "r", "w"]),
    #                                                            ([0.12, 0.2, 0.3], ["b", "r", "w"]),
    #                                                            ([0.075, 0.14, 0.3], ["b", "r", "w"])],
    #                                         Save_Figures = True, 
    #                                         Custom_fig_title = r"Gauss-Bonnet ($\gamma = 1.15$)")

    # Visualizer.save_console_log_to_file()

               
    # plt.close("all")
    plt.show()
