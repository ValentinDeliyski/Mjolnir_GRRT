from Visualizer_Class import Sim_Visualizer
import matplotlib.pyplot as plt
from Support_functions.Parsers import*

if __name__ == "__main__":
    
    plt.rcParams['axes.titlepad'] = 20

    EHT_Array           = ["ngEHT"]
    Sim_path            = "C:\\Users\\Valur\\Documents\\Papers\\Sim_Paper_2\\Wormhole_a_0.5_alpha_0\\"
    Sim_Frequency_Bins  = ["230"] # In units of [GHz]

    Visualizer = Sim_Visualizer(Sim_path           = Sim_path, 
                                Sim_Frequency_Bins = Sim_Frequency_Bins,
                                Array              = EHT_Array,
                                Font_size          = 32, 
                                Label_Pad          = 8, 
                                Common_file_name   = "Wormhole",
                                Respect_folder_structure = True)

    Visualizer.plot_ray_tracer_results(Export_data_for_Ehtim = False, 
                                       Save_Figures          = False, 
                                       Radiation_Component   = "Stokes I",
                                       Custom_fig_title      = r"Wormhole $a= 0.5$ $\alpha = 0$")

    # Visualizer.plot_EHTIM_results(Make_contour_plots = False,                                                      
    #                               Contour_specs      = None, 
    #                               Save_Figures       = False,
    #                               Plot_no_blur       = False,
    #                               Custom_fig_title = r"Wormhole $a= 0$ $\alpha = 2$") 
    
    # Visualizer.plot_EHTIM_results(Make_contour_plots = True,                                                      
    #                               Contour_specs      = [([0.16, 0.2, 0.3], ["r", "w", "k"])], 
    #                               Save_Figures       = True,
    #                               Plot_no_blur       = False,
    #                               Custom_fig_title = r"Wormhole $a a= 0.5$ $\alpha = 0$") 
    
    Visualizer.plot_EHTIM_results(Make_contour_plots = False,                                                      
                                  Contour_specs      = None, 
                                  Save_Figures       = False,
                                  Plot_no_blur       = True,
                                  Custom_fig_title = r"Wormhole $a= 0.5$ $\alpha = 0$") 

    Visualizer.plot_VIDA_style(Center_plot = False, 
                               Save_Figures = True, 
                               Custom_fig_title = r"Wormhole $a= 0$ $\alpha = 2$")

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

    plt.rcParams['text.usetex'] = True
    plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
    plt.show()
