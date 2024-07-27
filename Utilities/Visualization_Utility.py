from Visualizer_Class import Sim_Visualizer
import matplotlib.pyplot as plt
from Support_functions.Parsers import*

if __name__ == "__main__":
    
    plt.rcParams['axes.titlepad'] = 20

    EHT_Array           = ["2017"]
    # Sim_path            = "C:\\Users\\Valur\\Documents\\Repos\\Gravitational_Lenser\\Sim_Results\\Wormhole"
    Sim_path            = "C:\\Users\\Valur\\Documents\\Papers\\Sim_Paper_2\\Wormhole_a_0.9_alpha_2\\"
    Sim_Frequency_Bins  = ["230", "345"] # In units of [GHz]

    Visualizer = Sim_Visualizer(Sim_path           = Sim_path, 
                                Sim_Frequency_Bins = Sim_Frequency_Bins,
                                Array              = EHT_Array,
                                Font_size          = 24, 
                                Label_Pad          = 8,
                                Common_file_name   = "Wormhole",
                                Respect_folder_structure = True)

    Visualizer.plot_ray_tracer_results(Export_data_for_Ehtim = True, 
                                       Save_Figures          = False, 
                                       Radiation_Component   = "Stokes I",
                                       Custom_fig_title      = r"Kerr ($a / M = 0.5$)")        

    # Visualizer.plot_EHTIM_results(Make_contour_plots = False,                                                      
    #                               Contour_specs      = [([0.02, 0.1], ["r", "w"])], 
    #                               Save_Figures       = True,
    #                               Plot_no_blur       = False,
    #                               Custom_fig_title = r"Kerr ($a / M = 0.5$)") 
    
    # Visualizer.plot_EHTIM_results(Make_contour_plots = True,                                                      
    #                               Contour_specs      = [([0.02, 0.1], ["r", "w"])], 
    #                               Save_Figures       = True,
    #                               Plot_no_blur       = False,
    #                               Custom_fig_title = r"Kerr ($a / M = 0.5$)") 
    
    # Visualizer.plot_EHTIM_results(Make_contour_plots = False,                                                      
    #                               Contour_specs      = [([0.02, 0.1], ["r", "w"])], 
    #                               Save_Figures       = True,
    #                               Plot_no_blur       = True,
    #                               Custom_fig_title = r"Kerr ($a / M = 0.5$)") 

    # Visualizer.plot_VIDA_style(Center_plot = False, 
    #                            Save_Figures = False, 
    #                            Custom_fig_title = r"Kerr ($a / M = 0.5$)")

    # Visualizer.create_EHTIM_superposition()

    # Visualizer.plot_superposition(Center_plot = False,
    #                               Save_Figures = True,
    #                               Custom_fig_title = r"Kerr ($a / M = 0.5$)")

    # Visualizer.compare_superpos_w_single_freq(Contour_specs = [([0.004, 0.02, 0.1], ["b", "r", "w"]),
    #                                                            ([0.004, 0.02, 0.1], ["b", "r", "w"]),
    #                                                            ([0.004, 0.02, 0.1], ["b", "r", "w"])],
    #                                         Save_Figures = True, 
    #                                         Custom_fig_title = r"Kerr ($a / M = 0.5$)")

    # Visualizer.save_console_log_to_file()

    plt.rcParams['text.usetex'] = True
    plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
    plt.show()
