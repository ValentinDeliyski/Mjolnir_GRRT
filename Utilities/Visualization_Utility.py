from Visualizer_Class import Sim_Visualizer
import matplotlib.pyplot as plt
from Support_functions.Parsers import Units_class


if __name__ == "__main__":
    
    plt.rcParams['axes.titlepad'] = 20

    EHT_Array           = []
    Sim_path            = "C:\\Users\\Valur\\Documents\\Repos\\Gravitational_Lenser\\Sim_Results\\Gauss_Bonnet"
    Sim_Frequency_Bins  = ["230"] # In units of [GHz]

    Visualizer = Sim_Visualizer(Sim_path           = Sim_path, 
                                Sim_Frequency_Bins = Sim_Frequency_Bins,
                                Array              = EHT_Array,
                                Units_class_inst   = Units_class(),
                                Font_size = 16, 
                                Label_Pad = 8,
                                Respect_folder_structure = False)

    Visualizer.plot_ray_tracer_results(Export_data_for_Ehtim = False, 
                                       Save_Figures = True, 
                                       Stokes_component = "I",
                                       Custom_fig_title = r"Wormhole ($\gamma = 2,\,\,\, a = 0.5$)")    
    
    # Visualizer.plot_EHTIM_results(Make_contour_plots = False,                                                      
    #                               Contour_specs      = [([0.02, 0.1], ["r", "w"])], 
    #                               Save_Figures       = True,
    #                               Plot_no_blur       = False,
    #                               Custom_fig_title = r"Schwarzschild") 
    
    # Visualizer.plot_EHTIM_results(Make_contour_plots = True,                                                      
    #                               Contour_specs      = [([0.02, 0.1], ["r", "w"])], 
    #                               Save_Figures       = True,
    #                               Plot_no_blur       = False,
    #                               Custom_fig_title = r"Schwarzschild") 
    

    # Visualizer.plot_EHTIM_results(Make_contour_plots = False,                                                      
    #                               Contour_specs      = [([0.02, 0.1], ["r", "w"])], 
    #                               Save_Figures       = True,
    #                               Plot_no_blur       = True,
    #                               Custom_fig_title = r"Schwarzschild") 

    # Visualizer.plot_VIDA_style(Center_plot = False, 
    #                            Save_Figures = False, 
    #                            Custom_fig_title = r"Kerr ($\text{a} = 0.9$)")

    # Visualizer.create_EHTIM_superposition()

    # Visualizer.plot_superposition(Center_plot = False,
    #                               Save_Figures = False,
    #                               Custom_fig_title = r"Kerr ($\text{a} = 0.9$)")

    # Visualizer.compare_superpos_w_single_freq(Contour_specs = [([0.12, 0.15, 0.2, 0.3], ["k", "r", "b", 'w']),
    #                                                            ([0.12, 0.15, 0.2, 0.3], ["k", "r", "b", 'w']),
    #                                                            ([0.095, 0.14], ["r", "w"])],
    #                                         Save_Figures = False, 
    #                                         Custom_fig_title = r"Kerr ($\text{a} = 0.9$)")

    # Visualizer.save_console_log_to_file()

    # print(Visualizer.Sim_Parsers[0][0].Source_R_Coord)

    plt.rcParams['text.usetex'] = True
    plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
    plt.show()
