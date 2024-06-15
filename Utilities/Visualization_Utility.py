from Visualizer_Class import Sim_Visualizer
import matplotlib.pyplot as plt
from Support_functions.Parsers import*

if __name__ == "__main__":
    
    plt.rcParams['axes.titlepad'] = 20

    EHT_Array           = []
    Sim_path            = "C:\\Users\\Valur\\Documents\\Repos\\Gravitational_Lenser\\Sim_Results\\Wormhole"
    Sim_Frequency_Bins  = ["230"] # In units of [GHz]

    Visualizer = Sim_Visualizer(Sim_path           = Sim_path, 
                                Sim_Frequency_Bins = Sim_Frequency_Bins,
                                Array              = EHT_Array,
                                Units_class_inst   = Units_class(),
                                Font_size = 16, 
                                Label_Pad = 8,
                                Respect_folder_structure = False)

    Visualizer.plot_ray_tracer_results(Export_data_for_Ehtim = False, 
                                       Save_Figures = False, 
                                       Radiation_Component = "NT",
                                       Custom_fig_title = r"Wormhole ($\gamma = 2,\,\,\, a = 0$)")    
    
    # I_Intensity, Q_Intensity, U_Intensity, V_Intensity, NT_Flux, NT_Redshift, NT_Flux_Shifted = Simulation_Parser(Sim_path).get_plottable_sim_data()

    # plt.imshow(NT_Flux_Shifted, cmap = "hot")
    
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


    # """ =================================== Test 1 =================================== """

    # Test_Data = []

    # for test in [1,2,3,4,5]:

    #     Test_1_Sim_Path = "C:\\Users\\Valur\\Documents\\Repos\\Gravitational_Lenser\\Sim_Results\\{}\\Kerr".format(test)

    #     Visualizer = Sim_Visualizer(Sim_path           = Test_1_Sim_Path, 
    #                                 Sim_Frequency_Bins = Sim_Frequency_Bins,
    #                                 Array              = EHT_Array,
    #                                 Units_class_inst   = Units_class(),
    #                                 Font_size = 16, 
    #                                 Label_Pad = 8,
    #                                 Respect_folder_structure = False)
        
    #     I_Intensity_0, Q_Intensity_0, U_Intensity_0, V_Intensity_0, NT_Redshift_n0, NT_Flux_n0, NT_Flux_Shifted_n0 = Visualizer.Sim_Parsers[0][0].get_plottable_sim_data()
    #     I_Intensity_1, Q_Intensity_1, U_Intensity_1, V_Intensity_1, NT_Redshift_n1, NT_Flux_n1, NT_Flux_Shifted_n1 = Visualizer.Sim_Parsers[0][1].get_plottable_sim_data()
    #     I_Intensity_2, Q_Intensity_2, U_Intensity_2, V_Intensity_2, NT_Redshift_n2, NT_Flux_n2, NT_Flux_Shifted_n2 = Visualizer.Sim_Parsers[0][2].get_plottable_sim_data()
    #     I_Intensity_3, Q_Intensity_3, U_Intensity_3, V_Intensity_3, NT_Redshift_n3, NT_Flux_n3, NT_Flux_Shifted_n3 = Visualizer.Sim_Parsers[0][3].get_plottable_sim_data()


    #     Test_Data.append(I_Intensity_0 + I_Intensity_1 + I_Intensity_2 + I_Intensity_3)


    # Figure_Pattern, (Test_1_plot, Test_2_plot, Test_3_plot, Test_4_plot, Test_5_plot) = plt.subplots(1, 5, gridspec_kw = {'width_ratios': [1, 1, 1, 1, 1]}, constrained_layout = True)

    # axes_limits = np.array([(limit) for limit in Visualizer.Sim_Parsers[0][0].WINDOW_LIMITS]) * Visualizer.Sim_Parsers[0][0].OBS_DISTANCE / Visualizer.Units.SGRA_DISTANCE_GEOMETRICAL * Visualizer.Units.RAD_TO_MICRO_AS

    # # The literature (for some reason) has the X axis going positive to negative, 
    # # so I invert the X axis limits
    # axes_limits[0] = -axes_limits[0]
    # axes_limits[1] = -axes_limits[1]

    # Figure_Pattern.subplots_adjust()
    # Figure_Pattern.set_figwidth(40)
    # Figure_Pattern.set_figheight(8)

    # Cmap_max = max(np.abs(Test_Data[0].flatten()))
    # Cmap_min = 0
    # Test_1_plot.imshow(Test_Data[0], interpolation = 'bilinear', cmap = "plasma", extent = axes_limits, vmin = Cmap_min, vmax = Cmap_max)

    # Test_1_plot.set_xlabel(r'$\alpha_{rel}\,\,[\mu$as]', fontsize = 32)
    # Test_1_plot.set_ylabel(r'$\delta_{rel}\,\,[\mu$as]', fontsize = 32)   
    # Test_1_plot.tick_params(axis='x', labelsize=32)
    # Test_1_plot.tick_params(axis='y', labelsize=32)
    # Test_1_plot.set_title("Test 1", fontsize = 32)
    
    # Cmap_max = max(np.abs(Test_Data[1].flatten()))
    # Cmap_min = 0
    # Test_2_plot.imshow(Test_Data[1], interpolation = 'bilinear', cmap = "plasma", extent = axes_limits, vmin = Cmap_min, vmax = Cmap_max)

    # Test_2_plot.set_xlabel(r'$\alpha_{rel}\,\,[\mu$as]', fontsize = 32)
    # Test_2_plot.set_ylabel(r'$\delta_{rel}\,\,[\mu$as]', fontsize = 32)
    # Test_2_plot.tick_params(axis='x', labelsize=32)
    # Test_2_plot.tick_params(axis='y', labelsize=32)
    # Test_2_plot.set_title("Test 2", fontsize = 32)

    # Cmap_max = max(np.abs(Test_Data[2].flatten()))
    # Cmap_min = 0
    # Test_3_plot.imshow(Test_Data[2], interpolation = 'bilinear', cmap = "plasma", extent = axes_limits, vmin = Cmap_min, vmax = Cmap_max)

    # Test_3_plot.set_xlabel(r'$\alpha_{rel}\,\,[\mu$as]', fontsize = 32)
    # Test_3_plot.set_ylabel(r'$\delta_{rel}\,\,[\mu$as]', fontsize = 32)
    # Test_3_plot.tick_params(axis='x', labelsize=32)
    # Test_3_plot.tick_params(axis='y', labelsize=32)
    # Test_3_plot.set_title("Test 3", fontsize = 32)

    # Cmap_max = max(np.abs(Test_Data[3].flatten()))
    # Cmap_min = 0
    # Test_4_plot.imshow(Test_Data[3], interpolation = 'bilinear', cmap = "plasma", extent = axes_limits, vmin = Cmap_min, vmax = Cmap_max)

    # Test_4_plot.set_xlabel(r'$\alpha_{rel}\,\,[\mu$as]', fontsize = 32)
    # Test_4_plot.set_ylabel(r'$\delta_{rel}\,\,[\mu$as]', fontsize = 32)
    # Test_4_plot.tick_params(axis='x', labelsize=32)
    # Test_4_plot.tick_params(axis='y', labelsize=32)
    # Test_4_plot.set_title("Test 4", fontsize = 32)

    # Cmap_max = max(np.abs(Test_Data[4].flatten()))
    # Cmap_min = 0
    # Test_5_plot.imshow(Test_Data[4], interpolation = 'bilinear', cmap = "plasma", extent = axes_limits, vmin = Cmap_min, vmax = Cmap_max)

    # Test_5_plot.set_xlabel(r'$\alpha_{rel}\,\,[\mu$as]', fontsize = 32)
    # Test_5_plot.set_ylabel(r'$\delta_{rel}\,\,[\mu$as]', fontsize = 32)
    # Test_5_plot.tick_params(axis='x', labelsize=32)
    # Test_5_plot.tick_params(axis='y', labelsize=32)
    # Test_5_plot.set_title("Test 5", fontsize = 32)

    plt.rcParams['text.usetex'] = True
    plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
    plt.show()
