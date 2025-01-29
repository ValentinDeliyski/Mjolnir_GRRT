from Support_functions.Parsers import Units_class, Simulation_Parser
import matplotlib.pyplot as plt
import numpy as np
import os 
 
parent_directory = os.path.abspath('...')

Units = Units_class()

Total_intensity = np.zeros((512,512))
Total_flux = []
Azimuth = []
axes_limits = 0

for hotspot_azimuth_offset in range(50):
            
    """ Evaluate the simulataion results """
    Sim_parser_n0 = Simulation_Parser(parent_directory + "Reference_simulations\\Hotspot_Reference_Simulation_{}".format(hotspot_azimuth_offset) + "\\Kerr_n0")
    I_Intensity_0, _, _, _, _, _, _ = Sim_parser_n0.get_plottable_sim_data()
    Total_flux_n0 = Sim_parser_n0.get_total_flux(Units.SGRA_DISTANCE_GEOMETRICAL, unit = "mJy")
            
    Sim_parser_n1 = Simulation_Parser(parent_directory + "Reference_simulations\\Hotspot_Reference_Simulation_{}".format(hotspot_azimuth_offset) + "\\Kerr_n1")
    I_Intensity_1, _, _, _, _, _, _ = Sim_parser_n1.get_plottable_sim_data()
    Total_flux_n1 = Sim_parser_n1.get_total_flux(Units.SGRA_DISTANCE_GEOMETRICAL, unit = "mJy")
            
    Sim_parser_n2 = Simulation_Parser(parent_directory + "Reference_simulations\\Hotspot_Reference_Simulation_{}".format(hotspot_azimuth_offset) + "\\Kerr_n2")
    I_Intensity_2, _, _, _, _, _, _ = Sim_parser_n2.get_plottable_sim_data()
    Total_flux_n2 = Sim_parser_n2.get_total_flux(Units.SGRA_DISTANCE_GEOMETRICAL, unit = "mJy")
            
    Sim_parser_n3 = Simulation_Parser(parent_directory + "Reference_simulations\\Hotspot_Reference_Simulation_{}".format(hotspot_azimuth_offset) + "\\Kerr_n3")
    I_Intensity_3, _, _, _, _, _, _ = Sim_parser_n3.get_plottable_sim_data()
    Total_flux_n3 = Sim_parser_n3.get_total_flux(Units.SGRA_DISTANCE_GEOMETRICAL, unit = "mJy")
    
    Total_intensity = Total_intensity +  I_Intensity_0 + I_Intensity_1 + I_Intensity_2 + I_Intensity_3        
    Total_flux.append(Total_flux_n0 + Total_flux_n1 + Total_flux_n2 + Total_flux_n3)
    Azimuth.append(Sim_parser_n3.hotspot_azimuth)
    
    axes_limits = np.array([(limit) for limit in Sim_parser_n3.WINDOW_LIMITS]) / Units.SGRA_DISTANCE_GEOMETRICAL
    axes_limits = np.tan(axes_limits) * Units.RAD_TO_MICRO_AS

    # The literature (for some reason) has the X axis going positive to negative, 
    # so I invert the X axis limits
    axes_limits[0] = -axes_limits[0]
    axes_limits[1] = -axes_limits[1]
    
Hotspot_figure = plt.figure()

Subplot_position = Hotspot_figure.add_subplot(121)
Subplot_position.imshow(Total_intensity, interpolation = 'bilinear', cmap = "hot", extent = axes_limits, vmin = 0, vmax = max(np.abs(Total_intensity.flatten())))
Subplot_position.set_xlabel(r'$\alpha_{rel}\,\,[\mu$as]',fontsize = 24)
Subplot_position.set_ylabel(r'$\delta_{rel}\,\,[\mu$as]',fontsize = 24)

Total_flux = np.array(Total_flux)
Azimuth    = np.array(Azimuth)

Azimuth[Azimuth < 90] = Azimuth[Azimuth < 90] + 360

idx        = np.argsort(Azimuth)
Azimuth    = Azimuth[idx]
Total_flux = Total_flux[idx]

# Azimuth = Azimuth - 55

Intensity_plot = Hotspot_figure.add_subplot(122)
Intensity_plot.plot(Azimuth, Total_flux / 15)

Intensity_plot.set_ylim([0, 0.6])
Intensity_plot.set_ylabel(r"$F_{\text{hotspot}}/F_{\text{S2}}$",fontsize = 24)
Intensity_plot.set_xlabel(r"Azimuth [deg]",fontsize = 24)

plt.show()