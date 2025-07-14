from Support_functions.Parsers import Units_class, Simulation_Parser
import matplotlib.pyplot as plt
import numpy as np
import os 
 
parent_directory = os.path.abspath('...')

Units = Units_class()

metric = "Kerr"

Hotspot_figure = plt.figure()
Subplot_position = Hotspot_figure.add_subplot(151)
Intensity_plot = Hotspot_figure.add_subplot(152)
QU_plot = Hotspot_figure.add_subplot(153)
Rel_U_plot = Hotspot_figure.add_subplot(154)
EVPA_plot = Hotspot_figure.add_subplot(155)

Time = []
Total_flux_w_rings = []
Polarized_flux_w_rings = []
Total_flux_without_rings = []
LP_Fraction = []
EVPA = []

Total_intensity = np.zeros((128, 128))

Net_Q = []
Net_U = []

Rel_Q = []
Rel_U = []

axes_limits = 0

for hotspot_azimuth_offset in range(0, 40):
            
    try:    
        
        """ Evaluate the simulataion results """
        Sim_parser_n0 = Simulation_Parser(parent_directory + "Reference_simulations\\Hotspot_Reference_Simulation_{}".format(hotspot_azimuth_offset) + "\\" + metric + "_n0")
        I_Intensity_0, Q_Intensity_n0, U_Intensity_n0, _, _, _ = Sim_parser_n0.get_plottable_sim_data()
        Total_flux_n0 = Sim_parser_n0.get_total_flux(Units.SGRA_DISTANCE_GEOMETRICAL, unit = "mJy")
        Polarized_flux_n0 = Sim_parser_n0.get_polarized_flux(Units.SGRA_DISTANCE_GEOMETRICAL, unit = "mJy")
                
        Sim_parser_n1 = Simulation_Parser(parent_directory + "Reference_simulations\\Hotspot_Reference_Simulation_{}".format(hotspot_azimuth_offset) + "\\" + metric + "_n1")
        I_Intensity_1, Q_Intensity_n1, U_Intensity_n1, _, _, _ = Sim_parser_n1.get_plottable_sim_data()
        Total_flux_n1 = Sim_parser_n1.get_total_flux(Units.SGRA_DISTANCE_GEOMETRICAL, unit = "mJy")
        Polarized_flux_n1 = Sim_parser_n1.get_polarized_flux(Units.SGRA_DISTANCE_GEOMETRICAL, unit = "mJy")
        
        Sim_parser_n2 = Simulation_Parser(parent_directory + "Reference_simulations\\Hotspot_Reference_Simulation_{}".format(hotspot_azimuth_offset) + "\\" + metric + "_n2")
        I_Intensity_2, Q_Intensity_n2, U_Intensity_n2, _, _, _ = Sim_parser_n2.get_plottable_sim_data()
        Total_flux_n2 = Sim_parser_n2.get_total_flux(Units.SGRA_DISTANCE_GEOMETRICAL, unit = "mJy")
        Polarized_flux_n2 = Sim_parser_n2.get_polarized_flux(Units.SGRA_DISTANCE_GEOMETRICAL, unit = "mJy")
                
        Sim_parser_n3 = Simulation_Parser(parent_directory + "Reference_simulations\\Hotspot_Reference_Simulation_{}".format(hotspot_azimuth_offset) + "\\" + metric + "_n3")
        I_Intensity_3, Q_Intensity_n3, U_Intensity_n3, _, _, _ = Sim_parser_n3.get_plottable_sim_data()
        Total_flux_n3 = Sim_parser_n3.get_total_flux(Units.SGRA_DISTANCE_GEOMETRICAL, unit = "mJy")
        Polarized_flux_n3 = Sim_parser_n3.get_polarized_flux(Units.SGRA_DISTANCE_GEOMETRICAL, unit = "mJy")
        
        Total_intensity = Total_intensity + I_Intensity_0 + I_Intensity_1 + I_Intensity_2 + I_Intensity_3
  
        Total_flux_w_rings.append(Total_flux_n0 + Total_flux_n1 + Total_flux_n2 + Total_flux_n3)
        Polarized_flux_w_rings.append(Polarized_flux_n0)
        Total_flux_without_rings.append(Total_flux_n0)
        Time.append(float(Sim_parser_n0.Simulation_metadata["Observervation Time [M]"]))
        
        axes_limits = Sim_parser_n0.Simulation_metadata["Observation Window Dimentions (-X,+X,-Y,+Y) [M]"].split(",")
        axes_limits = np.array([float(Limit) for Limit in axes_limits]) / Units.SGRA_DISTANCE_GEOMETRICAL
        axes_limits = np.tan(axes_limits) * Units.RAD_TO_MICRO_AS
        
        # The literature (for some reason) has the X axis going positive to negative, 
        # so I invert the X axis limits
        axes_limits[0] = -axes_limits[0]
        axes_limits[1] = -axes_limits[1]
        
        X_resolution: int = int(Sim_parser_n0.Simulation_metadata["Simulation Resolutoin"].split(" ")[0])
        Y_resolution: int = int(Sim_parser_n0.Simulation_metadata["Simulation Resolutoin"].split(" ")[2])
        
        Net_Q.append(sum((Q_Intensity_n0).flatten()) * abs(axes_limits[0] - axes_limits[1]) * abs(axes_limits[2] - axes_limits[3]) / X_resolution / Y_resolution / Units.RAD_TO_MICRO_AS**2 * 1e3)
        Net_U.append(sum((U_Intensity_n0).flatten()) * abs(axes_limits[0] - axes_limits[1]) * abs(axes_limits[2] - axes_limits[3]) / X_resolution / Y_resolution / Units.RAD_TO_MICRO_AS**2 * 1e3)
  
        Rel_Q.append(sum((Q_Intensity_n0).flatten()) / sum((I_Intensity_0).flatten()))
        Rel_U.append(sum((U_Intensity_n0).flatten()) / sum((I_Intensity_0).flatten()))
        
        LP_Fraction.append(Polarized_flux_n0 / Total_flux_n0)
        EVPA.append(0.5 * np.arctan2(sum((U_Intensity_n0).flatten()), sum((Q_Intensity_n0).flatten())))
                     
    except:
        print("Could not read simulation {}.".format(hotspot_azimuth_offset))
    
Subplot_position.imshow(Total_intensity, interpolation = 'bilinear', cmap = "hot", extent = axes_limits, vmin = 0, vmax = max(np.abs(Total_intensity.flatten()))) # type: ignore
Subplot_position.set_xlabel(r'$\alpha_{rel}\,\,[\mu$as]', fontsize = 24)
Subplot_position.set_ylabel(r'$\delta_{rel}\,\,[\mu$as]', fontsize = 24)

Time = np.array(Time) * 30 / 85
Time = Time - Time[0]

Intensity_plot.plot(Total_flux_without_rings, "r--")
# Intensity_plot.plot(Time, Total_flux_w_rings, "b-")
Intensity_plot.plot(Polarized_flux_w_rings, "k-")

Intensity_plot.set_ylabel(r"$Flux [mJy]$", fontsize = 24)
Intensity_plot.set_xlabel(r"Time [min]", fontsize = 24)
LP_Fraction_plot = Intensity_plot.twinx()
LP_Fraction_plot.plot(LP_Fraction) # type: ignore
LP_Fraction_plot.set_ylim(0, 1)
# Intensity_plot.set_ylim(0,25)

QU_plot.plot(0, 0, "ko")
QU_plot.plot(Net_Q, Net_U)
QU_plot.scatter(Net_Q[0], Net_U[0])
# QU_plot.set_ylim(-25,25)
# QU_plot.set_xlim(-25,25)

Rel_U_plot.plot(Rel_Q, Rel_U)
Rel_U_plot.plot(0, 0, "ko")
# Rel_U_plot.set_ylim(-0.7,0.7)
# Rel_U_plot.set_xlim(-0.7,0.7)

EVPA_plot.plot(EVPA)

plt.show()