import numpy as np
from matplotlib import pyplot as plt

from Support_functions.Parsers import*
from Support_functions.Shadows import*

Units        = Units_class()
Sim_Parser_0 = Simulation_Parser("..\\Sim_results\\Kerr_n0")
Sim_Parser_1 = Simulation_Parser("..\\Sim_results\\Kerr_n1")
Sim_Parser_2 = Simulation_Parser("..\\Sim_results\\Kerr_n2")
Sim_Parser_3 = Simulation_Parser("..\\Sim_results\\Kerr_n3")

Obs_effective_distance = Units.M87_DISTANCE_GEOMETRICAL

Intensity_0, NT_Flux_0, NT_Redshift_0, NT_Flux_Shifted_0, Metadata_0 = Sim_Parser_0.get_plottable_sim_data()
Intensity_1, NT_Flux_1, NT_Redshift_1, NT_Flux_Shifted_1, Metadata_1 = Sim_Parser_1.get_plottable_sim_data()
Intensity_2, NT_Flux_2, NT_Redshift_2, NT_Flux_Shifted_2, Metadata_2 = Sim_Parser_2.get_plottable_sim_data()
Intensity_3, NT_Flux_3, NT_Redshift_3, NT_Flux_Shifted_3, Metadata_3 = Sim_Parser_3.get_plottable_sim_data()

#=============== PLot the Simulated Image ===============#

Main_Figure = plt.figure()

# Set the X and Y axis limits, rescaling them for an observer, located at "Obs_effective_distance", rather than the simulation OBS_DISTANCE, and converto to micro AS 
axes_limits = np.array([(limit) for limit in Metadata_0[2]]) * Sim_Parser_0.OBS_DISTANCE / Obs_effective_distance * Units.RAD_TO_MICRO_AS

Subplot      = Main_Figure.add_subplot(121)
Data_to_plot = Intensity_0 + Intensity_1 + Intensity_2 + Intensity_3
Image_norm   = max(Data_to_plot.flatten())

# Create the plot of the Simulated Image
Subplot.imshow(Data_to_plot, interpolation = 'bilinear', cmap = 'hot', extent = axes_limits, vmin = 0, vmax = Image_norm)
Subplot.set_title("Simulated Image")
Subplot.set_xlabel(r'$\alpha_{rel}\,\,[\mu$as]')
Subplot.set_ylabel(r'$\delta_{rel}\,\,[\mu$as]')

#=============== PLot the Brigtness Temperature at y = 0 of the Simulated Image ===============#

Subplot = Main_Figure.add_subplot(122)

# Convert the spectral density at y = 0 to brightness temperature, normalized to GK
T_Brightness = Units.Spectral_density_to_T(Data_to_plot[int(Sim_Parser_0.Y_PIXEL_COUNT / 2) ] / Units.W_M2_TO_JY, Sim_Parser_0.OBS_FREQUENCY) / Units.GIGA
T_Brightness_norm = max(T_Brightness)

# Rescale the celestial coordinates for an observer, located at "Obs_effective_distance", rather than the simulation OBS_DISTANCE, and converto to micro AS
x_coords  = np.linspace(Sim_Parser_0.WINDOW_LIMITS[0], Sim_Parser_0.WINDOW_LIMITS[1], Sim_Parser_0.X_PIXEL_COUNT) # These limits are in radians, for an observer located at the ray-tracer's OBS_DISTANCE
x_coords *= Sim_Parser_0.OBS_DISTANCE / Obs_effective_distance * Units.RAD_TO_MICRO_AS

# Set the aspect ratio of the figure to 1:1 (y:x)
Subplot.set_aspect(2 * x_coords[-1] / T_Brightness_norm )

# Create the plot of "T_b(alpha) | y = 0"
Subplot.plot(x_coords, T_Brightness)
Subplot.set_ylim([0, 1.1 * T_Brightness_norm])
Subplot.set_title("Brightness temperature at " + r'$\delta_{rel} = 0$')
Subplot.set_xlabel(r'$\alpha_{rel}\,\,[\mu$as]')
Subplot.set_ylabel(r'$T_b\,\,[10^9\, K]$')

Total_flux = (Sim_Parser_0.get_total_flux(Obs_effective_distance) + 
              Sim_Parser_1.get_total_flux(Obs_effective_distance) + 
              Sim_Parser_2.get_total_flux(Obs_effective_distance) + 
              Sim_Parser_3.get_total_flux(Obs_effective_distance))

if Sim_Parser_0.emission_model == " Phenomenological":

    Fig_title = (Sim_Parser_0.metric + ", " + Sim_Parser_0.disk_profile + 
                 ", Height Scale [M] = {}".format(Sim_Parser_0.height_scale) + 
                 ", Radial Scale [M] = {}".format(Sim_Parser_0.radial_scale) + 
                 ", Total Flux [Jy] = {}".format(np.round(Total_flux,4)) +
                 ", Emission constant = {}".format(Sim_Parser_0.Emission_Scale))
else:

    Fig_title = (Sim_Parser_0.metric + ", " + Sim_Parser_0.disk_profile + 
                 ", Opening Angle [deg] = {}".format(np.round(np.arctan(Sim_Parser_0.disk_opening_angle) * Units.RAD_TO_DEG, 2)) + 
                 ", R_0 [M] = {}".format(Sim_Parser_0.R_0) + 
                 ", Cutoff [M] = {}".format(Sim_Parser_0.R_Cutoff) + 
                 ", Total Flux [Jy] = {}".format(np.round(Total_flux,4)))


Main_Figure.suptitle(Fig_title)

# Sim_Parser_0.export_ehtim_data(data = Data_to_plot)

plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
plt.show()

