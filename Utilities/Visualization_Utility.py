import numpy as np
from matplotlib import pyplot as plt

from Support_functions.Parsers import*
from Support_functions.Shadows import*
from Support_functions.Image_processing import*

PARSE_SIMLATION             = True
PLOT_EHTIM_RESULTS          = True
PLOT_VIDA_TEMPLATE          = True
EXPORT_EHTIM_INPUT_FROM_SIM = False
CENTER_EHTIM_PLOT           = False

ARRAY = "2017"

Units = Units_class()

Sim_path = "..\\GB\\Sim_Results\\"

if PARSE_SIMLATION:

    Sim_Parser_0 = Simulation_Parser(Sim_path + "Ray_Tracer_output\\Gauss_Bonnet_n0")
    Sim_Parser_1 = Simulation_Parser(Sim_path + "Ray_Tracer_output\\Gauss_Bonnet_n1")
    Sim_Parser_2 = Simulation_Parser(Sim_path + "Ray_Tracer_output\\Gauss_Bonnet_n2")
    Sim_Parser_3 = Simulation_Parser(Sim_path + "Ray_Tracer_output\\Gauss_Bonnet_n3")

    Obs_effective_distance = Units.M87_DISTANCE_GEOMETRICAL

    Intensity_0, NT_Flux_0, NT_Redshift_0, NT_Flux_Shifted_0, Metadata_0 = Sim_Parser_0.get_plottable_sim_data()
    Intensity_1, NT_Flux_1, NT_Redshift_1, NT_Flux_Shifted_1, Metadata_1 = Sim_Parser_1.get_plottable_sim_data()
    Intensity_2, NT_Flux_2, NT_Redshift_2, NT_Flux_Shifted_2, Metadata_2 = Sim_Parser_2.get_plottable_sim_data()
    Intensity_3, NT_Flux_3, NT_Redshift_3, NT_Flux_Shifted_3, Metadata_3 = Sim_Parser_3.get_plottable_sim_data()

    #=============== PLot the Simulated Image ===============#

    Main_Figure = plt.figure(figsize = (20, 8))

    # Set the X and Y axis limits, rescaling them for an observer, located at "Obs_effective_distance", rather than the simulation OBS_DISTANCE, and converto to micro AS 
    axes_limits = np.array([(limit) for limit in Metadata_0[2]]) * Sim_Parser_0.OBS_DISTANCE / Obs_effective_distance * Units.RAD_TO_MICRO_AS

    # The literature (for some reason) has the X axis going positive to negative, 
    # so I invert the X axis limits
    axes_limits[0] = -axes_limits[0]
    axes_limits[1] = -axes_limits[1]
    
    Subplot      = Main_Figure.add_subplot(121)
    Data_to_plot = Intensity_0 + Intensity_1 + Intensity_2 + Intensity_3

    if EXPORT_EHTIM_INPUT_FROM_SIM:

        Sim_Parser_0.export_ehtim_data(Spacetime = Sim_Parser_0.metric, data = Data_to_plot, path = Sim_path)

    # Create the plot of the Simulated Image
    # I first convert the specific intensity to brightness temperature, so I can make a colorbar
    
    Data_to_plot = Units.Spectral_density_to_T(Data_to_plot / Units.W_M2_TO_JY, Sim_Parser_0.OBS_FREQUENCY) / Units.GIGA
    Image_norm   = max(Data_to_plot.flatten())

    Sim_subplot = Subplot.imshow(Data_to_plot, interpolation = 'bilinear', cmap = 'hot', extent = axes_limits, vmin = 0, vmax = Image_norm)

    colorbar = Main_Figure.colorbar(Sim_subplot, ax = Subplot, fraction=0.046, pad=0.04)
    colorbar.set_label(r"Brightness Temperature [$10^9$K]")

    Subplot.set_title("Simulated Image")
    Subplot.set_xlabel(r'$\alpha_{rel}\,\,[\mu$as]')
    Subplot.set_ylabel(r'$\delta_{rel}\,\,[\mu$as]')

    #=============== PLot the Brigtness Temperature at y = 0 of the Simulated Image ===============#

    Subplot = Main_Figure.add_subplot(122)

    # Convert the spectral density at y = 0 to brightness temperature, normalized to 10^9 Kelvin
    T_Brightness = Data_to_plot[int(Sim_Parser_0.Y_PIXEL_COUNT / 2)]
    T_Brightness_norm = max(T_Brightness)

    # Rescale the celestial coordinates for an observer, located at "Obs_effective_distance", rather than the simulation OBS_DISTANCE, and converto to micro AS
    x_coords  = np.linspace(Sim_Parser_0.WINDOW_LIMITS[0], Sim_Parser_0.WINDOW_LIMITS[1], Sim_Parser_0.X_PIXEL_COUNT) # These limits are in radians, for an observer located at the ray-tracer's OBS_DISTANCE
    x_coords *= Sim_Parser_0.OBS_DISTANCE / Obs_effective_distance * Units.RAD_TO_MICRO_AS

    # Set the aspect ratio of the figure to 1:1 (y:x)
    Subplot.set_aspect(2 * x_coords[-1] / T_Brightness_norm )

    # Create the plot of "T_b(alpha) | y = 0"
    Subplot.plot(-x_coords, T_Brightness)
    Subplot.invert_xaxis()
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

    # Main_Figure.suptitle(Fig_title)

    print(Fig_title)

if PLOT_EHTIM_RESULTS:

    #========================= EHTIM Parsing/Plotting =========================#

    Ehtim_Parser_no_blur = ehtim_Parser(Sim_path + "Ehtim_" + ARRAY + "\\Results_no_blur")
    Intensity_ehtim_no_blur, Ehtim_metadata_no_blur = Ehtim_Parser_no_blur.get_plottable_ehtim_data()

    Ehtim_Parser_blur = ehtim_Parser(Sim_path + "Ehtim_" + ARRAY + "\\Results_5")
    Intensity_ehtim_blur, Ehtim_metadata_blur = Ehtim_Parser_blur.get_plottable_ehtim_data()

    #========================= Plot the main EHTIM image =========================#

    # EHTIM saves the axis limits in arcsec - here I convert to micro-arcsec
    axes_limits = np.array([(limit) for limit in Ehtim_metadata_blur]) * Units.MEGA

    # The literature (for some reason) has the X axis going positive to negative, 
    # so I invert the X axis limits
    axes_limits[0] = -axes_limits[0]
    axes_limits[1] = -axes_limits[1]

    pixel_size = np.abs(axes_limits[0] - axes_limits[1]) / Ehtim_Parser_no_blur.X_PIXEL_COUNT / Units.MEGA * Units.ARCSEC_TO_RAD 

    Intensity_ehtim_no_blur = Units.Spectral_density_to_T(Intensity_ehtim_no_blur / pixel_size**2 / Units.W_M2_TO_JY, 230e9) / Units.GIGA
    Intensity_ehtim_blur    = Units.Spectral_density_to_T(Intensity_ehtim_blur    / pixel_size**2 / Units.W_M2_TO_JY, 230e9) / Units.GIGA

    # Create the plot of the Simulated Observations

    Ehtim_figure = plt.figure(figsize = (20, 7.5))

    Subplot = Ehtim_figure.add_subplot(121)
    pre_Convolution = Subplot.imshow(Intensity_ehtim_no_blur, interpolation = 'bilinear', cmap = 'hot', extent = axes_limits)
    Subplot.set_title("Simulated Observation pre- Clean Beam Convolution")
    Subplot.set_xlabel(r'$\alpha_{rel}\,\,[\mu$as]')
    Subplot.set_ylabel(r'$\delta_{rel}\,\,[\mu$as]')
    colorbar = Ehtim_figure.colorbar(pre_Convolution, ax = Subplot, fraction=0.046, pad=0.04)
    colorbar.set_label(r"Brightness Temperature [$10^9$K]")

    Subplot = Ehtim_figure.add_subplot(122)
    post_Convolution = Subplot.imshow(Intensity_ehtim_blur, interpolation = 'bilinear', cmap = 'hot', extent = axes_limits)
    Subplot.set_title("Simulated Observation post- Clean Beam Convolution")
    Subplot.set_xlabel(r'$\alpha_{rel}\,\,[\mu$as]')
    Subplot.set_ylabel(r'$\delta_{rel}\,\,[\mu$as]')
    colorbar = Ehtim_figure.colorbar(post_Convolution, ax = Subplot, fraction=0.046, pad=0.04)
    colorbar.set_label(r"Brightness Temperature [$10^9$K]")

    if PLOT_VIDA_TEMPLATE:

        # I then invert them back to not break the VIDA functions
        axes_limits[0] = -axes_limits[0]
        axes_limits[1] = -axes_limits[1]

        Ehtim_Vida_plot = plot_VIDA_templte(Intensity_ehtim_blur, axes_limits, Sim_path + "Ehtim_" + ARRAY + "\\fit_params", CENTER_EHTIM_PLOT, [50, 50, 50, 50])
        
Main_Figure.savefig(Sim_path + "Ehtim_" + ARRAY + "\\Sim_plot.png", bbox_inches = 'tight')
Ehtim_figure.savefig(Sim_path + "Ehtim_" + ARRAY + "\\Ehtim_plot.png", bbox_inches = 'tight')
Ehtim_Vida_plot.savefig(Sim_path + "Ehtim_" + ARRAY + "\\Ehtim_Vida_plot.png", bbox_inches = 'tight')

plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
plt.show()

