import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from astropy.io import fits
import os

from Support_functions.Parsers import Simulation_Parser, Units_class, ehtim_Parser, VIDA_params_Parser
from Support_functions.Image_processing import*

class Sim_Visualizer():

    def __init__(self, Sim_path: str, 
                 Sim_Frequency_Bins: list[str], 
                 Array: list[str], 
                 Font_size: int, 
                 Label_Pad: int,
                 Common_file_name: str,
                 Respect_folder_structure: bool):

        self.Sim_Parsers     = []
        self.Ehtim_Parsers   = []
        self.VIDA_Parsers    = []
        self.Metric          = Common_file_name
        self.Sim_path        = Sim_path
        self.Arrays          = Array
        self.Units           = Units_class()
        self.Frequency_Bins  = Sim_Frequency_Bins
        self.Total_flux_str  = ""
        self.Console_log_str = []
        self.Respect_folder_structure = Respect_folder_structure
        self.Font_size = Font_size
        self.Label_Pad = Label_Pad
        #========= Enums =========#

        self.NO_BLUR = 0
        self.BLUR    = 1

        #=========================#

        self.__make_paths(Sim_path)

        for Sim_number, _ in enumerate(self.Frequency_Bins):

            try:
                Sim_Parser_0 = Simulation_Parser(self.Ray_tracer_paths[Sim_number] + "_n0")
                Sim_Parser_1 = Simulation_Parser(self.Ray_tracer_paths[Sim_number] + "_n1")
                Sim_Parser_2 = Simulation_Parser(self.Ray_tracer_paths[Sim_number] + "_n2")
                Sim_Parser_3 = Simulation_Parser(self.Ray_tracer_paths[Sim_number] + "_n3")

                self.Sim_Parsers.append([Sim_Parser_0, Sim_Parser_1, Sim_Parser_2, Sim_Parser_3])
                
            except:
                print("Could not parse ray-tracer logs!")
                print("I looked at this path: {}".format(self.Ray_tracer_paths[Sim_number]))
                exit()

            if int(Sim_Parser_0.Simulation_metadata["Active Simulation Mode"]) != 2:

                Total_flux = (self.Sim_Parsers[Sim_number][0].get_total_flux(self.Units.M87_DISTANCE_GEOMETRICAL, unit = "mJy") +
                              self.Sim_Parsers[Sim_number][1].get_total_flux(self.Units.M87_DISTANCE_GEOMETRICAL, unit = "mJy") +
                              self.Sim_Parsers[Sim_number][2].get_total_flux(self.Units.M87_DISTANCE_GEOMETRICAL, unit = "mJy") +
                              self.Sim_Parsers[Sim_number][3].get_total_flux(self.Units.M87_DISTANCE_GEOMETRICAL, unit = "mJy"))
                
                self.Total_flux_str = self.Total_flux_str + "Total flux at {}GHz = {} [mJy]\n".format(float(self.Sim_Parsers[Sim_number][0].Simulation_metadata["Observation Frequency [Hz]"]) / 1e9, np.round(Total_flux, 4))

        if int(Sim_Parser_0.Simulation_metadata["Active Simulation Mode"]) != 2: # type: ignore

            for msg in self.Total_flux_str.split("\n"):
                
                print("=" * len(msg))
                print(msg)

        for Array_num, _ in enumerate(self.Ehtim_paths):
            try:
                Ehtim_Parser_no_blur = ehtim_Parser(self.Ehtim_paths[Array_num] + "Results") 
                Ehtim_Parser_blur    = ehtim_Parser(self.Ehtim_paths[Array_num] + "Results_blur") 

                self.Ehtim_Parsers.append([Ehtim_Parser_no_blur, Ehtim_Parser_blur])
                self.Console_log_str.append(self.Total_flux_str)

            except:
                print("Could not parse ehtim logs!")
                print("I looked at this path: {}".format(self.Ehtim_paths[Array_num] + "Results"))

            try:
                VIDA_parser = VIDA_params_Parser(self.Ehtim_paths[Array_num] + "fit_params")
                self.VIDA_Parsers.append(VIDA_parser)

            except:
                print("Could not parse VIDA template!")
                print("I looked at this path: {}".format(self.Ehtim_paths[Array_num] + "fit_params"))

    def __make_paths(self, Sim_path: str):

        self.Ray_tracer_paths = []
        self.Ehtim_paths      = []
        
        if self.Respect_folder_structure:

            for freq in self.Frequency_Bins:

                self.Ray_tracer_paths.append(Sim_path + freq + "GHz\\" + "Sim_Results\\Ray_tracer_output\\" + self.Metric)

        else:

            for freq in self.Frequency_Bins:
                
                self.Ray_tracer_paths.append(Sim_path)  

        for array in self.Arrays:
            for freq in self.Frequency_Bins:
                
                self.Ehtim_paths.append(Sim_path + freq + "GHz\\" + "Sim_Results\\Ehtim_" + array + "\\")
        
    def plot_ray_tracer_results(self, 
                                Export_data_for_Ehtim: bool, 
                                Radiation_Component: str,
                                Save_Figures: bool,
                                Custom_fig_title: str):

        Obs_effective_distance = self.Units.M87_DISTANCE_GEOMETRICAL
        Frequency_str_addon    = ""

        if len(self.Frequency_Bins) == 1:

            Main_Figure = plt.figure(figsize = (20, 8))
            
        else:
            Main_Figure = plt.figure(figsize = (20, 16))

        Main_Figure.suptitle(Custom_fig_title, fontsize = self.Font_size)

        for Sim_number, Freq_str in enumerate(self.Frequency_Bins):          
            
            Obs_frequency: float = float(self.Sim_Parsers[Sim_number][0].Simulation_metadata["Observation Frequency [Hz]"])
            X_resolution: int = int(self.Sim_Parsers[Sim_number][0].Simulation_metadata["Simulation Resolutoin"].split(" ")[0])
            Y_resolution: int = int(self.Sim_Parsers[Sim_number][0].Simulation_metadata["Simulation Resolutoin"].split(" ")[2])

            I_Intensity_0, Q_Intensity_0, U_Intensity_0, V_Intensity_0, Disk_redshift_n0, Disk_flux_n0 = self.Sim_Parsers[Sim_number][0].get_plottable_sim_data()
            I_Intensity_1, Q_Intensity_1, U_Intensity_1, V_Intensity_1, Disk_redshift_n1, Disk_flux_n1 = self.Sim_Parsers[Sim_number][1].get_plottable_sim_data()
            I_Intensity_2, Q_Intensity_2, U_Intensity_2, V_Intensity_2, Disk_redshift_n2, Disk_flux_n2 = self.Sim_Parsers[Sim_number][2].get_plottable_sim_data()
            I_Intensity_3, Q_Intensity_3, U_Intensity_3, V_Intensity_3, Disk_redshift_n3, Disk_flux_n3 = self.Sim_Parsers[Sim_number][3].get_plottable_sim_data()

            #=============== PLot the Simulated Image ===============#

            Fig_title = "Simulated Image at {}GHz".format(int(Obs_frequency / 1e9))
            X_Slice_tile = "Brightness temperature at " + r'$\delta_{\text{rel}} = 0$'
            X_Slice_y_label = r'$T_b\,\,[10^9\, K]$'

            # Set the X and Y axis limits, rescaling them for an observer, located at "Obs_effective_distance", rather than the simulation "Observer Distance [M]", and conver to to micro AS 
            axes_limits = self.Sim_Parsers[Sim_number][0].Simulation_metadata["Observation Window Dimentions (-X,+X,-Y,+Y) [M]"].split(",")
            axes_limits = np.array([float(Limit) for Limit in axes_limits]) / Obs_effective_distance
            axes_limits = np.arctan(axes_limits) * self.Units.RAD_TO_MICRO_AS

            # The literature (for some reason) has the X axis going positive to negative, 
            # so I invert the X axis limits
            axes_limits[0] = -axes_limits[0]
            axes_limits[1] = -axes_limits[1]
            
            Subplot_count = 100 * len(self.Frequency_Bins)

            Subplot  = Main_Figure.add_subplot(Subplot_count + 20 + (2 * Sim_number + 1))
            Colormap = "seismic"

            Total_intensity = I_Intensity_0 + I_Intensity_1 + I_Intensity_2 + I_Intensity_3
            Intensity_normalization = max(np.abs(Total_intensity.flatten()))

            if Radiation_Component == "Stokes I":
                Data_to_plot_Intensity = I_Intensity_0 + I_Intensity_1 + I_Intensity_2 + I_Intensity_3
                Data_to_plot = self.Units.Spectral_density_to_T(Data_to_plot_Intensity / self.Units.W_M2_TO_JY, 
                                                                Obs_frequency) / self.Units.GIGA
                                                                
                Cbar_label = r"Brightness Temperature [$10^9$K]"
                Colormap = "hot"

                Cmap_max = max(np.abs(Data_to_plot.flatten()))
                Cmap_min = 0

            elif Radiation_Component == "Stokes Q":
                
                Data_to_plot = (Q_Intensity_0 + Q_Intensity_1 + Q_Intensity_2 + Q_Intensity_3) / Intensity_normalization * 100
                Cbar_label = r"Q Fractional Intensity [\%]"
                Cmap_max = max(np.abs(Data_to_plot.flatten())) 
                Cmap_min = -Cmap_max

            elif Radiation_Component == "Stokes U":
                Data_to_plot = (U_Intensity_0 + + U_Intensity_1 + U_Intensity_2 + U_Intensity_3) / Intensity_normalization * 100
                Cbar_label = r"U Fractional Intensity [\%]"
                Cmap_max = max(np.abs(Data_to_plot.flatten()))
                Cmap_min = -Cmap_max

            elif Radiation_Component == "Stokes V":
                Data_to_plot = (V_Intensity_0 + V_Intensity_1 + V_Intensity_2 + V_Intensity_3) / Intensity_normalization * 100
                Cbar_label = r"V Fractional Intensity [\%]"
                Cmap_max = max(np.abs(Data_to_plot.flatten()))
                Cmap_min = -Cmap_max

            elif Radiation_Component == "LP Fraction":
                
                U_intensity = U_Intensity_0 + U_Intensity_1 + U_Intensity_2 + U_Intensity_3
                Q_intensity = Q_Intensity_0 + Q_Intensity_1 + Q_Intensity_2 + Q_Intensity_3
                
                Cbar_label   = r"LP fraction [\%]"
                
                Data_to_plot = sqrt(U_intensity**2 + Q_intensity**2) / Intensity_normalization * 100
                # Data_to_plot[Total_intensity < Intensity_normalization / 100] = 0
                Cmap_max = max(Data_to_plot.flatten())
                Cmap_min = 0
                
            elif Radiation_Component == "NT":
                
                Shifted_flux_n0: NDArray = Disk_flux_n0 * (Disk_redshift_n0)**4
                Shifted_flux_n1: NDArray = Disk_flux_n1 * (Disk_redshift_n1)**4
                Shifted_flux_n2: NDArray = Disk_flux_n2 * (Disk_redshift_n2)**4
                Shifted_flux_n3: NDArray = Disk_flux_n3 * (Disk_redshift_n3)**4
            
                Data_to_plot = Shifted_flux_n0
                Data_to_plot = Data_to_plot + Shifted_flux_n1 * (Data_to_plot == 0)
                Data_to_plot = Data_to_plot + Shifted_flux_n2 * (Data_to_plot == 0)
                Data_to_plot = Data_to_plot + Shifted_flux_n3 * (Data_to_plot == 0)
                Data_to_plot = Data_to_plot / 1e-6

                Cmap_max = max(np.abs(Data_to_plot.flatten()))
                Cmap_min = 0

                Colormap = "hot"

                Subplot.set_aspect(X_resolution / Y_resolution)

                Cbar_label   = r"Intensity [$10^{-6}\dot{M}M^{-2}$]"
                X_Slice_tile = r"Intensity at $\delta_{\text{rel}} = 0$"

                Fig_title       = "Simulated Image"
                X_Slice_y_label = r"Intensity at $\delta_{\text{rel}} = 0$ $[10^{-6}\dot{M}M^{-2}]$"

            else:
                print("Incorrect Radiation Component!")
                return

            if Export_data_for_Ehtim:
                self.Sim_Parsers[Sim_number][0].export_ehtim_data(Spacetime = self.Sim_Parsers[Sim_number][0].metric, 
                                                                  data = Data_to_plot_Intensity,  # type: ignore
                                                                  path = self.Sim_path)

            # Create the plot of the Simulated Image
            Sim_subplot = Subplot.imshow(Data_to_plot, interpolation = 'bilinear', cmap = Colormap, extent = tuple(axes_limits), vmin = Cmap_min, vmax = Cmap_max)

            colorbar = Main_Figure.colorbar(Sim_subplot, ax = Subplot, fraction = 0.046, pad = 0.04)
            colorbar.set_label(Cbar_label, fontsize = self.Font_size, labelpad = self.Label_Pad)
            colorbar.ax.tick_params(labelsize = self.Font_size)

            Subplot.set_title(Fig_title, fontsize = self.Font_size)
            Subplot.set_xlabel(r'$\alpha_{rel}\,\,[\mu$as]', fontsize = self.Font_size)
            Subplot.set_ylabel(r'$\delta_{rel}\,\,[\mu$as]', fontsize = self.Font_size)

            plt.xticks(fontsize = self.Font_size)
            plt.yticks(fontsize = self.Font_size)

            #=============== PLot the Brigtness Temperature at y = 0 of the Simulated Image ===============#

            Subplot = Main_Figure.add_subplot(Subplot_count + 20 + (2 * Sim_number + 2))

            # Convert the spectral density at y = 0 to brightness temperature, normalized to 10^9 Kelvin
            T_Brightness = Data_to_plot[int(X_resolution / 2)]
            T_Brightness_norm     = max(T_Brightness)
            T_Brightness_min_norm = min(T_Brightness)
            x_coords = np.linspace(axes_limits[0], axes_limits[1], X_resolution)

            # Set the aspect ratio of the figure to 1:1 (y:x)
            # Subplot.set_aspect(2 * x_coords[-1] / T_Brightness_norm)

            # Create the plot of "T_b(alpha) | y = 0"
            Subplot.plot(x_coords, T_Brightness)
            Subplot.invert_xaxis()
            Subplot.set_ylim(1.1 * T_Brightness_min_norm, 1.1 * T_Brightness_norm)
            Subplot.set_title(X_Slice_tile, fontsize = self.Font_size)
            Subplot.set_xlabel(r'$\alpha_{rel}\,\,[\mu$as]', fontsize = self.Font_size)
            Subplot.set_ylabel(X_Slice_y_label, fontsize = self.Font_size, labelpad = self.Label_Pad)

            Frequency_str_addon += Freq_str + "_"

            plt.xticks(fontsize = self.Font_size)
            plt.yticks(fontsize = self.Font_size)

        Frequency_str_addon = Frequency_str_addon[:len(Frequency_str_addon) - 1]

        Main_Figure.tight_layout()

        if Save_Figures:

            if self.Respect_folder_structure:
                Figures_folder_path = self.Sim_path + "Figures\\"
            else:
                Figures_folder_path = self.Sim_path + "\\Figures\\"

            if not os.path.exists(Figures_folder_path):
                os.makedirs(Figures_folder_path)

            Main_Figure.savefig(Figures_folder_path + 
                                "Ray_tracer_plot_" + 
                                Frequency_str_addon +
                                "_" + Radiation_Component +
                                ".png", bbox_inches = 'tight')

    def plot_EHTIM_results(self, Make_contour_plots: bool, Contour_specs: list, Save_Figures: bool, Plot_no_blur: bool, Custom_fig_title: str):

        for Array_num, Array_str in enumerate(self.Arrays):

            if Plot_no_blur:
                if len(self.Frequency_Bins) == 1:            
                    Ehtim_figure = plt.figure(figsize = (25, 12))

                else:
                    Ehtim_figure = plt.figure(figsize = (25, 24))

            else:
                if len(self.Frequency_Bins) == 1:            
                    Ehtim_figure = plt.figure(figsize = (14, 13))

                else:
                    Ehtim_figure = plt.figure(figsize = (25, 12))


            for Sim_number, _ in enumerate(self.Frequency_Bins):

                Index = Sim_number + Array_num * len(self.Frequency_Bins)

                # The nested lists are getting out of hand, so Im abbreviating this here 
                Ehtim_Parser_no_Blur = self.Ehtim_Parsers[Index][self.NO_BLUR]
                Ehtim_Parser_Blur    = self.Ehtim_Parsers[Index][self.BLUR]

                #========================= EHTIM Parsing/Plotting =========================#

                Intensity_ehtim_no_blur_jy, _                = Ehtim_Parser_no_Blur.get_plottable_ehtim_data()
                Intensity_ehtim_blur_jy, Ehtim_metadata_blur = Ehtim_Parser_Blur.get_plottable_ehtim_data()

                #========================= Plot the main EHTIM image =========================#

                # EHTIM saves the axis limits in arcsec - here I convert to micro-arcsec
                axes_limits = np.array([(limit) for limit in Ehtim_metadata_blur]) * self.Units.MEGA

                # The literature (for some reason) has the X axis going positive to negative, 
                # so I invert the X axis limits
                axes_limits[0] = -axes_limits[0]
                axes_limits[1] = -axes_limits[1]

                pixel_size = np.abs(axes_limits[0] - axes_limits[1]) / Ehtim_Parser_no_Blur.X_PIXEL_COUNT / self.Units.MEGA * self.Units.ARCSEC_TO_RAD 

                Intensity_ehtim_no_blur_T = self.Units.Spectral_density_to_T(Intensity_ehtim_no_blur_jy / pixel_size**2 / self.Units.W_M2_TO_JY, Ehtim_Parser_no_Blur.OBS_FREQUENCY * self.Units.GIGA) / self.Units.GIGA
                Intensity_ehtim_blur_T    = self.Units.Spectral_density_to_T(Intensity_ehtim_blur_jy    / pixel_size**2 / self.Units.W_M2_TO_JY, Ehtim_Parser_Blur.OBS_FREQUENCY * self.Units.GIGA) / self.Units.GIGA

                # Create the plot of the Simulated Observations
                if Plot_no_blur:
                    Subplot = Ehtim_figure.add_subplot(len(self.Frequency_Bins) * 100 + 20 + (2 * Sim_number + 1))

                    Subplot.set_title("Pre-Clean Beam Convolution at {}GHz".format(int(Ehtim_Parser_no_Blur.OBS_FREQUENCY)), fontsize = self.Font_size)
                    Subplot.set_xlabel(r'$\alpha_{rel}\,\,[\mu$as]', fontsize = self.Font_size)
                    Subplot.set_ylabel(r'$\delta_{rel}\,\,[\mu$as]', fontsize = self.Font_size)

                    pre_Convolution_T = Subplot.imshow(Intensity_ehtim_no_blur_T,  interpolation = 'bilinear', cmap = 'hot', extent = tuple(axes_limits))
                    colorbar = Ehtim_figure.colorbar(pre_Convolution_T, ax = Subplot, fraction=0.046, pad=0.04)
                    colorbar.set_label(r"Brightness Temperature [$10^9$K]", fontsize = self.Font_size, labelpad = self.Label_Pad)
                    colorbar.ax.tick_params(labelsize = self.Font_size)

                    plt.xticks(fontsize = self.Font_size)
                    plt.yticks(fontsize = self.Font_size)
                    
                    Subplot = Ehtim_figure.add_subplot(len(self.Frequency_Bins) * 100 + 20 + (2 * Sim_number + 2))

                else:
                    Subplot = Ehtim_figure.add_subplot(100 + len(self.Frequency_Bins) * 10 + (Sim_number + 1))

                # Make the contour plots
                if Make_contour_plots:
                    self.plot_contours([Ehtim_Parser_Blur], self.VIDA_Parsers[Index], Subplot, Contour_specs[Sim_number])

                Subplot.set_title("Post Clean Beam Convolution at {}GHz".format(int(Ehtim_Parser_Blur.OBS_FREQUENCY)), fontsize = self.Font_size)
                Subplot.set_xlabel(r'$\alpha_{rel}\,\,[\mu$as]', fontsize = self.Font_size)
                Subplot.set_ylabel(r'$\delta_{rel}\,\,[\mu$as]', fontsize = self.Font_size)

                post_Convolution_T = Subplot.imshow(Intensity_ehtim_blur_T, interpolation = 'bilinear', cmap = 'hot', extent = tuple(axes_limits))
                colorbar = Ehtim_figure.colorbar(post_Convolution_T, ax = Subplot, fraction = 0.046, pad = 0.04)
                colorbar.set_label(r"Brightness Temperature [$10^9$K]", fontsize = self.Font_size, labelpad = self.Label_Pad)
                colorbar.ax.tick_params(labelsize = self.Font_size)
                
                plt.xticks(fontsize = self.Font_size)
                plt.yticks(fontsize = self.Font_size)

            if Custom_fig_title != None:

                if Array_str in ["2017", "2022", "2025"]:
                    Ehtim_figure.suptitle("{}, viewed by EHT {}".format(Custom_fig_title, Array_str), fontsize = 1.2 * self.Font_size) 

                else:
                    Ehtim_figure.suptitle("{}, viewed by {}".format(Custom_fig_title, Array_str), fontsize = 1.2 * self.Font_size)

            else:
                Ehtim_figure.suptitle("Using Array {}".format(Array_str), fontsize = 1.2 * self.Font_size)

            Ehtim_figure.tight_layout()

            if Save_Figures:

                if not os.path.exists(self.Sim_path + "Figures\\"):
                    os.makedirs(self.Sim_path + "Figures\\")

                fig_title = "Ehtim_plot_" + Array_str

                if Plot_no_blur:
                    fig_title += "_no_blur"

                if Make_contour_plots:
                    fig_title += "_contour"  

                if Array_str == "ngEHT" and len(self.Frequency_Bins) == 1:
                    fig_title += "_" + str(int(Ehtim_Parser_Blur.OBS_FREQUENCY)) # type: ignore

                fig_title += ".png"

                Ehtim_figure.savefig(self.Sim_path + 
                                    "Figures\\" + 
                                    fig_title, bbox_inches = 'tight')

    def plot_VIDA_style(self, 
                        Center_plot: bool, 
                        Save_Figures: bool,
                        Custom_fig_title: str):

        for Array_num, Array_str in enumerate(self.Arrays):
            
            for Sim_number, Freq_str in enumerate(self.Frequency_Bins):

                Index = Sim_number + Array_num * len(self.Frequency_Bins)
                    
                # The nested lists are getting out of hand, so Im abbreviating this here 
                Ehtim_Parser_Blur = self.Ehtim_Parsers[Index][self.BLUR]
                Vida_Parser = self.VIDA_Parsers[Index]

                Ehtim_Vida_plot, Brightness_ratio_str = self.plot_VIDA_templte(Ehtim_Parsers    = [Ehtim_Parser_Blur], 
                                                                               VIDA_parser      = Vida_Parser, 
                                                                               CROP             = Center_plot, 
                                                                               crop_rel_rage    = [50, 50, 50, 50],
                                                                               Plot_Brihtness_T = True,
                                                                               Custom_fig_title = Custom_fig_title,
                                                                               Array_str        = Array_str)
                
                Ehtim_Vida_plot.tight_layout()

                print(Brightness_ratio_str)

                self.Console_log_str[Array_num] += "\n" + Brightness_ratio_str

                if Save_Figures:

                    if not os.path.exists(self.Sim_path + "Figures\\"):
                        os.makedirs(self.Sim_path + "Figures\\")

                    Ehtim_Vida_plot.savefig(self.Sim_path + 
                                        "Figures\\" + 
                                        "Ehtim_Vida_plot_" + 
                                            Array_str + 
                                            "_" +
                                            Freq_str + 
                                            ".png", bbox_inches = 'tight')

    def create_EHTIM_superposition(self):
        
        if len(self.Frequency_Bins) < 2:

            print("Not enough frequency bins to make a superposition!")
            return

        else:

            Ehtim_outputs = []

            for Array_num, Array_str in enumerate(self.Arrays):
                
                for Sim_number, _ in enumerate(self.Frequency_Bins):

                    Ehtim_outputs.append(fits.open(self.Ehtim_paths[Sim_number + Array_num * (len(self.Arrays) - 1)] + "Results_blur.fits"))

                Superposition_File = Ehtim_outputs[0].__deepcopy__()

                for Sim_number, _ in enumerate(self.Frequency_Bins):

                    Superposition_File[0].data += Ehtim_outputs[Sim_number + Array_num * (len(self.Arrays) - 1)][0].data

                Superposition_File[0].data -= Ehtim_outputs[0][0].data

                Superposition_File.writeto("Superposition_{}.fits".format(Array_str), overwrite = True)

                print("Superposition of {} Array exported to fits file!".format(Array_str))

    def plot_superposition(self, Center_plot: bool, Save_Figures: bool, Custom_fig_title: str):

        for Array_num, Array_str in enumerate(self.Arrays):

            if len(self.Frequency_Bins) > 1:
                Ehtim_Parsers = [self.Ehtim_Parsers[Sim_num + Array_num * (len(self.Arrays) - 1)][self.BLUR] 
                                 for Sim_num, _ in enumerate(self.Frequency_Bins)]
                    
            else:
                print("Not enough frequency bins to make a superposition!")
                return
            
            try:
                VIDA_Parser = VIDA_params_Parser(self.Sim_path + "\\Superposition\\fit_params_superposition")

            except:
                print("Could not parse VIDA template!")
                print("I looked at this path: {}".format(self.Sim_path + "\\Superposition\\fit_params_superposition"))
                return 
            
            Ehtim_Vida_plot, Brightness_ratio_str = self.plot_VIDA_templte(Ehtim_Parsers    = Ehtim_Parsers, 
                                                                           VIDA_parser      = VIDA_Parser, 
                                                                           CROP             = Center_plot, 
                                                                           crop_rel_rage    = [50, 50, 50, 50],
                                                                           Plot_Brihtness_T = False,
                                                                           Array_str        = Array_str,
                                                                           Custom_fig_title = Custom_fig_title)
            
            print(Brightness_ratio_str)

            self.Console_log_str[Array_num] += "\n" + Brightness_ratio_str

            Ehtim_Vida_plot.tight_layout()
  
            if Save_Figures:

                if not os.path.exists(self.Sim_path + "Figures\\"):
                    os.makedirs(self.Sim_path + "Figures\\")

                Ehtim_Vida_plot.savefig(self.Sim_path + 
                                        "Figures\\" + 
                                        "Ehtim_Vida_Superposition_plot_" + 
                                         Array_str +  
                                        ".png", bbox_inches = 'tight')
                
    def save_console_log_to_file(self):

        for Array_num, Array_str in enumerate(self.Arrays):

            with open(self.Sim_path + "Figures\\Flux_ratios_" + Array_str + ".csv", "w") as file:
                    print("=" * len(self.Total_flux_str.split("\n")[0]), file = file)
                    print(self.Console_log_str[Array_num], file = file)
                    print("=" * len(self.Total_flux_str.split("\n")[0]), file = file)

        print("=" * len(self.Total_flux_str.split("\n")[0]))

    def plot_contours(self, Ehtim_Parsers: list, VIDA_parser: VIDA_params_Parser, Subplot, Contour_specs: tuple[list, list]):
        
        Contour_levels, Contour_colors = Contour_specs

        Ehtim_Parser    = Ehtim_Parsers[0]
        axes_limits     = [limit * self.Units.MEGA for limit in Ehtim_Parser.WINDOW_LIMITS ]

        Intensity_ehtim_jy, _, _ = self.get_plottable_intensity_from_parsers(Ehtim_Parsers)

        # The literature (for some reason) has the X axis going positive to negative, 
        # so I invert the X axis limits
        axes_limits[0] = -axes_limits[0]
        axes_limits[1] = -axes_limits[1]

        max_value = np.max(Intensity_ehtim_jy)

        # Im not even sure what is going on with the axis limits at this points - TODO: figure out the axis inversion
        x_axis = np.linspace(axes_limits[0],axes_limits[1], Ehtim_Parser.X_PIXEL_COUNT)
        y_axis = np.linspace(axes_limits[3],axes_limits[2], Ehtim_Parser.Y_PIXEL_COUNT)

        ring_mask, dark_spot_mask = get_template_pixel_mask(VIDA_parser = VIDA_parser, 
                                                            FOV         = np.abs(axes_limits[0] - axes_limits[1]), 
                                                            N_pixels    = Ehtim_Parser.X_PIXEL_COUNT,
                                                            std_scale   = 0.5)
                    
        # Cast to a numpy array, so I can scale it by max_value
        Contour_levels = np.array(Contour_levels)

        format = {}
        for label_idx, string in zip(max_value * Contour_levels, Contour_levels):
            format[label_idx] = str(string)
        
        Contour_mask = np.ma.array(Intensity_ehtim_jy, 
                                    mask = np.logical_and(np.logical_not(dark_spot_mask), np.logical_not(ring_mask)))
                    
        Contour = Subplot.contour(x_axis, y_axis, Contour_mask, levels = max_value * Contour_levels, colors = Contour_colors)
        Labels  = Subplot.clabel(Contour, inline = True, fontsize = 12, fmt = format)
                    
        for label in Labels:
            label.set_rotation(0)
    
    def plot_VIDA_templte(self, 
                          Ehtim_Parsers: list[ehtim_Parser],
                          VIDA_parser: VIDA_params_Parser, 
                          CROP: bool, 
                          crop_rel_rage: list[float], 
                          Plot_Brihtness_T: bool,
                          Custom_fig_title: str = "",
                          Array_str: str = "") -> tuple:

        Ehtim_Parser = Ehtim_Parsers[0]

        axes_limits     = [limit * self.Units.MEGA for limit in Ehtim_Parser.WINDOW_LIMITS]
        Ehtim_image_FOV = np.abs(axes_limits[0] - axes_limits[1])  # Units of [uas]
        Ehtim_image_res = Ehtim_Parser.X_PIXEL_COUNT

        Intensity_ehtim_jy, Intensity_ehtim_T, Frequency_str = self.get_plottable_intensity_from_parsers(Ehtim_Parsers)

        template = generate_general_gaussian_template(Ehtim_image_res, VIDA_parser, Ehtim_image_FOV)

        if Plot_Brihtness_T:
            Intensity_ehtim = Intensity_ehtim_T
            colorbar_legend = r"Brightness Temperature [$10^9$K]"

        else:
            Intensity_ehtim = Intensity_ehtim_jy * self.Units.KILO
            colorbar_legend = r"Flux Per Pixel [mJy]"

        template_x_slice, slice_x_offset, template_y_slice, slice_y_offset = get_template_slices(Ehtim_image_res, template, VIDA_parser, Ehtim_image_FOV)
        Ehtim_x_slice, _, Ehtim_y_slice, _ = get_template_slices(Ehtim_image_res, Intensity_ehtim, VIDA_parser, Ehtim_image_FOV)


        # This figure is a bit large, and with the giant fontsize its going to need to be readable
        # on an A4 paper, it needs a areal big figsige parameter to look like anything
        template_fig = plt.figure(figsize = (45, 10))

        #========================= Ehtim Image Plot =========================#

        if CROP:

            x_crop_range = np.array([(-(slice_x_offset - Ehtim_image_FOV / 2) - crop_rel_rage[0]), (-(slice_x_offset - Ehtim_image_FOV / 2) + crop_rel_rage[1])])
            y_crop_range = np.array([(-(slice_y_offset - Ehtim_image_FOV / 2) - crop_rel_rage[2]), (-(slice_y_offset - Ehtim_image_FOV / 2) + crop_rel_rage[3])])
            
            # The desired crop window could "cut" outside the simulated window
            FOV_overshoot_x = max(np.absolute(x_crop_range)) - Ehtim_image_FOV / 2
            FOV_overshoot_y = max(np.absolute(y_crop_range)) - Ehtim_image_FOV / 2

            FOV_overshoot = max(FOV_overshoot_x, FOV_overshoot_y)
            
            axes_limits = [-crop_rel_rage[0], 
                            crop_rel_rage[1], 
                           -crop_rel_rage[2], 
                            crop_rel_rage[3]]

            # If it does, crop to the end of the simulated window, while keeping the aspec ratio
            if FOV_overshoot > 0:

                x_crop_range = x_crop_range - np.array([-FOV_overshoot, FOV_overshoot])
                y_crop_range = y_crop_range - np.array([-FOV_overshoot, FOV_overshoot])

                axes_limits = [-crop_rel_rage[0] + FOV_overshoot, 
                                crop_rel_rage[1] - FOV_overshoot, 
                               -crop_rel_rage[2] + FOV_overshoot, 
                                crop_rel_rage[3] - FOV_overshoot]

            x_crop_idx = (x_crop_range / (Ehtim_image_FOV / 2) + 1) / 2 * Ehtim_image_res
            y_crop_idx = (y_crop_range / (Ehtim_image_FOV / 2) + 1) / 2 * Ehtim_image_res

            x_crop_idx = x_crop_idx.astype(int)
            y_crop_idx = y_crop_idx.astype(int)

        else:
        
            x_crop_idx = [0, (Ehtim_image_res - 1)]
            y_crop_idx = [0, (Ehtim_image_res - 1)]

        crop_res_x = x_crop_idx[1] - x_crop_idx[0]
        crop_res_y = y_crop_idx[1] - y_crop_idx[0]

        # The literature (for some reason) has the X axis going positive to negative, 
        # so I invert the X axis limits

        axes_limits[0] = -axes_limits[0]
        axes_limits[1] = -axes_limits[1]
        
        # So the type checker does not complain at the imshow() call
        axes_limits = np.array(axes_limits)

        x_coords = np.linspace(axes_limits[0], axes_limits[1], crop_res_x)
        y_coords = np.linspace(axes_limits[2], axes_limits[3], crop_res_y)

        Subplot = template_fig.add_subplot(141)
        Ehtim_crop        = Intensity_ehtim[y_crop_idx[0] : y_crop_idx[1], x_crop_idx[0] : x_crop_idx[1]]
        Ehtim_crop_figure = Subplot.imshow(Ehtim_crop, cmap = "hot", extent = tuple(axes_limits))

        if CROP:

            Subplot.plot(np.zeros(crop_res_y), y_coords, "r", linewidth = 4)
            Subplot.plot(x_coords, np.zeros(crop_res_x), "b", linewidth = 4)

        else:

            Subplot.plot((slice_x_offset - Ehtim_image_FOV / 2) * np.ones(crop_res_y), y_coords, "r", linewidth = 4)
            Subplot.plot(x_coords, (slice_y_offset - Ehtim_image_FOV / 2) * np.ones(crop_res_x), "b", linewidth = 4)

        Subplot.set_xlabel(r'$\alpha_{rel}\,\,[\mu$as]', fontsize = self.Font_size)
        Subplot.set_ylabel(r'$\delta_{rel}\,\,[\mu$as]', fontsize = self.Font_size)
        plt.xticks(fontsize = self.Font_size)
        plt.yticks(fontsize = self.Font_size)

        if Custom_fig_title != "":

            if Array_str in ["2017", "2022", "2025"]:
                template_fig.suptitle("{}, viewed by EHT {}".format(Custom_fig_title, Array_str), fontsize = 1.2 * self.Font_size)
            else:
                template_fig.suptitle("{}, viewed by {}".format(Custom_fig_title, Array_str), fontsize = 1.2 * self.Font_size)

        if Array_str != "":

            if Array_str in ["2017", "2022", "2025"]:
                    Subplot.set_title("EHT {} at {}GHz".format(Array_str, Frequency_str.split(".")[0]), fontsize = self.Font_size)

            else:
                Subplot.set_title("ngEHT at {}GHz".format(Frequency_str.split(".")[0]), fontsize = self.Font_size)

        else:
            Subplot.set_title("EHT ?????? at {}GHz".format(Array_str, Frequency_str.split(".")[0]), fontsize = self.Font_size)

        colorbar = template_fig.colorbar(Ehtim_crop_figure, ax = Subplot, fraction=0.046, pad=0.04)
        colorbar.set_label(colorbar_legend, fontsize = self.Font_size, labelpad = self.Label_Pad)
        colorbar.ax.tick_params(labelsize = self.Font_size)

        #========================= Template Plot =========================#

        Subplot = template_fig.add_subplot(142)
        Subplot.imshow(template[y_crop_idx[0] : y_crop_idx[1], x_crop_idx[0] : x_crop_idx[1]], cmap = "hot", extent = tuple(axes_limits))

        if CROP:

            Subplot.plot(np.zeros(crop_res_y), y_coords, "r--", linewidth = 4)
            Subplot.plot(x_coords, np.zeros(crop_res_x), "b--", linewidth = 4)

        else:

            Subplot.plot((slice_x_offset - Ehtim_image_FOV / 2) * np.ones(crop_res_y), y_coords, "r--", linewidth = 4)
            Subplot.plot(x_coords, (slice_y_offset - Ehtim_image_FOV / 2) * np.ones(crop_res_x), "b--", linewidth = 4)

        Subplot.set_xlabel(r'$\alpha_{rel}\,\,[\mu$as]', fontsize = self.Font_size)
        Subplot.set_ylabel(r'$\delta_{rel}\,\,[\mu$as]', fontsize = self.Font_size)
        plt.xticks(fontsize = self.Font_size)
        plt.yticks(fontsize = self.Font_size)

        ring_mask, dark_spot_mask = get_template_pixel_mask(VIDA_parser, Ehtim_image_FOV, Ehtim_image_res)
        Subplot.set_title("VIDA Template", fontsize = self.Font_size)

        colorbar = template_fig.colorbar(Ehtim_crop_figure, ax = Subplot, fraction=0.046, pad=0.04)
        colorbar.set_label(colorbar_legend, fontsize = self.Font_size, labelpad = self.Label_Pad)
        colorbar.ax.tick_params(labelsize = self.Font_size)

        #------------------------ Y Slice Plot ------------------------#

        Subplot = template_fig.add_subplot(143)
        Subplot.plot(y_coords, Ehtim_y_slice[Ehtim_image_res - y_crop_idx[1] : Ehtim_image_res - y_crop_idx[0]], "r", linewidth = 4)
        Subplot.plot(y_coords, template_y_slice[Ehtim_image_res - y_crop_idx[1] : Ehtim_image_res - y_crop_idx[0]], "r--", linewidth = 4)

        Subplot.set_ylim(0, 1)
        Subplot.set_xlim(axes_limits[0], axes_limits[1])

        Subplot.set_aspect(np.absolute(axes_limits[1] - axes_limits[0]))
        Subplot.tick_params(left = True, right = False, labelleft = True,
                            labelbottom = True, bottom = True)

        Subplot.set_xlabel(r'$\delta_{rel}\,\,[\mu$as]', fontsize = self.Font_size)
        Subplot.set_ylabel('Relative Intensity', fontsize = self.Font_size)
        Subplot.set_title("Y Intensity Slice", fontsize = self.Font_size)
        plt.xticks(fontsize = self.Font_size)
        plt.yticks(fontsize = self.Font_size)

        # Subplot.imshow(ring_mask, cmap = "hot", extent = axes_limits)

        #------------------------ X Slice Plot ------------------------#

        Subplot = template_fig.add_subplot(144)
        Subplot.plot(x_coords, Ehtim_x_slice[x_crop_idx[0] : x_crop_idx[1]], "b", linewidth = 4)
        Subplot.plot(x_coords, template_x_slice[x_crop_idx[0] : x_crop_idx[1]], "b--", linewidth = 4)
        
        Subplot.set_ylim(0, 1)
        Subplot.set_xlim(axes_limits[0], axes_limits[1])

        Subplot.set_aspect(np.absolute(axes_limits[1] - axes_limits[0]))
        Subplot.tick_params(left = True, right = False, labelleft = True,
                                labelbottom = True, bottom = True)
            
        Subplot.set_xlabel(r'$\alpha_{rel}\,\,[\mu$as]', fontsize = self.Font_size)
        Subplot.set_ylabel('Relative Intensity', fontsize = self.Font_size)
        Subplot.set_title("X Intensity Slice", fontsize = self.Font_size)
        plt.xticks(fontsize = self.Font_size)
        plt.yticks(fontsize = self.Font_size)

        template_fig.subplots_adjust(bottom=0.22)

        # Subplot.imshow(dark_spot_mask, cmap = "hot", extent = axes_limits)
        
        Brightness_ratio_str = "f at {}GHz = {}".format(Frequency_str, 
                                                        get_brigness_depression_ratio(ring_mask, dark_spot_mask, Intensity_ehtim_jy))
        
        return template_fig, Brightness_ratio_str
    
    def get_plottable_intensity_from_parsers(self, Ehtim_Parsers: list):

        Ehtim_Parser    = Ehtim_Parsers[0]
        Ehtim_image_res = Ehtim_Parser.X_PIXEL_COUNT

        axes_limits     = np.array([(limit) for limit in Ehtim_Parser.WINDOW_LIMITS ]) * self.Units.MEGA
        Ehtim_image_FOV = np.abs(axes_limits[0] - axes_limits[1])  # Units of [uas]
        pixel_size      = Ehtim_image_FOV / Ehtim_image_res / self.Units.MEGA * self.Units.ARCSEC_TO_RAD  

        Intensity_ehtim_jy = np.zeros((Ehtim_image_res, Ehtim_image_res))
        Intensity_ehtim_T  = np.zeros((Ehtim_image_res, Ehtim_image_res))

        for Parser in Ehtim_Parsers:

            temp_Intensity_ehtim_jy, _ = Parser.get_plottable_ehtim_data()
            # Convert the Flux from [Jy] to brightness temperature in Giga [K]
            temp_Intensity_ehtim_T = self.Units.Spectral_density_to_T(temp_Intensity_ehtim_jy / pixel_size**2 / self.Units.W_M2_TO_JY, Parser.OBS_FREQUENCY * self.Units.GIGA) / self.Units.GIGA
            
            Intensity_ehtim_jy += temp_Intensity_ehtim_jy
            Intensity_ehtim_T  += temp_Intensity_ehtim_T

        if len(Ehtim_Parsers) > 1:
            Frequency_str = "{"

            for Freq_num, Parser in enumerate(Ehtim_Parsers):

                Frequency_str += str(int(Parser.OBS_FREQUENCY))

                if Freq_num < len(Ehtim_Parsers) - 1:

                    Frequency_str += ", "

            Frequency_str += "}"

        else:
            Frequency_str = str(Ehtim_Parser.OBS_FREQUENCY)

        return Intensity_ehtim_jy, Intensity_ehtim_T, Frequency_str

    def compare_superpos_w_single_freq(self, 
                                       Contour_specs: list,
                                       Save_Figures: bool,
                                       Custom_fig_title: str):
        
        for Array_num, Array_str in enumerate(self.Arrays):
            
            if len(self.Frequency_Bins) > 1:
                Ehtim_Parsers = [self.Ehtim_Parsers[Sim_num + Array_num * (len(self.Arrays) - 1)][self.BLUR] 
                                for Sim_num, _ in enumerate(self.Frequency_Bins)]
                        
            else:
                print("Not enough frequency bins to make a superposition!")
                return
                
            try:
                VIDA_Parser = VIDA_params_Parser(self.Sim_path + "\\Superposition\\fit_params_superposition")

            except:
                print("Could not parse VIDA template!")
                print("I looked at this path: {}".format(self.Sim_path + "\\Superposition\\fit_params_superposition"))
                return 

            Intensity_ehtim_jy, _, Frequency_str = self.get_plottable_intensity_from_parsers(Ehtim_Parsers)
            axes_limits     = np.array([(limit) for limit in Ehtim_Parsers[0].WINDOW_LIMITS ]) * self.Units.MEGA

            # The literature (for some reason) has the X axis going positive to negative, 
            # so I invert the X axis limits
            axes_limits[0] = -axes_limits[0]
            axes_limits[1] = -axes_limits[1]

            # Convert the flux to [mJy]
            Intensity_ehtim_mjy = Intensity_ehtim_jy * 1e3

            Superposition_w_contour_fig = plt.figure(figsize = (30, 9))

            Subplot = Superposition_w_contour_fig.add_subplot(131)
            Superposition_subplot = Subplot.imshow(Intensity_ehtim_mjy, cmap = "hot", extent = tuple(axes_limits), interpolation = 'bilinear')

            Subplot.set_xlabel(r'$\alpha_{rel}\,\,[\mu$as]', fontsize = self.Font_size)
            Subplot.set_ylabel(r'$\delta_{rel}\,\,[\mu$as]', fontsize = self.Font_size)
            plt.xticks(fontsize = self.Font_size)
            plt.yticks(fontsize = self.Font_size)

            Subplot.set_title("ngEHT at {}GHz".format(Frequency_str.split(".")[0]), fontsize = self.Font_size)

            if Custom_fig_title != None:
                Superposition_w_contour_fig.suptitle("{}, viewed by {}".format(Custom_fig_title, Array_str), fontsize = 1.2 * self.Font_size)

            colorbar = Superposition_w_contour_fig.colorbar(Superposition_subplot, ax = Subplot, fraction=0.046, pad=0.04)
            colorbar.set_label(r"Flux Per Pixel [mJy]", fontsize = self.Font_size, labelpad = self.Label_Pad)
            colorbar.ax.tick_params(labelsize = self.Font_size)

            self.plot_contours(Ehtim_Parsers, VIDA_Parser, Subplot, Contour_specs[0])
            
            for Sim_number, _ in enumerate(self.Frequency_Bins):
                
                Subplot = Superposition_w_contour_fig.add_subplot(132 + Sim_number)

                Index = Sim_number + Array_num * len(self.Frequency_Bins)
                Ehtim_Parser_Blur = self.Ehtim_Parsers[Index][self.BLUR]

                Intensity_ehtim_blur_jy, _ = Ehtim_Parser_Blur.get_plottable_ehtim_data()
                Intensity_ehtim_blur_mjy = Intensity_ehtim_blur_jy * self.Units.KILO

                Subplot.set_title("ngEHT at {}GHz".format(int(Ehtim_Parser_Blur.OBS_FREQUENCY)), fontsize = self.Font_size)
                Subplot.set_xlabel(r'$\alpha_{rel}\,\,[\mu$as]', fontsize = self.Font_size)
                Subplot.set_ylabel(r'$\delta_{rel}\,\,[\mu$as]', fontsize = self.Font_size)
                plt.xticks(fontsize = self.Font_size)
                plt.yticks(fontsize = self.Font_size)

                post_Convolution_T = Subplot.imshow(Intensity_ehtim_blur_mjy, interpolation = 'bilinear', cmap = 'hot', extent = tuple(axes_limits))
                colorbar = Superposition_w_contour_fig.colorbar(post_Convolution_T, ax = Subplot, fraction = 0.046, pad = 0.04)
                colorbar.set_label(r"Flux Per Pixel [mJy]", fontsize = self.Font_size, labelpad = self.Label_Pad)
                colorbar.ax.tick_params(labelsize = self.Font_size)

                self.plot_contours([self.Ehtim_Parsers[Sim_number + Array_num * (len(self.Arrays) - 1)][self.BLUR]], 
                                   VIDA_Parser, 
                                   Subplot, 
                                   Contour_specs[1 + Sim_number])
                
            Superposition_w_contour_fig.tight_layout()

            if Save_Figures:

                if not os.path.exists(self.Sim_path + "Figures\\"):
                    os.makedirs(self.Sim_path + "Figures\\")

                fig_title = "Superpos_Compare.png"

                Superposition_w_contour_fig.savefig(self.Sim_path + 
                                                    "Figures\\" + 
                                                    fig_title, bbox_inches = 'tight')
                
    def make_PIL_image(self, arr):
        
        from PIL import Image

        size = len(arr)

        # Initialize an empty image in RGB mode
        img = Image.new("RGB", (size, size))

        # Map the array to colors (0 -> white, 1 -> black)
        # Using NumPy's broadcasting for efficiency
        color_array = np.zeros((size, size, 3), dtype = np.uint8)  # Create a blank RGB array
        color_array[arr == 0] = [255, 255, 255]  # White for 0
        color_array[arr != 0] = [155, 0, 0]  # Black for non-0 values

        # Convert the NumPy array to a PIL image
        img = Image.fromarray(color_array)

        # Save the image as a PNG file
        img.save("output.png")






































# #     Ehtim_Parser_no_blur = ehtim_Parser("Results_specIDX_blur")
# #     Intensity_ehtim_no_blur, Ehtim_metadata_no_blur = Ehtim_Parser_no_blur.get_plottable_ehtim_data()

# #     Spectral_idx_figure = plt.figure()

# #     Spectral_idx_axis = Spectral_idx_figure.add_subplot(111)
# #     Spectral_idx_plot = Spectral_idx_axis.imshow(-Intensity_ehtim_no_blur, interpolation = 'bilinear', cmap = 'jet', extent = axes_limits)

# #     colorbar = Spectral_idx_figure.colorbar(Spectral_idx_plot, ax = Spectral_idx_axis, fraction=0.046, pad=0.04)
# #     colorbar.set_label(r"Spectral Index")

# #     print(np.mean(-Intensity_ehtim_no_blur))
