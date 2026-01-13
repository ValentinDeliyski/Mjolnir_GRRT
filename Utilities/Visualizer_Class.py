from numpy.typing import NDArray
from numpy import float64, bool_
from numpy import array, arctan, zeros, abs, linspace, sqrt, pi, full, ma, logical_and, logical_not, absolute, ones, swapaxes, argmax

from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib.image import AxesImage 
from matplotlib.colorbar import Colorbar 

from matplotlib import pyplot as plt
from astropy.io import fits
import os

from Support_functions.Parsers import Simulation_Parser, Units_class, ehtim_Parser, VIDA_params_Parser
from Support_functions.Image_processing import generate_general_gaussian_template, get_template_pixel_mask, get_brigness_depression_ratio, get_template_slices

class Sim_Visualizer():

    def __init__(self, Sim_path: str, 
                 Sim_Frequency_Bins: list[str], 
                 Array: list[str], 
                 Font_size: int, 
                 Label_Pad: int,
                 Common_file_name: str,
                 Respect_folder_structure: bool):

        self.Sim_Parsers: list[Simulation_Parser]    = []
        self.Ehtim_Parsers: list[list[ehtim_Parser]] = []
        self.VIDA_Parsers: list[VIDA_params_Parser]  = []
        self.Metric: str                = Common_file_name
        self.Sim_path: str              = Sim_path
        self.Arrays: list[str]          = Array
        self.Units: Units_class         = Units_class()
        self.Frequency_Bins: list[str]  = Sim_Frequency_Bins
        self.Total_flux_str: str        = ""
        self.Console_log_str: list[str] = []
        self.Respect_folder_structure   = Respect_folder_structure
        self.Font_size: int = Font_size
        self.Label_Pad: int = Label_Pad
        
        #========= Enums =========#

        self.NO_BLUR: int = 0
        self.BLUR: int    = 1

        #=========================#

        self.__make_paths(Sim_path)

        for Sim_number, _ in enumerate(self.Frequency_Bins):

            try:
                Sim_Parser = Simulation_Parser(self.Ray_tracer_paths[Sim_number])
                
                if 2 == int(Sim_Parser.Simulation_metadata["Active Simulation Mode"]):
                    print("Simulation with \"Simulation mode = 2\" aren't meant to be visualized with this script!")
                    exit()

                self.Sim_Parsers.append(Sim_Parser)
                
            except:
                print("Could not parse ray-tracer logs!")
                print("I looked at this path: {}".format(self.Ray_tracer_paths[Sim_number]))
                exit()

            if int(Sim_Parser.Simulation_metadata["Active Simulation Mode"]) != 2:

                Total_flux = Sim_Parser.get_total_flux(self.Units.M87_DISTANCE_GEOMETRICAL, unit = "mJy")
                self.Total_flux_str = self.Total_flux_str + "Total flux at {}GHz = {} [mJy]\n".format(float64(Sim_Parser.Simulation_metadata["Observation Frequency [Hz]"]) / 1e9, round(Total_flux, 4))

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

    def __make_paths(self, Sim_path: str) -> None:

        self.Ray_tracer_paths = []
        self.Ehtim_paths      = []
        
        if not self.Respect_folder_structure:
            self.Ray_tracer_paths.append(Sim_path)  
            return
            
        for freq in self.Frequency_Bins:
            self.Ray_tracer_paths.append(Sim_path + freq + "GHz\\" + "Sim_Results\\Ray_tracer_output\\" + self.Metric)

        for array in self.Arrays:
            for freq in self.Frequency_Bins:
                self.Ehtim_paths.append(Sim_path + freq + "GHz\\" + "Sim_Results\\Ehtim_" + array + "\\")
        
    def get_celestial_sphere_pattern(self, Celestial_Theta: NDArray, Celestial_Phi: NDArray) -> NDArray[float64]:
        
        def find_nearest(array, value):
            idx = (abs(array - value)).argmin()
            return array[idx]
             
        X_resolution: int = int(self.Sim_Parsers[0].Simulation_metadata["Simulation Resolution"].split(" ")[0])
        Y_resolution: int = int(self.Sim_Parsers[0].Simulation_metadata["Simulation Resolution"].split(" ")[2])

        Celestial_sphere_pattern = zeros((X_resolution, Y_resolution, 3))
        
        N_stripes: int = 40
        Stripe_width: float = pi / 400
        
        Stripe_center_phi = linspace(-pi, pi, N_stripes + 1)
        Stripe_center_theta = linspace(0 , pi, int(N_stripes / 2) + 1)
        
        for px in range(X_resolution):

            for py in range(Y_resolution):
                
                # Black rays
                if Celestial_Phi[py][px] > 1e99 or Celestial_Theta[py][px] > 1e99:
                    Celestial_sphere_pattern[px, py] = 0, 0, 0
                    
                # The red quadrant
                elif Celestial_Phi[py][px] < 0 and Celestial_Theta[py][px] < pi / 2:
                    Celestial_sphere_pattern[px, py] = 255, 0, 0
                    
                # The yellow quadrant
                elif Celestial_Phi[py][px] < 0 and Celestial_Theta[py][px] > pi / 2: 
                    Celestial_sphere_pattern[px, py] = 255, 255, 0
                    
                # The green quadrant
                elif Celestial_Phi[py][px] > 0 and Celestial_Theta[py][px] < pi / 2:
                    Celestial_sphere_pattern[px, py] = 0, 255, 0
        
                # The blue quadrant
                elif Celestial_Phi[py][px] > 0 and Celestial_Theta[py][px] > pi / 2:
                    Celestial_sphere_pattern[px, py] = 0, 0, 255
                                                     
                if (abs(Celestial_Phi[py][px] - find_nearest(Stripe_center_phi, Celestial_Phi[py][px])) < Stripe_width or
                    abs(Celestial_Theta[py][px] - find_nearest(Stripe_center_theta, Celestial_Theta[py][px])) < Stripe_width):
                            Celestial_sphere_pattern[px, py] = 0, 0, 0
                        
        return swapaxes(Celestial_sphere_pattern, 0, 1)
        
    def plot_ray_tracer_results(self, 
                                Export_data_for_Ehtim: bool, 
                                Radiation_Component: str,
                                Save_Figures: bool,
                                Custom_fig_title: str,
                                Obs_effective_distance: float,
                                Colormap: str = "seismic") -> None:

        Frequency_str_addon: str = ""

        if len(self.Frequency_Bins) == 1:
            Main_Figure: Figure = plt.figure(figsize = (20, 8))
            
        else:
            Main_Figure: Figure = plt.figure(figsize = (20, 16))

        Main_Figure.suptitle(Custom_fig_title, fontsize = self.Font_size)

        for Sim_number, Freq_str in enumerate(self.Frequency_Bins):          
            
            Obs_frequency: float = float(self.Sim_Parsers[Sim_number].Simulation_metadata["Observation Frequency [Hz]"])
            
            I_Intensity, Q_Intensity, U_Intensity, V_Intensity, Disk_redshift, Disk_flux, _, _, Celestial_theta, Celestial_phi = self.Sim_Parsers[Sim_number].get_plottable_sim_data()

            # =============== PLot the Simulated Image =============== #
            
            Fig_title: str       = "Simulated Image at {}GHz".format(int(Obs_frequency / 1e9))
            X_Slice_tile: str    = "Brightness temperature at " + r'$\delta_{\text{rel}} = 0$'
            X_Slice_y_label: str = r'$T_b\,\,[10^9\, K]$'

            # Set the X and Y axis limits, rescaling them for an observer, located at "Obs_effective_distance", rather than the simulation "Observer Distance [M]", and conver to to micro AS 
            axes_limits: NDArray[float64] = self.Sim_Parsers[Sim_number].Simulation_metadata["Observation Window Dimentions (-X,+X,-Y,+Y) [M]"].split(",")
            axes_limits = array([float64(Limit) for Limit in axes_limits])
            # axes_limits = arctan(axes_limits) * self.Units.RAD_TO_MICRO_AS

            # The literature (for some reason) has the X axis going positive to negative, 
            # so I invert the X axis limits
            axes_limits[0] = -axes_limits[0]
            axes_limits[1] = -axes_limits[1]
            
            Image_Subplot: Axes = Main_Figure.add_subplot(100 * len(self.Frequency_Bins) + 20 + (2 * Sim_number + 1))
            
            match Radiation_Component:
            
                case "Stokes I":
                    Data_to_plot: NDArray[float64] = self.Units.Spectral_density_to_T(I_Intensity / self.Units.W_M2_TO_JY, Obs_frequency) / self.Units.GIGA    
                    
                    Cmap_max: float = max(abs(Data_to_plot.flatten()))
                    Cmap_min: float = 0
                                
                    Cbar_label: str = r"Brightness Temperature [$10^9$K]"

                case "Stokes Q":
                    Data_to_plot: NDArray[float64] = Q_Intensity / max(I_Intensity.flatten()) * 100
                    
                    Cmap_max: float = max(abs(Data_to_plot.flatten())) 
                    Cmap_min: float = -Cmap_max
                    
                    Cbar_label: str = r"Q Fractional Intensity [\%]"

                case "Stokes U":
                    Data_to_plot = U_Intensity / max(I_Intensity.flatten()) * 100
                    
                    Cmap_max: float = max(abs(Data_to_plot.flatten()))
                    Cmap_min: float = -Cmap_max

                    Cbar_label: str = r"U Fractional Intensity [\%]"

                case "Stokes V":
                    Data_to_plot: NDArray[float64] = V_Intensity / max(I_Intensity.flatten()) * 100
                                       
                    Cmap_max: float = max(abs(Data_to_plot.flatten()))
                    Cmap_min: float = -Cmap_max
                    
                    Cbar_label: str = r"V Fractional Intensity [\%]"

                case "LP Fraction":
                    Data_to_plot: NDArray[float64] = sqrt(U_Intensity**2 + Q_Intensity**2) / max(abs(I_Intensity.flatten())) * 100
 
                    Cmap_max: float = max(Data_to_plot.flatten())
                    Cmap_min: float = 0
              
                    Cbar_label: str = r"LP fraction [\%]"
                    
                case "NT":
                    Data_to_plot: NDArray[float64] = Disk_flux * Disk_redshift**4 / 1e-5 
                    
                    Cmap_max: float = max((Data_to_plot.flatten()))
                    Cmap_min: float = 0
                    
                    Cbar_label: str      = r"Flux [$10^{-5}\dot{M}M^{-2}$]"
                    Fig_title: str       = r"Simulated Image"
                    X_Slice_tile: str    = r"Flux at $\delta_{\text{rel}} = 0$"
                    X_Slice_y_label: str = r"Flux $[10^{-5}\dot{M}M^{-2}]$"
                    
                    idx = argmax(Data_to_plot.flatten())
                    print(min(Data_to_plot.flatten()))
                    print(Data_to_plot.flatten()[idx])
                    print(1 / Disk_redshift.flatten()[idx] - 1)

                case "Pattern":
                    Data_to_plot = self.get_celestial_sphere_pattern(Celestial_Theta = Celestial_theta, Celestial_Phi = Celestial_phi)
                    
                    Cmap_max: float = 1.0
                    Cmap_min: float = 0.0
                                
                    Cbar_label: str = r"Brightness Temperature [$10^9$K]"

                case _:
                    print("Incorrect Radiation Component!")
                    return

            if Export_data_for_Ehtim:
                self.Sim_Parsers[Sim_number].export_ehtim_data(Spacetime = self.Sim_Parsers[Sim_number].Simulation_metadata["Spacetime [-]"], 
                                                               data = I_Intensity,
                                                               path = self.Sim_path)

            # Create the plot of the Simulated Image
            Image: AxesImage = Image_Subplot.imshow(Data_to_plot, interpolation = 'bilinear', cmap = Colormap, extent = tuple(axes_limits), vmin = Cmap_min, vmax = Cmap_max)

            colorbar: Colorbar = Main_Figure.colorbar(Image, ax = Image_Subplot, fraction = 0.046, pad = 0.04)
            colorbar.set_label(Cbar_label, fontsize = self.Font_size, labelpad = self.Label_Pad)
            colorbar.ax.tick_params(labelsize = self.Font_size)

            Image_Subplot.set_title(Fig_title, fontsize = self.Font_size)
            Image_Subplot.set_xlabel(r'$\alpha_{rel}\,\,[\mu$as]', fontsize = self.Font_size)
            Image_Subplot.set_ylabel(r'$\delta_{rel}\,\,[\mu$as]', fontsize = self.Font_size)

            plt.xticks(fontsize = self.Font_size)
            plt.yticks(fontsize = self.Font_size)
            
            if "Pattern" != Radiation_Component:

                #=============== PLot the Brigtness Temperature at y = 0 of the Simulated Image ===============#

                T_Brightness_Subplot: Axes = Main_Figure.add_subplot(100 * len(self.Frequency_Bins) + 20 + (2 * Sim_number + 2))

                # Convert the spectral density at y = 0 to brightness temperature, normalized to 10^9 Kelvin
                X_resolution: int = int(self.Sim_Parsers[Sim_number].Simulation_metadata["Simulation Resolution"].split(" ")[0])
                
                T_Brightness: NDArray[float64] = Data_to_plot[int(X_resolution / 2) - 1]
                T_Brightness_norm: float     = max(T_Brightness)
                T_Brightness_min_norm: float = min(T_Brightness)
                x_coords: NDArray[float64]   = linspace(axes_limits[0], axes_limits[1], X_resolution)

                # Create the plot of "T_b(alpha) | y = 0"
                T_Brightness_Subplot.plot(x_coords, T_Brightness)
                T_Brightness_Subplot.invert_xaxis()
                T_Brightness_Subplot.set_ylim(1.1 * T_Brightness_min_norm, 1.1 * T_Brightness_norm)
                T_Brightness_Subplot.set_title(X_Slice_tile, fontsize = self.Font_size)
                T_Brightness_Subplot.set_xlabel(r'$\alpha_{rel}\,\,[\mu$as]', fontsize = self.Font_size)
                T_Brightness_Subplot.set_ylabel(X_Slice_y_label, fontsize = self.Font_size, labelpad = self.Label_Pad)

            Frequency_str_addon += Freq_str

            plt.xticks(fontsize = self.Font_size)
            plt.yticks(fontsize = self.Font_size)

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
                axes_limits = array([(limit) for limit in Ehtim_metadata_blur]) * self.Units.MEGA

                # The literature (for some reason) has the X axis going positive to negative, 
                # so I invert the X axis limits
                axes_limits[0] = -axes_limits[0]
                axes_limits[1] = -axes_limits[1]

                pixel_size = abs(axes_limits[0] - axes_limits[1]) / Ehtim_Parser_no_Blur.X_PIXEL_COUNT / self.Units.MEGA * self.Units.ARCSEC_TO_RAD 

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

        max_value = max(Intensity_ehtim_jy.flatten())

        # Im not even sure what is going on with the axis limits at this points - TODO: figure out the axis inversion
        x_axis = linspace(axes_limits[0],axes_limits[1], Ehtim_Parser.X_PIXEL_COUNT)
        y_axis = linspace(axes_limits[3],axes_limits[2], Ehtim_Parser.Y_PIXEL_COUNT)

        ring_mask, dark_spot_mask = get_template_pixel_mask(VIDA_parser = VIDA_parser, 
                                                            FOV         = abs(axes_limits[0] - axes_limits[1]), 
                                                            N_pixels    = Ehtim_Parser.X_PIXEL_COUNT,
                                                            std_scale   = 0.5)
                    
        # Cast to a numpy array, so I can scale it by max_value
        Contour_levels = array(Contour_levels)

        format = {}
        for label_idx, string in zip(max_value * Contour_levels, Contour_levels):
            format[label_idx] = str(string)
        
        Contour_mask = ma.array(Intensity_ehtim_jy, 
                                mask = logical_and(logical_not(dark_spot_mask), logical_not(ring_mask)))
                    
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
        Ehtim_image_FOV = abs(axes_limits[0] - axes_limits[1])  # Units of [uas]
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

            x_crop_range = array([(-(slice_x_offset - Ehtim_image_FOV / 2) - crop_rel_rage[0]), (-(slice_x_offset - Ehtim_image_FOV / 2) + crop_rel_rage[1])])
            y_crop_range = array([(-(slice_y_offset - Ehtim_image_FOV / 2) - crop_rel_rage[2]), (-(slice_y_offset - Ehtim_image_FOV / 2) + crop_rel_rage[3])])
            
            # The desired crop window could "cut" outside the simulated window
            FOV_overshoot_x = max(absolute(x_crop_range)) - Ehtim_image_FOV / 2
            FOV_overshoot_y = max(absolute(y_crop_range)) - Ehtim_image_FOV / 2

            FOV_overshoot = max(FOV_overshoot_x, FOV_overshoot_y)
            
            axes_limits = [-crop_rel_rage[0], 
                            crop_rel_rage[1], 
                           -crop_rel_rage[2], 
                            crop_rel_rage[3]]

            # If it does, crop to the end of the simulated window, while keeping the aspec ratio
            if FOV_overshoot > 0:

                x_crop_range = x_crop_range - array([-FOV_overshoot, FOV_overshoot])
                y_crop_range = y_crop_range - array([-FOV_overshoot, FOV_overshoot])

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
        axes_limits = array(axes_limits)

        x_coords = linspace(axes_limits[0], axes_limits[1], crop_res_x)
        y_coords = linspace(axes_limits[2], axes_limits[3], crop_res_y)

        Subplot = template_fig.add_subplot(141)
        Ehtim_crop        = Intensity_ehtim[y_crop_idx[0] : y_crop_idx[1], x_crop_idx[0] : x_crop_idx[1]]
        Ehtim_crop_figure = Subplot.imshow(Ehtim_crop, cmap = "hot", extent = tuple(axes_limits))

        if CROP:

            Subplot.plot(zeros(crop_res_y), y_coords, "r", linewidth = 4)
            Subplot.plot(x_coords, zeros(crop_res_x), "b", linewidth = 4)

        else:

            Subplot.plot((slice_x_offset - Ehtim_image_FOV / 2) * ones(crop_res_y), y_coords, "r", linewidth = 4)
            Subplot.plot(x_coords, (slice_y_offset - Ehtim_image_FOV / 2) * ones(crop_res_x), "b", linewidth = 4)

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

            Subplot.plot(zeros(crop_res_y), y_coords, "r--", linewidth = 4)
            Subplot.plot(x_coords, zeros(crop_res_x), "b--", linewidth = 4)

        else:

            Subplot.plot((slice_x_offset - Ehtim_image_FOV / 2) * ones(crop_res_y), y_coords, "r--", linewidth = 4)
            Subplot.plot(x_coords, (slice_y_offset - Ehtim_image_FOV / 2) * ones(crop_res_x), "b--", linewidth = 4)

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

        Subplot.set_aspect(absolute(axes_limits[1] - axes_limits[0]))
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

        Subplot.set_aspect(absolute(axes_limits[1] - axes_limits[0]))
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

        axes_limits     = array([(limit) for limit in Ehtim_Parser.WINDOW_LIMITS ]) * self.Units.MEGA
        Ehtim_image_FOV = abs(axes_limits[0] - axes_limits[1])  # Units of [uas]
        pixel_size      = Ehtim_image_FOV / Ehtim_image_res / self.Units.MEGA * self.Units.ARCSEC_TO_RAD  

        Intensity_ehtim_jy = zeros((Ehtim_image_res, Ehtim_image_res))
        Intensity_ehtim_T  = zeros((Ehtim_image_res, Ehtim_image_res))

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
            axes_limits     = array([(limit) for limit in Ehtim_Parsers[0].WINDOW_LIMITS ]) * self.Units.MEGA

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