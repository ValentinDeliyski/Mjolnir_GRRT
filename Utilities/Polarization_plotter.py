import sys
import os

""" Add the parent directory of this file to the search path, 
    so this script can be ran from the "Utilities" folder """
    
parent_directory = os.path.abspath('...')
sys.path.append(parent_directory)

from Support_functions.Parsers import Simulation_Parser
from numpy import array, sqrt, flip, arctan, nan, float64, arctan2, argsort, pi, argmin, unique, argmax, arccos, append
from numpy.typing import NDArray
        
from enum import Enum

import matplotlib.pyplot as plt
from matplotlib import colormaps
from matplotlib import colors
from matplotlib.ticker import MaxNLocator
            
class Polarization_visuzlier():
    
    class Metric_enums(Enum):
        
        Numerical = 0
        Kerr = 1
    
    def __init__(self, Sim_path: str, Fig_output_path: str, fontsize: int, Numerical_title: str, Colormap: str = "plasma", Colormap_bad_color: str = "k") -> None:
        
        self.Output_path = Fig_output_path
        
        """ ================================================= Numerical metric results parsing ================================================= """
    
        self.Numerical_sim_parser = Simulation_Parser(Sim_path + "\\Numerical_results\\All_Segments_Results")
        _, _, _, _, Numerical_Disk_redshift, _, Numerical_Pol_vec_x, Numerical_Pol_vec_y, _, _ = self.Numerical_sim_parser.get_plottable_sim_data()
        
        self.Solution_mass: float = float(self.Numerical_sim_parser.Simulation_metadata["ADM Mass [M]"])
        
        X_resolution: int = int(self.Numerical_sim_parser.Simulation_metadata["Simulation Resolution"].split(" ")[0])
        Y_resolution: int = int(self.Numerical_sim_parser.Simulation_metadata["Simulation Resolution"].split(" ")[2])          
        
        Numerical_X_coords = self.Numerical_sim_parser.X_coords.reshape(X_resolution, Y_resolution)
        Numerical_X_coords = flip(Numerical_X_coords, axis = 0) / self.Solution_mass
        
        Numerical_Y_coords = self.Numerical_sim_parser.Y_coords.reshape(X_resolution, Y_resolution)
        Numerical_Y_coords = flip(Numerical_Y_coords, axis = 0) / self.Solution_mass
        
        Numerical_axes_limits = self.Numerical_sim_parser.Simulation_metadata["Observation Window Dimentions (-X,+X,-Y,+Y) [M]"].split(",")
        Numerical_axes_limits = array([float(Limit) / self.Solution_mass for Limit in Numerical_axes_limits])
        
        """ =================================================== Kerr analog results parsing =================================================== """
        
        self.Kerr_sim_parser = Simulation_Parser(Sim_path + "\\Kerr_analog\\Kerr")
        _, _, _, _, Kerr_Disk_redshift, _, Kerr_Pol_vec_x, Kerr_Pol_vec_y, _, _ = self.Kerr_sim_parser.get_plottable_sim_data()
        
        self.X_resolution: int = int(self.Kerr_sim_parser.Simulation_metadata["Simulation Resolution"].split(" ")[0])
        self.Y_resolution: int = int(self.Kerr_sim_parser.Simulation_metadata["Simulation Resolution"].split(" ")[2])          
        
        Kerr_X_coords = self.Kerr_sim_parser.X_coords.reshape(X_resolution, Y_resolution)
        Kerr_X_coords = flip(Kerr_X_coords, axis = 0) / self.Solution_mass
        
        Kerr_Y_coords = self.Kerr_sim_parser.Y_coords.reshape(X_resolution, Y_resolution)
        Kerr_Y_coords = flip(Kerr_Y_coords, axis = 0) / self.Solution_mass
        
        Kerr_axes_limits = self.Kerr_sim_parser.Simulation_metadata["Observation Window Dimentions (-X,+X,-Y,+Y) [M]"].split(",")
        Kerr_axes_limits = array([float(Limit) / self.Solution_mass for Limit in Kerr_axes_limits])
        
        """ =================================================================================================================================== """

        self.Redshift_data = [Numerical_Disk_redshift, Kerr_Disk_redshift]
        self.Polarized_Intensity_data = array([Numerical_Disk_redshift**4 * (Numerical_Pol_vec_x**2 + Numerical_Pol_vec_y**2), Kerr_Disk_redshift**4 * (Kerr_Pol_vec_x**2 + Kerr_Pol_vec_y**2)])
        
        self.Pol_vec_x_data = [Numerical_Pol_vec_x, Kerr_Pol_vec_x]
        self.Pol_vec_y_data = [Numerical_Pol_vec_y, Kerr_Pol_vec_y]
        
        self.X_coords_data = [Numerical_X_coords, Kerr_X_coords]
        self.Y_coords_data = [Numerical_Y_coords, Kerr_Y_coords]
        
        self.Axes_limits_data = [Numerical_axes_limits, Kerr_axes_limits]
        
        """ =================================================================================================================================== """
        
        self.Slice_idx = int(self.X_resolution / 2) - 1
               
        self.B_r = float(self.Numerical_sim_parser.Simulation_metadata["Disk Magnetic field geometry"].split(" ")[0][1:])
        self.B_theta = float(self.Numerical_sim_parser.Simulation_metadata["Disk Magnetic field geometry"].split(" ")[1])
        self.B_phi = float(self.Numerical_sim_parser.Simulation_metadata["Disk Magnetic field geometry"].split(" ")[2][:-1])

        self.fontsize = fontsize
        
        self.Numerical_title = Numerical_title
        
        self.Colormap = Colormap
        self.Colormap_bad_color = Colormap_bad_color

    def Plot_polarization_ticks(self, Pixel_skip_step: int, Scale_factor: float = 1) -> None:
        
        """ ========================================================== Figure setup ========================================================== """
        
        Main_Figure = plt.figure(figsize = (20, 8), layout = 'compressed')
  
        Main_Figure.suptitle(r"$i = {}$, $B_r = {}$, $B_\theta = {}$, $B_\phi = {}$".format(float(self.Numerical_sim_parser.Simulation_metadata["Observer Inclination [Deg]"]),
                                                                                            self.B_r,
                                                                                            self.B_theta,
                                                                                            self.B_phi),
                             bbox = dict(facecolor = 'none', edgecolor = 'black', boxstyle = 'round, pad = 0.2'),
                             fontsize = self.fontsize + 4,
                             y = 1)
        
        Numerical_subplot = Main_Figure.add_subplot(141)
        Numerical_subplot.set_aspect(1)
        Numerical_subplot.set_title(self.Numerical_title, fontsize = self.fontsize)
        Numerical_subplot.set_ylabel(r"$y\,[M_{\text{ADM}}]$", fontsize = self.fontsize)
        Numerical_subplot.set_xlabel(r"$x\,[M_{\text{ADM}}]$", fontsize = self.fontsize)
        
        Kerr_subplot = Main_Figure.add_subplot(142)
        Kerr_subplot.set_aspect(1)
        Kerr_subplot.set_title(r"Kerr Analog", fontsize = self.fontsize)
        Kerr_subplot.set_xlabel(r"$x\,[M_{\text{ADM}}]$", fontsize = self.fontsize)
        
        Subplot_list = [Numerical_subplot, Kerr_subplot] 

        """ ============================================================================================================================================================================================================ """
           
        Pol_intensity_plot = Main_Figure.add_subplot(143)
        
        Numerical_intensity = self.Polarized_Intensity_data[self.Metric_enums.Numerical.value][self.Slice_idx]
        Kerr_intensity = self.Polarized_Intensity_data[self.Metric_enums.Kerr.value][self.Slice_idx]
        
        Pol_intensity_plot.plot(self.X_coords_data[self.Metric_enums.Numerical.value][self.Slice_idx], Numerical_intensity, color = "m")
        Pol_intensity_plot.plot(self.X_coords_data[self.Metric_enums.Kerr.value][self.Slice_idx], Kerr_intensity, color = "C0")
        
        Intensity_max = max(max(Numerical_intensity.flatten()), max(Kerr_intensity.flatten()))
        Intensity_min = min(min(Numerical_intensity.flatten()), min(Kerr_intensity.flatten()))
        
        Pol_intensity_plot.set_aspect((max(self.X_coords_data[self.Metric_enums.Numerical.value][self.Slice_idx].flatten()) - min(self.X_coords_data[self.Metric_enums.Numerical.value][self.Slice_idx].flatten())) / (1.1 * (Intensity_max - Intensity_min)))
        Pol_intensity_plot.set_title(r"Intensity at $y = 0$ [-]", fontsize = self.fontsize)
        Pol_intensity_plot.set_xlabel(r"$x\,[M_{\text{ADM}}]$", fontsize = self.fontsize)
        Pol_intensity_plot.set_ylim(0, 1.1 * Intensity_max)
        Pol_intensity_plot.set_xlim(self.Axes_limits_data[0][0], self.Axes_limits_data[0][1])
    
        Pol_intensity_plot.minorticks_on()
        Pol_intensity_plot.tick_params(axis = "both", direction = "out")
        Pol_intensity_plot.tick_params(which = 'minor', length = 4, labelsize = self.fontsize)
        Pol_intensity_plot.tick_params(which = 'major', length = 8, labelsize = self.fontsize) 
               
        Pol_intensity_plot.legend([self.Numerical_title, "Kerr Analog"], loc = "upper left", fontsize = self.fontsize - 4)
        
        """ ============================================================================================================================================================================================================ """
        
        Colorbar_norm = colors.Normalize(vmin = float(min(self.Polarized_Intensity_data.flatten())), vmax = float(max(self.Polarized_Intensity_data.flatten())))
        Colormap = colormaps[self.Colormap]
        Colormap.set_bad(self.Colormap_bad_color)

        for _, (Intensity, Redshift, Pol_x, Pol_y, X_coord, Y_coord, Axes_limits, Subplot) in enumerate(zip(self.Polarized_Intensity_data, self.Redshift_data, self.Pol_vec_x_data, self.Pol_vec_y_data, 
                                                                                                            self.X_coords_data, self.Y_coords_data, self.Axes_limits_data, Subplot_list)):
        
            """ Set the empty pixels to nan, so the they get colored with the "bad" color """
            Intensity[Intensity == 0] = nan
            
            Disk_image = Subplot.imshow(Intensity, extent = tuple(Axes_limits), cmap = Colormap, interpolation = "nearest", norm = Colorbar_norm)
            
            Subplot.set_xlim([Axes_limits[0], Axes_limits[1]])
            Subplot.set_ylim([Axes_limits[2], Axes_limits[3]])
            
            Subplot.minorticks_on()
            Subplot.tick_params(axis = "both", direction = "out")
            Subplot.tick_params(which = 'minor', length = 4, labelsize = self.fontsize)
            Subplot.tick_params(which = 'major', length = 8, labelsize = self.fontsize)

            Pol_tick_scale = Scale_factor / sqrt(Intensity)

            X_coords_to_plot = []
            Y_coords_to_plot = []
            
            Pol_vec_x_to_plot = []
            Pol_vec_y_to_plot = []
            
            for x_idx in range(0, self.Y_resolution, Pixel_skip_step):
                
                for y_idx in range(0, self.Y_resolution, Pixel_skip_step):
                    
                    X_coords_to_plot.append(X_coord[x_idx][y_idx] - Pol_tick_scale[x_idx][y_idx] * Redshift[x_idx][y_idx]**2 * Pol_x[x_idx][y_idx] / 2)
                    Y_coords_to_plot.append(Y_coord[x_idx][y_idx] - Pol_tick_scale[x_idx][y_idx] * Redshift[x_idx][y_idx]**2 * Pol_y[x_idx][y_idx] / 2)
                    
                    Pol_vec_x_to_plot.append(Pol_tick_scale[x_idx][y_idx] * Redshift[x_idx][y_idx]**2 * Pol_x[x_idx][y_idx])
                    Pol_vec_y_to_plot.append(Pol_tick_scale[x_idx][y_idx] * Redshift[x_idx][y_idx]**2 * Pol_y[x_idx][y_idx])
                    
            Subplot.quiver(X_coords_to_plot,
                           Y_coords_to_plot,
                           Pol_vec_x_to_plot,
                           Pol_vec_y_to_plot,
                           headwidth = 0,
                           headlength = 0,
                           headaxislength = 0,
                           angles = 'xy', 
                           scale_units = 'xy',
                           scale = 1,
                           color = "black",
                           width = 0.005)
            
        colorbar = Main_Figure.colorbar(Disk_image, ax = Subplot_list, orientation = "horizontal", extend = "both", aspect = 50) # type: ignore
        colorbar.set_label(r"Intensity [-]", fontsize = self.fontsize)
        colorbar.ax.tick_params(labelsize = self.fontsize)
        
        """ ============================================================================================================================================================================================================ """

        EVPA_plot = Main_Figure.add_subplot(144)
  
        (Left_Numerical_Pol_x_slice, Right_Numerical_Pol_x_slice, 
         Left_Numerical_X_coord, Right_Numerical_X_coord) = self.Split_curve_across_brightness_depression(self.Pol_vec_x_data[self.Metric_enums.Numerical.value][self.Slice_idx], 
                                                                                                          self.X_coords_data[self.Metric_enums.Numerical.value][self.Slice_idx])
    
        
        Left_Numerical_Pol_y_slice, Right_Numerical_Pol_y_slice, _, _ = self.Split_curve_across_brightness_depression(self.Pol_vec_y_data[self.Metric_enums.Numerical.value][self.Slice_idx], 
                                                                                                                      self.X_coords_data[self.Metric_enums.Numerical.value][self.Slice_idx])
    
        
        Left_Numerical_EVPA = arctan(-Left_Numerical_Pol_x_slice / Left_Numerical_Pol_y_slice)
        Right_Numerical_EVPA = arctan(-Right_Numerical_Pol_x_slice / Right_Numerical_Pol_y_slice)
        
          
        (Left_Kerr_Pol_x_slice, Right_Kerr_Pol_x_slice, 
         Left_Kerr_X_coord, Right_Kerr_X_coord) = self.Split_curve_across_brightness_depression(self.Pol_vec_x_data[self.Metric_enums.Kerr.value][self.Slice_idx], 
                                                                                                self.X_coords_data[self.Metric_enums.Kerr.value][self.Slice_idx])
    
        
        Left_Kerr_Pol_y_slice, Right_Kerr_Pol_y_slice, _, _ = self.Split_curve_across_brightness_depression(self.Pol_vec_y_data[self.Metric_enums.Kerr.value][self.Slice_idx], 
                                                                                                            self.X_coords_data[self.Metric_enums.Kerr.value][self.Slice_idx])
    
        Left_Kerr_EVPA = arctan(-Left_Kerr_Pol_x_slice / Left_Kerr_Pol_y_slice)
        Right_Kerr_EVPA = arctan(-Right_Kerr_Pol_x_slice / Right_Kerr_Pol_y_slice)
        
        """ ====================================================================================================== """

        EVPA_plot.set_title(r"$\text{EVPA at}\,y = 0\,[\text{rad}]$", fontsize = self.fontsize)
        EVPA_plot.set_xlabel(r"$x\,[M_{\text{ADM}}]$", fontsize = self.fontsize)

        EVPA_plot.plot([max(Left_Numerical_X_coord), max(Left_Numerical_X_coord)], [-10, 10], "k")
        EVPA_plot.plot([min(Right_Numerical_X_coord), min(Right_Numerical_X_coord)], [-10, 10], "k")
        
        EVPA_plot.plot(Left_Numerical_X_coord, Left_Numerical_EVPA, color = "m")
        EVPA_plot.plot(Right_Numerical_X_coord, Right_Numerical_EVPA, color = "m")

        EVPA_plot.plot([max(Left_Kerr_X_coord), max(Left_Kerr_X_coord)], [-10, 10], "k--")
        EVPA_plot.plot([min(Right_Kerr_X_coord), min(Right_Kerr_X_coord)], [-10, 10], "k--")
        
        EVPA_plot.plot(Left_Kerr_X_coord, Left_Kerr_EVPA, "C0")
        EVPA_plot.plot(Right_Kerr_X_coord, Right_Kerr_EVPA, "C0")
        
        Numerical_EVPA_max = max(max(Left_Numerical_EVPA.flatten()), max(Right_Numerical_EVPA.flatten()))
        Numerical_EVPA_min = min(min(Left_Numerical_EVPA.flatten()), min(Right_Numerical_EVPA.flatten()))
        
        Kerr_EVPA_max = max(max(Left_Kerr_EVPA.flatten()), max(Right_Kerr_EVPA.flatten()))
        Kerr_EVPA_min = min(min(Left_Kerr_EVPA.flatten()), min(Right_Kerr_EVPA.flatten()))
        
        EVPA_max = max(Kerr_EVPA_max, Numerical_EVPA_max)
        EVPA_min = min(Kerr_EVPA_min, Numerical_EVPA_min)

        X_coord_max = max(max(Left_Numerical_X_coord.flatten()), max(Right_Numerical_X_coord.flatten()))
        X_coord_min = min(min(Left_Numerical_X_coord.flatten()), min(Right_Numerical_X_coord.flatten()))
        
        if EVPA_min < 0:
            EVPA_plot.set_ylim(1.1 * EVPA_min, 1.1 * EVPA_max)
            EVPA_plot.set_aspect((X_coord_max - X_coord_min) / (1.1 * (EVPA_max - EVPA_min)))
            
        else:
            EVPA_plot.set_ylim(0.9 * EVPA_min, 1.1 * EVPA_max)
            EVPA_plot.set_aspect((X_coord_max - X_coord_min) / (1.1 * EVPA_max - 0.9 * EVPA_min))      
            
        EVPA_plot.set_xlim(self.Axes_limits_data[0][0], self.Axes_limits_data[0][1])
                        
        EVPA_plot.minorticks_on()
        EVPA_plot.tick_params(axis = "both", direction = "out")
        EVPA_plot.tick_params(which = 'minor', length = 4, labelsize = self.fontsize)
        EVPA_plot.tick_params(which = 'major', length = 8, labelsize = self.fontsize) 
        
        if not os.path.exists(self.Output_path):
                os.makedirs(self.Output_path)

        Main_Figure.savefig(self.Output_path + "\\Whole_disk_plot_B_{}_{}_{}.pdf".format(self.B_r,
                                                                                         self.B_theta,
                                                                                         self.B_phi), 
                            bbox_inches = 'tight', 
                            pad_inches = 0.3)
        
    def Split_curve_across_brightness_depression(self, Curve_to_split: NDArray, Curve_X_coords: NDArray) -> tuple[NDArray[float64], NDArray[float64], NDArray[float64], NDArray[float64]]:

        Left_slice  = Curve_to_split[Curve_X_coords < 0]
        Right_slice = Curve_to_split[Curve_X_coords > 0]
          
        Left_X_coord  = Curve_X_coords[Curve_X_coords < 0]
        Right_X_coord = Curve_X_coords[Curve_X_coords > 0]
        
        Left_X_coord = Left_X_coord[Left_slice != 0]
        Right_X_coord = Right_X_coord[Right_slice != 0]
        
        Left_slice = Left_slice[Left_slice != 0]
        Right_slice = Right_slice[Right_slice != 0]
        
        return Left_slice, Right_slice, Left_X_coord, Right_X_coord

    def get_image_at_fixed_source_radius_numerical_only(self, r_source: float, tolerance: float, num_points: int, Scale_factor: float, axis_limits: list[float]) -> None:
         
        """ ====================================== Kerr image extraction ====================================== """
    
        X_coords = self.Numerical_sim_parser.X_coords[abs(self.Numerical_sim_parser.Source_r - r_source) < tolerance]
        Y_coords = self.Numerical_sim_parser.Y_coords[abs(self.Numerical_sim_parser.Source_r - r_source) < tolerance]

        Pol_x = self.Numerical_sim_parser.Polarization_vec_X[abs(self.Numerical_sim_parser.Source_r - r_source) < tolerance]
        Pol_y = self.Numerical_sim_parser.Polarization_vec_Y[abs(self.Numerical_sim_parser.Source_r - r_source) < tolerance]
        
        Redshift = self.Numerical_sim_parser.Disk_redshift[abs(self.Numerical_sim_parser.Source_r - r_source) < tolerance]

        Image_azimuth = arctan2(Y_coords, X_coords)

        Sorting_idx = argsort(Image_azimuth)
        X_coords = X_coords[Sorting_idx]
        Y_coords = Y_coords[Sorting_idx]
        redshift = Redshift[Sorting_idx]
        Pol_x = Pol_x[Sorting_idx]
        Pol_y = Pol_y[Sorting_idx]

        Current_length = 0
        Total_length = 0
        Final_idx_list = [0]

        for idx, _ in enumerate(Image_azimuth):
            
            if idx > 0:
                Total_length = Total_length + sqrt(X_coords[idx]**2 + Y_coords[idx]**2) * abs(Image_azimuth[idx] - Image_azimuth[idx - 1])
            
        for idx, _ in enumerate(Image_azimuth):
            
            if idx > 0:
                Current_length = Current_length + sqrt(X_coords[idx]**2 + Y_coords[idx]**2) * abs(Image_azimuth[idx] - Image_azimuth[idx - 1])
            
            if Current_length >= Total_length / num_points:
                Current_length = 0
                Final_idx_list.append(idx)
                
        Pol_x = Pol_x[Final_idx_list]
        Pol_y = Pol_y[Final_idx_list]
        X_coords = X_coords[Final_idx_list] / self.Solution_mass
        Y_coords = Y_coords[Final_idx_list] / self.Solution_mass
        redshift = redshift[Final_idx_list]
        Image_azimuth = Image_azimuth[Final_idx_list]

        Colormap = colormaps[self.Colormap] 
        Colormap.set_bad(self.Colormap_bad_color)
        
        Intensity = self.Polarized_Intensity_data[self.Metric_enums.Numerical.value]
        Intensity[Intensity == 0] = nan
        
        Axes_limits = self.Numerical_sim_parser.Simulation_metadata["Observation Window Dimentions (-X,+X,-Y,+Y) [M]"].split(",")
        Axes_limits = array([float(Limit) / self.Solution_mass for Limit in Axes_limits])
        
        Figure = plt.figure()
        Subplot = Figure.add_subplot(111)
        
        Subplot.imshow(Intensity, extent = tuple(Axes_limits), cmap = Colormap, interpolation = "nearest", vmin = 0)
    
        Final_tick_scale = Scale_factor / sqrt(Pol_x**2 + Pol_y**2)
        Subplot.quiver(X_coords - Final_tick_scale * Pol_x / 2,
                       Y_coords - Final_tick_scale * Pol_y / 2,
                       Final_tick_scale * Pol_x,
                       Final_tick_scale * Pol_y,
                       headwidth = 0,
                       headlength = 0,
                       headaxislength = 0,
                       angles = 'xy', 
                       scale_units = 'xy',
                       scale = 1,
                       color = "grey",
                       width = 0.005)
        
        Subplot.set_xlim(axis_limits[0], axis_limits[1])
        Subplot.set_ylim(axis_limits[2], axis_limits[3])       

    def __get_image_at_fixed_kerr_coords(self, r_source: float, tolerance: float, num_points: int, Curve_visualization: bool):

        """ ====================================== Kerr image extraction ====================================== """
    
        Kerr_X_coords = self.Kerr_sim_parser.X_coords[abs(self.Kerr_sim_parser.Source_r - r_source) < tolerance]
        Kerr_Y_coords = self.Kerr_sim_parser.Y_coords[abs(self.Kerr_sim_parser.Source_r - r_source) < tolerance]

        Kerr_Pol_x = self.Kerr_sim_parser.Polarization_vec_X[abs(self.Kerr_sim_parser.Source_r - r_source) < tolerance]
        Kerr_Pol_y = self.Kerr_sim_parser.Polarization_vec_Y[abs(self.Kerr_sim_parser.Source_r - r_source) < tolerance]
        
        Kerr_redshift = self.Kerr_sim_parser.Disk_redshift[abs(self.Kerr_sim_parser.Source_r - r_source) < tolerance]

        Image_azimuth = arctan2(Kerr_Y_coords, Kerr_X_coords)
        
        Image_azimuth[Image_azimuth < 0] = Image_azimuth[Image_azimuth < 0] + 2 * pi
        Sorting_idx = argsort(Image_azimuth)
        
        Image_azimuth = Image_azimuth[Sorting_idx]
        Kerr_X_coords = Kerr_X_coords[Sorting_idx]
        Kerr_Y_coords = Kerr_Y_coords[Sorting_idx]
        Kerr_redshift = Kerr_redshift[Sorting_idx]
        Kerr_Pol_x = Kerr_Pol_x[Sorting_idx]
        Kerr_Pol_y = Kerr_Pol_y[Sorting_idx]

        Total_length = 0
        Current_length = 0
        Final_idx_list = [0]

        for idx, _ in enumerate(Image_azimuth):
            
            if idx > 0:
                Total_length = Total_length + sqrt(Kerr_X_coords[idx]**2 + Kerr_Y_coords[idx]**2) * abs(Image_azimuth[idx] - Image_azimuth[idx - 1])
            
        for idx, _ in enumerate(Image_azimuth):
            
            if idx > 0:
                Current_length = Current_length + sqrt(Kerr_X_coords[idx]**2 + Kerr_Y_coords[idx]**2) * abs(Image_azimuth[idx] - Image_azimuth[idx - 1])
            
            if Current_length >= Total_length / num_points or Curve_visualization:
                Current_length = 0
                Final_idx_list.append(idx)
               
        Kerr_Pol_x = Kerr_Pol_x[Final_idx_list]
        Kerr_Pol_y = Kerr_Pol_y[Final_idx_list]
        Kerr_redshift = Kerr_redshift[Final_idx_list]
        Image_azimuth = Image_azimuth[Final_idx_list]
                
        """ ============================================ Numerical image extraction ============================================ """
    
        Kerr_X_coords = Kerr_X_coords[Final_idx_list]
        Kerr_Y_coords = Kerr_Y_coords[Final_idx_list]
    
        Nominal_X_grid = self.Numerical_sim_parser.X_coords[0:self.X_resolution]
        Nominal_Y_grid = self.Numerical_sim_parser.Y_coords[::self.Y_resolution]

        Numerical_X_coords = []
        Numerical_Y_coords = []
        Numerical_Pol_x = []
        Numerical_Pol_y = []
        Numerical_redshift = []
        
        for X_target, Y_target in zip(Kerr_X_coords, Kerr_Y_coords):
        
            X_idx = argmin(abs(Nominal_X_grid - X_target))
            Y_idx = argmin(abs(Nominal_Y_grid - Y_target))

            Numerical_X_coords.append(self.Numerical_sim_parser.X_coords[X_idx + self.Y_resolution * Y_idx])
            Numerical_Y_coords.append(self.Numerical_sim_parser.Y_coords[X_idx + self.Y_resolution * Y_idx])
            Numerical_Pol_x.append(self.Numerical_sim_parser.Polarization_vec_X[X_idx + self.Y_resolution * Y_idx])
            Numerical_Pol_y.append(self.Numerical_sim_parser.Polarization_vec_Y[X_idx + self.Y_resolution * Y_idx])
            Numerical_redshift.append(self.Numerical_sim_parser.Disk_redshift[X_idx + self.Y_resolution * Y_idx])

        Numerical_X_coords = array(Numerical_X_coords) / self.Solution_mass
        Numerical_Y_coords = array(Numerical_Y_coords) / self.Solution_mass
        Kerr_X_coords = array(Kerr_X_coords) / self.Solution_mass
        Kerr_Y_coords = array(Kerr_Y_coords) / self.Solution_mass
        
        Numerical_redshift = array(Numerical_redshift)
        Numerical_Pol_x = array(Numerical_Pol_x)
        Numerical_Pol_y = array(Numerical_Pol_y)
        
        return Image_azimuth, Kerr_Pol_x, Kerr_Pol_y, Kerr_redshift, Kerr_X_coords, Kerr_Y_coords, Numerical_Pol_x, Numerical_Pol_y, Numerical_redshift, Numerical_X_coords, Numerical_Y_coords

    def __get_image_at_fixed_numerical_coords(self, r_source: float, tolerance: float, num_points: int, Curve_visualization: bool):

        """ ====================================== Numerical image extraction ====================================== """
    
        Numerical_X_coords = self.Numerical_sim_parser.X_coords[abs(self.Numerical_sim_parser.Source_r - r_source) < tolerance]
        Numerical_Y_coords = self.Numerical_sim_parser.Y_coords[abs(self.Numerical_sim_parser.Source_r - r_source) < tolerance]

        Numerical_Pol_x = self.Numerical_sim_parser.Polarization_vec_X[abs(self.Numerical_sim_parser.Source_r - r_source) < tolerance]
        Numerical_Pol_y = self.Numerical_sim_parser.Polarization_vec_Y[abs(self.Numerical_sim_parser.Source_r - r_source) < tolerance]
        
        Numerical_redshift = self.Numerical_sim_parser.Disk_redshift[abs(self.Numerical_sim_parser.Source_r - r_source) < tolerance]

        Image_azimuth = arctan2(Numerical_Y_coords, Numerical_X_coords)
        
        Image_azimuth[Image_azimuth < 0] = Image_azimuth[Image_azimuth < 0] + 2 * pi

        Sorting_idx = argsort(Image_azimuth)
        Image_azimuth = Image_azimuth[Sorting_idx]
        Numerical_X_coords = Numerical_X_coords[Sorting_idx]
        Numerical_Y_coords = Numerical_Y_coords[Sorting_idx]
        Numerical_redshift = Numerical_redshift[Sorting_idx]
        Numerical_Pol_x = Numerical_Pol_x[Sorting_idx]
        Numerical_Pol_y = Numerical_Pol_y[Sorting_idx]

        Final_idx_list = []
        Total_length = 0
        Current_length = 0

        for idx, _ in enumerate(Image_azimuth):
            
            if idx > 0:
                Total_length = Total_length + sqrt(Numerical_X_coords[idx]**2 + Numerical_Y_coords[idx]**2) * abs(Image_azimuth[idx] - Image_azimuth[idx - 1])
        
        for idx, _ in enumerate(Image_azimuth):
            
            if idx > 0:
                Current_length = Current_length + sqrt(Numerical_X_coords[idx]**2 + Numerical_Y_coords[idx]**2) * abs(Image_azimuth[idx] - Image_azimuth[idx - 1])
            
            if Current_length > Total_length / num_points or Curve_visualization:
                Current_length = 0
                Final_idx_list.append(idx)
               
        Numerical_Pol_x = Numerical_Pol_x[Final_idx_list]
        Numerical_Pol_y = Numerical_Pol_y[Final_idx_list]
        Numerical_redshift = Numerical_redshift[Final_idx_list]
        Image_azimuth = Image_azimuth[Final_idx_list]
                
        """ ============================================ Kerr image extraction ============================================ """
    
        Numerical_X_coords = Numerical_X_coords[Final_idx_list]
        Numerical_Y_coords = Numerical_Y_coords[Final_idx_list]
    
        Nominal_X_grid = self.Kerr_sim_parser.X_coords[0:self.X_resolution]
        Nominal_Y_grid = self.Kerr_sim_parser.Y_coords[::self.Y_resolution]
        
        Kerr_X_coords = []
        Kerr_Y_coords = []
        Kerr_Pol_x = []
        Kerr_Pol_y = []
        Kerr_redshift = []
        
        for X_target, Y_target in zip(Numerical_X_coords, Numerical_Y_coords):
        
            X_idx = argmin(abs(Nominal_X_grid - X_target))
            Y_idx = argmin(abs(Nominal_Y_grid - Y_target))

            Kerr_X_coords.append(self.Kerr_sim_parser.X_coords[X_idx + self.Y_resolution * Y_idx])
            Kerr_Y_coords.append(self.Kerr_sim_parser.Y_coords[X_idx + self.Y_resolution * Y_idx])
            Kerr_Pol_x.append(self.Kerr_sim_parser.Polarization_vec_X[X_idx + self.Y_resolution * Y_idx])
            Kerr_Pol_y.append(self.Kerr_sim_parser.Polarization_vec_Y[X_idx + self.Y_resolution * Y_idx])
            Kerr_redshift.append(self.Kerr_sim_parser.Disk_redshift[X_idx + self.Y_resolution * Y_idx])

        Kerr_X_coords = array(Kerr_X_coords) / self.Solution_mass
        Kerr_Y_coords = array(Kerr_Y_coords) / self.Solution_mass
        Numerical_X_coords = array(Numerical_X_coords) / self.Solution_mass
        Numerical_Y_coords = array(Numerical_Y_coords) / self.Solution_mass
        
        Kerr_redshift = array(Kerr_redshift)
        Kerr_Pol_x = array(Kerr_Pol_x)
        Kerr_Pol_y = array(Kerr_Pol_y)
        
        return Image_azimuth, Kerr_Pol_x, Kerr_Pol_y, Kerr_redshift, Kerr_X_coords, Kerr_Y_coords, Numerical_Pol_x, Numerical_Pol_y, Numerical_redshift, Numerical_X_coords, Numerical_Y_coords

    def get_image_at_fixed_source_radius(self, r_source: float, tolerance: float, num_points: int, Scale_factor: float, axis_limits: list[float], Tick_visualization: bool, Curve_visualization: bool, Subplots: list, Source_coords: str = "Kerr", Plot_disk: bool = True) -> None:
         
        if Source_coords == "Kerr":
            Image_azimuth, Kerr_Pol_x, Kerr_Pol_y, Kerr_redshift, Kerr_X_coords, Kerr_Y_coords, Numerical_Pol_x, Numerical_Pol_y, Numerical_redshift, Numerical_X_coords, Numerical_Y_coords = self.__get_image_at_fixed_kerr_coords(r_source, tolerance, num_points, Curve_visualization)
        else:
            Image_azimuth, Kerr_Pol_x, Kerr_Pol_y, Kerr_redshift, Kerr_X_coords, Kerr_Y_coords, Numerical_Pol_x, Numerical_Pol_y, Numerical_redshift, Numerical_X_coords, Numerical_Y_coords = self.__get_image_at_fixed_numerical_coords(r_source, tolerance, num_points, Curve_visualization)

        Image_azimuth = Image_azimuth - Image_azimuth[0]
        Image_azimuth = append(Image_azimuth, Image_azimuth + 2 * pi)

        if Tick_visualization:

            """ ==================================================== Kerr Plotting ==================================================== """

            Colormap = colormaps[self.Colormap]
            Colormap.set_bad(self.Colormap_bad_color)
                
            if Subplots[1] != None:

                Kerr_Intensity = self.Polarized_Intensity_data[self.Metric_enums.Kerr.value]
                Kerr_Intensity[Kerr_Intensity == 0] = nan
                
                Kerr_axes_limits = self.Kerr_sim_parser.Simulation_metadata["Observation Window Dimentions (-X,+X,-Y,+Y) [M]"].split(",")
                Kerr_axes_limits = array([float(Limit) / self.Solution_mass for Limit in Kerr_axes_limits])
                
                if Plot_disk:
                    Subplots[1].imshow(Kerr_Intensity, extent = tuple(Kerr_axes_limits), cmap = Colormap, interpolation = "nearest", vmin = 0)
        
                Final_tick_scale = Scale_factor / sqrt(Kerr_Pol_x**2 + Kerr_Pol_y**2 + 1e-10)
                Subplots[1].quiver(Kerr_X_coords - Final_tick_scale * Kerr_Pol_x / 2,
                                    Kerr_Y_coords - Final_tick_scale * Kerr_Pol_y / 2,
                                    Final_tick_scale * Kerr_Pol_x,
                                    Final_tick_scale * Kerr_Pol_y,
                                    headwidth = 0,
                                    headlength = 0,
                                    headaxislength = 0,
                                    angles = 'xy', 
                                    scale_units = 'xy',
                                    scale = 1,
                                    color = "grey",
                                    width = 0.005)
                
                Subplots[1].set_xlim(axis_limits[0], axis_limits[1])
                Subplots[1].set_ylim(axis_limits[2], axis_limits[3])
                
            """ ==================================================== Numerical Plotting ==================================================== """
   
            if Subplots[0] != None:
    
                Numerical_Intensity = self.Polarized_Intensity_data[self.Metric_enums.Numerical.value]
                Numerical_Intensity[Numerical_Intensity == 0] = nan
                
                Numerical_axes_limits = self.Numerical_sim_parser.Simulation_metadata["Observation Window Dimentions (-X,+X,-Y,+Y) [M]"].split(",")
                Numerical_axes_limits = array([float(Limit) / self.Solution_mass for Limit in Numerical_axes_limits])
                
                if Plot_disk:
                    Subplots[0].imshow(Numerical_Intensity, extent = tuple(Numerical_axes_limits), cmap = Colormap, interpolation = "nearest", vmin = 0)
                
                Final_tick_scale = Scale_factor / sqrt(Numerical_Pol_x**2 + Numerical_Pol_y**2 + 1e-10)
                Subplots[0].quiver(Numerical_X_coords - Final_tick_scale * Numerical_Pol_x / 2,
                                    Numerical_Y_coords - Final_tick_scale * Numerical_Pol_y / 2,
                                    Final_tick_scale * Numerical_Pol_x,
                                    Final_tick_scale * Numerical_Pol_y,
                                    headwidth = 0,
                                    headlength = 0,
                                    headaxislength = 0,
                                    angles = 'xy', 
                                    scale_units = 'xy',
                                    scale = 1,
                                    color = "grey",
                                    width = 0.005)
                
                Subplots[0].set_xlim(axis_limits[0], axis_limits[1])
                Subplots[0].set_ylim(axis_limits[2], axis_limits[3])
        
        if Curve_visualization:
            
            if Subplots[1] != None:
    
                Kerr_X_coords = append(Kerr_X_coords, Kerr_X_coords[0])
                Kerr_Y_coords = append(Kerr_Y_coords, Kerr_Y_coords[0])
    
                Subplots[1].plot(Kerr_X_coords, Kerr_Y_coords, "r--")
                Subplots[1].set_ylabel(r"y [$M_\text{ADM}$]", fontsize = self.fontsize)
                Subplots[1].set_xlabel(r"x [$M_\text{ADM}$]", fontsize = self.fontsize)
                Subplots[1].set_title(r'Kerr Analog', usetex=True, fontsize = self.fontsize)
                
                Subplots[1].minorticks_on()
                Subplots[1].tick_params(axis = "both", direction = "out")
                Subplots[1].tick_params(which = 'minor', length = 4, labelsize = self.fontsize)
                Subplots[1].tick_params(which = 'major', length = 8, labelsize = self.fontsize) 
                
                Subplots[1].xaxis.set_major_locator(MaxNLocator(nbins = 5))
                Subplots[1].yaxis.set_major_locator(MaxNLocator(nbins = 5))

            if Subplots[0] != None:

                Numerical_X_coords = append(Numerical_X_coords, Numerical_X_coords[0])
                Numerical_Y_coords = append(Numerical_Y_coords, Numerical_Y_coords[0])

                Subplots[0].plot(Numerical_X_coords, Numerical_Y_coords, "r--")
                Subplots[0].set_ylabel(r"y [$M_\text{ADM}$]", fontsize = self.fontsize)
                Subplots[0].set_xlabel(r"x [$M_\text{ADM}$]", fontsize = self.fontsize)
                Subplots[0].set_title(self.Numerical_title, fontsize = self.fontsize)
                
                Subplots[0].minorticks_on()
                Subplots[0].tick_params(axis = "both", direction = "out")
                Subplots[0].tick_params(which = 'minor', length = 4, labelsize = self.fontsize)
                Subplots[0].tick_params(which = 'major', length = 8, labelsize = self.fontsize)
                # Subplots[0].set_xticks(xy_ticks)
                # Subplots[0].set_xticklabels(xy_tick_labels)
                # Subplots[0].set_yticks(xy_ticks)
                # Subplots[0].set_yticklabels(xy_tick_labels)
                
                Subplots[0].xaxis.set_major_locator(MaxNLocator(nbins = 5))
                Subplots[0].yaxis.set_major_locator(MaxNLocator(nbins = 5))
            
            """" ================================================================================================================================================== """
            
            x_ticks = [0, pi / 2, pi, 3 * pi / 2, 2 * pi]
            x_tick_labels = ["0", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2\pi$"]
            
            if Subplots[2] != None:
                
                Numerical_Intensity = Numerical_redshift**4 * (Numerical_Pol_x**2 + Numerical_Pol_y**2)
                Numerical_Intensity = append(Numerical_Intensity, Numerical_Intensity)
                
                Kerr_Intensity = Kerr_redshift**4 * (Kerr_Pol_x**2 + Kerr_Pol_y**2)
                Kerr_Intensity = append(Kerr_Intensity, Kerr_Intensity)
                
                Subplots[2].plot(Image_azimuth, Numerical_Intensity, "m")
                Subplots[2].plot(Image_azimuth, Kerr_Intensity, "C0")
                
                Intensity_max = max(max(Kerr_Intensity), max(Numerical_Intensity))
                
                Fig_padding = 0.1 * Intensity_max
                
                Subplots[2].set_ylim(-Fig_padding, Intensity_max + 3 * Fig_padding)
                Subplots[2].set_aspect(0.91 * 2 * pi  / ((Intensity_max + 4 * Fig_padding)))
                Subplots[2].yaxis.get_offset_text().set_fontsize(self.fontsize - 8)
                
                Subplots[2].minorticks_on()
                Subplots[2].tick_params(axis = "both", direction = "out")
                Subplots[2].tick_params(which = 'minor', length = 4, labelsize = self.fontsize)
                Subplots[2].tick_params(which = 'major', length = 8, labelsize = self.fontsize) 

                Subplots[2].set_xlim(0, 2 * pi)
                Subplots[2].set_title(r"Intensity [-]", fontsize = self.fontsize)
                Subplots[2].set_xlabel(r"Image azimuth [rad]", fontsize = self.fontsize)
                Subplots[2].legend([self.Numerical_title, "Kerr Analog"], fontsize = self.fontsize - 10, loc = 'upper center', ncols = 2)
                
                Subplots[2].set_xticks(x_ticks)
                Subplots[2].set_xticklabels(x_tick_labels)
                
                Subplots[2].ticklabel_format(style = 'sci', axis = 'y', scilimits=(0, 0))
                Subplots[2].yaxis.set_major_locator(MaxNLocator(nbins = 6))
                
            """" ================================================================================================================================================== """
            
            if Subplots[2] != None:
                    
                EVPA_Kerr = arctan(-Kerr_Pol_x / Kerr_Pol_y)
                EVPA_Kerr = append(EVPA_Kerr, EVPA_Kerr)
                
                EVPA_Numerical = arctan(-Numerical_Pol_x / Numerical_Pol_y)
                EVPA_Numerical = append(EVPA_Numerical, EVPA_Numerical)
                
                EVPA_Kerr_idx_split = self.split_EVPA(EVPA = EVPA_Kerr)
                EVPA_Numerical_idx_split = self.split_EVPA(EVPA = EVPA_Numerical)
                
                for idx in range(len(EVPA_Kerr_idx_split) - 1):
                    Subplots[3].plot(Image_azimuth[EVPA_Kerr_idx_split[idx] : EVPA_Kerr_idx_split[idx + 1]] - Image_azimuth[0], EVPA_Kerr[EVPA_Kerr_idx_split[idx] : EVPA_Kerr_idx_split[idx + 1]], "C0")       
                    Subplots[3].plot(Image_azimuth[EVPA_Kerr_idx_split[idx + 1] - 1 : EVPA_Kerr_idx_split[idx + 1] + 1] - Image_azimuth[0], EVPA_Kerr[EVPA_Kerr_idx_split[idx + 1] - 1 : EVPA_Kerr_idx_split[idx + 1] + 1], "C0--")               
                    
                for idx in range(len(EVPA_Numerical_idx_split) - 1):        
                    Subplots[3].plot(Image_azimuth[EVPA_Numerical_idx_split[idx] : EVPA_Numerical_idx_split[idx + 1]] - Image_azimuth[0], EVPA_Numerical[EVPA_Numerical_idx_split[idx] : EVPA_Numerical_idx_split[idx + 1]], "m")
                    Subplots[3].plot(Image_azimuth[EVPA_Numerical_idx_split[idx + 1] - 1 : EVPA_Numerical_idx_split[idx + 1] + 1] - Image_azimuth[0], EVPA_Numerical[EVPA_Numerical_idx_split[idx + 1] - 1 : EVPA_Numerical_idx_split[idx + 1] + 1], "m--")
                    
                EVPA_max = max(max(EVPA_Kerr), max(EVPA_Numerical))
                EVPA_min = min(min(EVPA_Kerr), min(EVPA_Numerical))
                
                Fig_padding = 0.1 * abs(EVPA_max - EVPA_min)
                
                Subplots[3].set_ylim(EVPA_min - Fig_padding, EVPA_max + Fig_padding)
                Subplots[3].set_aspect(0.91 * 2 * pi  / (EVPA_max - EVPA_min + 2 * Fig_padding))

                Subplots[3].minorticks_on()
                Subplots[3].tick_params(axis = "both", direction = "out")
                Subplots[3].tick_params(which = 'minor', length = 4, labelsize = self.fontsize)
                Subplots[3].tick_params(which = 'major', length = 8, labelsize = self.fontsize) 

                Subplots[3].set_xlim(0, 2 * pi)
                Subplots[3].set_title(r"EVPA [rad]", fontsize = self.fontsize)
                Subplots[3].set_xlabel(r"Image azimuth [rad]", fontsize = self.fontsize)
                Subplots[3].set_xticks(x_ticks)
                Subplots[3].set_xticklabels(x_tick_labels)
                
                y_ticks = [-pi / 2, -pi / 4, 0, pi / 4, pi / 2]
                y_tick_labels = [r"$-\frac{\pi}{2}$", r"$-\frac{\pi}{4}$", "0", r"$\frac{\pi}{4}$", r"$\frac{\pi}{2}$"]
                
                Subplots[3].set_yticks(y_ticks)
                Subplots[3].set_yticklabels(y_tick_labels)
            
            """" ================================================================================================================================================== """

            if Subplots[4] != None:   
                                        
                Numerical_Intensity = Numerical_redshift**4 * (Numerical_Pol_x**2 + Numerical_Pol_y**2)
                Numerical_Intensity = append(Numerical_Intensity, Numerical_Intensity)
                
                Kerr_Intensity = Kerr_redshift**4 * (Kerr_Pol_x**2 + Kerr_Pol_y**2)
                Kerr_Intensity = append(Kerr_Intensity, Kerr_Intensity)
                
                Subplots[4].plot(Image_azimuth, Numerical_Intensity - Kerr_Intensity)
                Delta_I_max = max(Numerical_Intensity - Kerr_Intensity)
                Delta_I_min = min(Numerical_Intensity - Kerr_Intensity)
                
                Fig_padding = 0.1 * abs(Delta_I_max - Delta_I_min)
                
                Subplots[4].set_ylim(Delta_I_min - Fig_padding, Delta_I_max + Fig_padding)
                Subplots[4].set_aspect(0.91 * 2 * pi  / (Delta_I_max - Delta_I_min + 2 * Fig_padding))
                Subplots[4].yaxis.get_offset_text().set_fontsize(self.fontsize - 8)
                
                Subplots[4].minorticks_on()
                Subplots[4].tick_params(axis = "both", direction = "out")
                Subplots[4].tick_params(which = 'minor', length = 4, labelsize = self.fontsize)
                Subplots[4].tick_params(which = 'major', length = 8, labelsize = self.fontsize) 
                
                Subplots[4].set_xlim(0, 2 * pi)
                Subplots[4].set_title(r"$\Delta I$ [-]", fontsize = self.fontsize)
                Subplots[4].set_xlabel(r"Image azimuth [rad]", fontsize = self.fontsize)
                Subplots[4].set_xticks(x_ticks)
                Subplots[4].set_xticklabels(x_tick_labels)
                
                Subplots[4].ticklabel_format(style = 'sci', axis = 'y', scilimits=(0, 0))
                Subplots[4].yaxis.set_major_locator(MaxNLocator(nbins = 6))
            
            """" ================================================================================================================================================== """
            
            if Subplots[5] != None:
                
                Delta_EVPA = arctan((Numerical_Pol_x * Kerr_Pol_y - Numerical_Pol_y * Kerr_Pol_x) / (Numerical_Pol_x * Kerr_Pol_x + Numerical_Pol_y * Kerr_Pol_y))
                Delta_EVPA = append(Delta_EVPA, Delta_EVPA)
                idx = argmax(abs(Delta_EVPA))
                
                print("X = {}".format(Kerr_X_coords[idx]))
                print("Y = {}".format(Kerr_Y_coords[idx]))
                
                Subplots[0].plot(Kerr_X_coords[idx], Kerr_Y_coords[idx], "ro")
                
                Delta_EVPA_idx_split = self.split_EVPA(EVPA = Delta_EVPA)
                
                for idx in range(len(Delta_EVPA_idx_split) - 1):
                    Subplots[5].plot(Image_azimuth[Delta_EVPA_idx_split[idx] : Delta_EVPA_idx_split[idx + 1]] - Image_azimuth[0], Delta_EVPA[Delta_EVPA_idx_split[idx] : Delta_EVPA_idx_split[idx + 1]], "C0")       
                    Subplots[5].plot(Image_azimuth[Delta_EVPA_idx_split[idx + 1] - 1 : Delta_EVPA_idx_split[idx + 1] + 1] - Image_azimuth[0], Delta_EVPA[Delta_EVPA_idx_split[idx + 1] - 1 : Delta_EVPA_idx_split[idx + 1] + 1], "C0--")               
                           
                Delta_EVPA_max = max(Delta_EVPA)
                Delta_EVPA_min = min(Delta_EVPA)
                
                Fig_padding = 0.1 * abs(Delta_EVPA_max - Delta_EVPA_min)
                
                Subplots[5].set_ylim(Delta_EVPA_min - Fig_padding, Delta_EVPA_max + Fig_padding)
                Subplots[5].set_aspect(0.91 * 2 * pi  / (Delta_EVPA_max - Delta_EVPA_min + 2 * Fig_padding))
                
                Subplots[5].minorticks_on()
                Subplots[5].tick_params(axis = "both", direction = "out")
                Subplots[5].tick_params(which = 'minor', length = 4, labelsize = self.fontsize)
                Subplots[5].tick_params(which = 'major', length = 8, labelsize = self.fontsize) 
            
                Subplots[5].set_xlim(0, 2 * pi)
                Subplots[5].set_title(r"$\Delta$EVPA [rad]", fontsize = self.fontsize)
                Subplots[5].set_xlabel(r"Image azimuth [rad]", fontsize = self.fontsize)
                Subplots[5].set_xticks(x_ticks)
                Subplots[5].set_xticklabels(x_tick_labels)
                
                Subplots[5].yaxis.set_major_locator(MaxNLocator(nbins = 6))
                
                idx = argmax(abs(Delta_EVPA))
                print("Azimuth at max Delta EVPA = {} pi".format(round(Image_azimuth[idx] / pi, 2)))
                print("Delta EVPA max = {} pi".format(round(Delta_EVPA[idx] / pi, 4)))
        
    def split_EVPA(self, EVPA):

        DISCONTINUITY_TRESHOLD = pi / 2

        branch_index = []
        branch_index.append(0)

        for index, _ in enumerate(EVPA):

            if index > 1 and abs(EVPA[index] - EVPA[index - 1]) > DISCONTINUITY_TRESHOLD:

                branch_index.append(index)

        branch_index.append(len(EVPA) - 1)

        return branch_index

        
if __name__ == "__main__":
    
    plt.rcParams['axes.titlepad'] = 15
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = 'serif' # or 'sans-serif' or 'monospace'
    plt.rcParams['font.serif'] = 'cmr10'
    plt.rcParams['font.sans-serif'] = 'cmss10'
    plt.rcParams['font.monospace'] = 'cmtt10'
    plt.rcParams["axes.formatter.use_mathtext"] = True
    plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}' + '\n' + r'\usepackage{xcolor}'
    
    horizon_radii = [0.01, 0.2] 
    
    # for idx, Model in enumerate(["II", "V"]):
        
    #     for B_field in [[0.0, 1.0, 0.0]]:
        
    #         for inc in [17]:
    
    #             Sim_path = "E:\\Numerical_metric_runs\\Zero_curvature\\Config_{}_{}_deg_B_{}_{}_{}_min_order_0_max_order_0".format(Model, inc, B_field[0], B_field[1], B_field[2])
                
    #             Visualizer_instance = Polarization_visuzlier(Sim_path = Sim_path, 
    #                                                         fontsize = 22, 
    #                                                         Fig_output_path = "E:\\Numerical_metric_runs\\Zero_curvature\\Figures\\Config_{}_{}_deg_min_order_0_max_order_0".format(Model, inc),
    #                                                         Numerical_title = r"Model {}$^0_{{{}}}$".format(Model, horizon_radii[idx]))
                
    #             Visualizer_instance.Plot_polarization_ticks(Pixel_skip_step = 32, Scale_factor = 0.75)
                
    Sim_path = "C:/Users/Valur/Documents/Repos/Mjolnir_GRRT/Utilities/Numerical_metric_runs/Zero_curvature/Config_V_70_deg_B_0.0_1.0_0.0_min_order_0_max_order_0_zoom"
                
    Visualizer_instance = Polarization_visuzlier(Sim_path = Sim_path, 
                                                            fontsize = 22, 
                                                            Fig_output_path = "E:\\Numerical_metric_runs\\Zero_curvature\\Figures\\Config_II_70_deg_min_order_0_max_order_0",
                                                            Numerical_title = r"Model V$^0_{0.2}$")
                    
    Figure = plt.figure(figsize = (13, 10.5), layout = "compressed")
    
    Figure.suptitle(r"$i = {}$, $B_r = {}$, $B_\theta = {}$, $B_\phi = {}$".format(float(Visualizer_instance.Numerical_sim_parser.Simulation_metadata["Observer Inclination [Deg]"]),
                                                                                         Visualizer_instance.B_r,
                                                                                         Visualizer_instance.B_theta,
                                                                                         Visualizer_instance.B_phi),
                    fontsize = Visualizer_instance.fontsize + 4,
                    bbox = dict(facecolor = 'none', edgecolor = 'black', boxstyle = 'round, pad = 0.3'),
                    y = 0.99)
    
    Numerical_subplot = Figure.add_subplot(231)
    Kerr_subplot = Figure.add_subplot(234)
    Intensity_subplot = Figure.add_subplot(232)
    Delta_Intensity_subplot = Figure.add_subplot(235)
    EVPA_subplot = Figure.add_subplot(233)
    Delta_EVPA_subplot = Figure.add_subplot(236)
    
    Subplots = [Numerical_subplot, Kerr_subplot, Intensity_subplot, EVPA_subplot, Delta_Intensity_subplot, Delta_EVPA_subplot]
    
    for orbit_mult in [1]:
        
        Visualizer_instance.get_image_at_fixed_source_radius(r_source = orbit_mult * 1.2121, 
                                                            tolerance = 0.02, 
                                                            num_points = 65, 
                                                            Scale_factor = 1.25,  
                                                            axis_limits = [-3, 4.5, -2, 4.5],
                                                            Tick_visualization = True,
                                                            Curve_visualization = False,
                                                            Subplots = Subplots)
        
        Visualizer_instance.get_image_at_fixed_source_radius(r_source = orbit_mult * 1.2121, 
                                                            tolerance = 0.0003, 
                                                            num_points = 65, 
                                                            Scale_factor = 1.25, 
                                                            axis_limits = [-4, 4, -2, 4],
                                                            Tick_visualization = False,
                                                            Curve_visualization = True,
                                                            Subplots = Subplots)
    
    # # Visualizer_instance.get_image_at_fixed_source_radius(r_source = 5 * 1.2121, 
    # #                                                      tolerance = 0.02, 
    # #                                                      length_step = 0.35, 
    # #                                                      Scale_factor = 1.75, 
    # #                                                      axis_limits = [-15, 15, -5, 11],
    # #                                                      Tick_visualization = True,
    # #                                                      Curve_visualization = False,
    # #                                                      Subplots = Subplots,
    # #                                                      Source_coords = "Kerr",
    # #                                                      Plot_disk = False)
    
    # # Visualizer_instance.get_image_at_fixed_source_radius(r_source = 5 * 1.2121, 
    # #                                                      tolerance = 0.005,
    # #                                                      length_step = 0.006, 
    # #                                                      Scale_factor = 1.75, 
    # #                                                      axis_limits = [-15, 15, -5, 11],
    # #                                                      Tick_visualization = False,
    # #                                                      Curve_visualization = True,
    # #                                                      Subplots = Subplots,
    # #                                                      Source_coords = "Kerr",
    # #                                                      Plot_disk = False)
 
    # # Visualizer_instance.get_image_at_fixed_source_radius_numerical_only(r_source = 0.24029573045223954249, 
    # #                                                                     tolerance = 0.002, 
    # #                                                                     length_step = 0.21, 
    # #                                                                     Scale_factor = 0.75, 
    # #                                                                     axis_limits = [-4, 4, -4, 4])
    
    # # if not os.path.exists(Visualizer_instance.Output_path):
    # #     os.makedirs(Visualizer_instance.Output_path)
            
    # Figure.savefig(Visualizer_instance.Output_path + "\\ISCO_zoom_B_{}_{}_{}_v3.pdf".format(Visualizer_instance.B_r,
    #                                                                                          Visualizer_instance.B_theta,
    #                                                                                          Visualizer_instance.B_phi), bbox_inches = 'tight',
    #                                                                                          pad_inches = 0.15)
    
    plt.show()