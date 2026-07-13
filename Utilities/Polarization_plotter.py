import sys
import os

""" Add the parent directory of this file to the search path, 
    so this script can be ran from the "Utilities" folder """
    
parent_directory = os.path.abspath('...')
sys.path.append(parent_directory)

from dataclasses import dataclass

from Support_functions.Parsers import Simulation_Parser
from numpy import array, sqrt, flip, arctan, nan, float64, arctan2, argsort, pi, argmin, unique, argmax, arccos, append, logical_and, logical_or, linspace
from numpy.typing import NDArray
        
from enum import Enum

import matplotlib.pyplot as plt
from matplotlib import colormaps
from matplotlib import colors
from matplotlib.ticker import MaxNLocator

@dataclass
class Sim_data:
    Sim_name: str
    
    Redshift: NDArray[float64]
    Pol_vec_x: NDArray[float64]
    Pol_vec_y: NDArray[float64]
    
    X_res: int
    Y_res: int
    
    X_coords: NDArray[float64]
    Y_coords: NDArray[float64]
    
    Axes_limits: NDArray[float64]
    
    B_field: NDArray[float64]
    
    ADM_mass: float
    
    Pixel_skip: int
            
class Polarization_visuzlier():
    
    class Panel_Enums(Enum):
        
        Lower = 1
        Upper = 0
    
    def __init__(self, Panel_A_path: str, Panel_B_path: str, Left_Panel_sim_name: str, Right_panel_sim_name: str, Fig_output_path: str, fontsize: int, Colormap: str = "plasma", Colormap_bad_color: str = "k") -> None:
        
        self.Output_path = Fig_output_path
        
        self.Left_Panel_sim_name = Left_Panel_sim_name
        self.Right_panel_sim_name = Right_panel_sim_name
        
        """ ================================================= Left/Upper panel results parsing ================================================= """
    
        self.Panel_A_parser = Simulation_Parser(Panel_A_path)
        _, _, _, _, Redshift, _, Pol_vec_x, Pol_vec_y, _, _, _, _, _, _ = self.Panel_A_parser.get_plottable_sim_data()
        
        X_res: int = int(self.Panel_A_parser.Simulation_metadata["Simulation Resolution"].split(" ")[0])
        Y_res: int = int(self.Panel_A_parser.Simulation_metadata["Simulation Resolution"].split(" ")[2])          
        
        ADM_mass = float(self.Panel_A_parser.Simulation_metadata["ADM Mass [M]"])
                
        X_coords = self.Panel_A_parser.X_coords.reshape(X_res, Y_res)
        X_coords = flip(X_coords, axis = 0) / ADM_mass
        
        Y_coords = self.Panel_A_parser.Y_coords.reshape(X_res, Y_res)
        Y_coords = flip(Y_coords, axis = 0) / ADM_mass
                
        Axes_limits = self.Panel_A_parser.Simulation_metadata["Observation Window Dimentions (-X,+X,-Y,+Y) [M]"].split(",")
        Axes_limits = array([float(Limit) / ADM_mass for Limit in Axes_limits])
        
        B_field = array([float(self.Panel_A_parser.Simulation_metadata["Disk Magnetic field geometry"].split(" ")[0][1:]), 
                         float(self.Panel_A_parser.Simulation_metadata["Disk Magnetic field geometry"].split(" ")[1]), 
                         float(self.Panel_A_parser.Simulation_metadata["Disk Magnetic field geometry"].split(" ")[2][:-1])])
        
        self.Panel_A_context = Sim_data(Sim_name = Left_Panel_sim_name,
                                        Redshift = Redshift,
                                        Pol_vec_x = Pol_vec_x,
                                        Pol_vec_y = Pol_vec_y,
                                        X_res = X_res,
                                        Y_res = Y_res,
                                        Pixel_skip = int(Y_res / 32),
                                        X_coords = X_coords,
                                        Y_coords = Y_coords,
                                        Axes_limits = Axes_limits,
                                        B_field = B_field,
                                        ADM_mass = ADM_mass)
        
        """ =================================================== Right/Lower analog results parsing =================================================== """
        
        self.Panel_B_parser = Simulation_Parser(Panel_B_path)
        _, _, _, _, Redshift, _, Pol_vec_x, Pol_vec_y, _, _, _, _, _, _ = self.Panel_B_parser.get_plottable_sim_data()
        
        X_res: int = int(self.Panel_B_parser.Simulation_metadata["Simulation Resolution"].split(" ")[0])
        Y_res: int = int(self.Panel_B_parser.Simulation_metadata["Simulation Resolution"].split(" ")[2])          
        
        ADM_mass = float(self.Panel_B_parser.Simulation_metadata["ADM Mass [M]"])
                
        X_coords = self.Panel_B_parser.X_coords.reshape(X_res, Y_res)
        X_coords = flip(X_coords, axis = 0) / ADM_mass
        
        Y_coords = self.Panel_B_parser.Y_coords.reshape(X_res, Y_res)
        Y_coords = flip(Y_coords, axis = 0) / ADM_mass
                
        Axes_limits = self.Panel_B_parser.Simulation_metadata["Observation Window Dimentions (-X,+X,-Y,+Y) [M]"].split(",")
        Axes_limits = array([float(Limit) / ADM_mass for Limit in Axes_limits])
        
        B_field = array([float(self.Panel_B_parser.Simulation_metadata["Disk Magnetic field geometry"].split(" ")[0][1:]), 
                         float(self.Panel_B_parser.Simulation_metadata["Disk Magnetic field geometry"].split(" ")[1]), 
                         float(self.Panel_B_parser.Simulation_metadata["Disk Magnetic field geometry"].split(" ")[2][:-1])])
        
        self.Panel_B_context = Sim_data(Sim_name = Left_Panel_sim_name,
                                        Redshift = Redshift,
                                        Pol_vec_x = Pol_vec_x,
                                        Pol_vec_y = Pol_vec_y,
                                        X_res = X_res,
                                        Y_res = Y_res,
                                        Pixel_skip = int(Y_res / 64),
                                        X_coords = X_coords,
                                        Y_coords = Y_coords,
                                        Axes_limits = Axes_limits,
                                        B_field = B_field,
                                        ADM_mass = ADM_mass)

        """ =================================================================================================================================== """
        
        self.fontsize = fontsize
        
        self.Colormap = Colormap
        self.Colormap_bad_color = Colormap_bad_color

    def Plot_polarization_ticks(self, Scale_factor: float = 1, Pattern_only: bool = False) -> None:
        
        """ ========================================================== Figure setup ========================================================== """
        
        if Pattern_only:
            Main_Figure = plt.figure(figsize = (10, 8), layout = 'compressed')
        else:
            Main_Figure = plt.figure(figsize = (20, 8), layout = 'compressed')
  
        Main_Figure.suptitle(r"$i = {}$, $B_r = {}$, $B_\theta = {}$, $B_\phi = {}$".format(float(self.Panel_A_parser.Simulation_metadata["Observer Inclination [Deg]"]),
                                                                                            self.Panel_A_context.B_field[0],
                                                                                            self.Panel_A_context.B_field[1],
                                                                                            self.Panel_A_context.B_field[2]),
                             bbox = dict(facecolor = 'none', edgecolor = 'black', boxstyle = 'round, pad = 0.2'),
                             fontsize = self.fontsize + 4,
                             y = 1)
        
        if Pattern_only:
            Subplot_A = Main_Figure.add_subplot(121)
            Subplot_B = Main_Figure.add_subplot(122)
        else:
            Subplot_A = Main_Figure.add_subplot(141)
            Subplot_B = Main_Figure.add_subplot(142)
            
        Subplot_A.set_aspect(1)
        Subplot_A.set_title(self.Left_Panel_sim_name, fontsize = self.fontsize)
        Subplot_A.set_ylabel(r"$y\,[M_{\text{ADM}}]$", fontsize = self.fontsize)
        Subplot_A.set_xlabel(r"$x\,[M_{\text{ADM}}]$", fontsize = self.fontsize)
        
        Subplot_B.set_aspect(1)
        Subplot_B.set_title(self.Right_panel_sim_name, fontsize = self.fontsize)
        Subplot_B.set_xlabel(r"$x\,[M_{\text{ADM}}]$", fontsize = self.fontsize)

        """ ============================================================================================================================================================================================================ """
        
        Panel_A_Intensity = self.Panel_A_context.Redshift**4 * (self.Panel_A_context.Pol_vec_x**2 + self.Panel_A_context.Pol_vec_y**2)
        Panel_B_Intensity = self.Panel_B_context.Redshift**4 * (self.Panel_B_context.Pol_vec_x**2 + self.Panel_B_context.Pol_vec_y**2)
        
        Max_X_range = max(max(self.Panel_A_context.X_coords.flatten()), max(self.Panel_B_context.X_coords.flatten()))
        Min_X_range = min(min(self.Panel_A_context.X_coords.flatten()), min(self.Panel_B_context.X_coords.flatten()))
            
        if not Pattern_only:
            
            Pol_intensity_plot = Main_Figure.add_subplot(143)
        
            Panel_A_Slice_idx = int(self.Panel_A_context.Y_res / 2) - 1
            Pol_intensity_plot.plot(self.Panel_A_context.X_coords[Panel_A_Slice_idx], Panel_A_Intensity[Panel_A_Slice_idx], color = "C0")
            
            Panel_B_Slice_idx = int(self.Panel_B_context.Y_res / 2) - 1
            Pol_intensity_plot.plot(self.Panel_B_context.X_coords[Panel_B_Slice_idx], Panel_B_Intensity[Panel_B_Slice_idx], color = "m")
    
            Intensity_max = max(max(Panel_A_Intensity.flatten()), max(Panel_B_Intensity.flatten()))
            Intensity_min = min(min(Panel_A_Intensity.flatten()), min(Panel_B_Intensity.flatten()))

            Pol_intensity_plot.set_aspect((Max_X_range - Min_X_range) / (1.1 * (Intensity_max - Intensity_min)))
            Pol_intensity_plot.set_title(r"Intensity at $y = 0$ [-]", fontsize = self.fontsize)
            Pol_intensity_plot.set_xlabel(r"$x\,[M_{\text{ADM}}]$", fontsize = self.fontsize)
            Pol_intensity_plot.set_ylim(0, 1.1 * Intensity_max)
            Pol_intensity_plot.set_xlim(float(Min_X_range), float(Max_X_range))
        
            Pol_intensity_plot.minorticks_on()
            Pol_intensity_plot.tick_params(axis = "both", direction = "out")
            Pol_intensity_plot.tick_params(which = 'minor', length = 4, labelsize = self.fontsize)
            Pol_intensity_plot.tick_params(which = 'major', length = 8, labelsize = self.fontsize) 
                
            Pol_intensity_plot.legend([self.Left_Panel_sim_name, self.Right_panel_sim_name], loc = "upper left", fontsize = self.fontsize - 6)
        
        """ ============================================================================================================================================================================================================ """
        
        Colorbar_norm = colors.Normalize(vmin = float(min(Panel_A_Intensity.flatten())), vmax = float(max(Panel_A_Intensity.flatten())))
        Colormap = colormaps[self.Colormap]
        Colormap.set_bad(self.Colormap_bad_color)

        for _, (Sim_context, Subplot) in enumerate(zip([self.Panel_A_context, self.Panel_B_context], 
                                                       [Subplot_A, Subplot_B])):
            
            Intensity = Sim_context.Redshift**4 * (Sim_context.Pol_vec_x**2 + Sim_context.Pol_vec_y**2)
            
            """ Set the empty pixels to nan, so the they get colored with the "bad" color """
            Intensity[Intensity == 0] = nan
            
            Disk_image = Subplot.imshow(Intensity, extent = tuple(Sim_context.Axes_limits), cmap = Colormap, interpolation = "nearest", norm = Colorbar_norm)
            
            Subplot.set_xlim((Sim_context.Axes_limits[0], Sim_context.Axes_limits[1]))
            Subplot.set_ylim((Sim_context.Axes_limits[2], Sim_context.Axes_limits[3]))
            
            Subplot.minorticks_on()
            Subplot.tick_params(axis = "both", direction = "out")
            Subplot.tick_params(which = 'minor', length = 4, labelsize = self.fontsize)
            Subplot.tick_params(which = 'major', length = 8, labelsize = self.fontsize)
            
            X_coords_to_plot = []
            Y_coords_to_plot = []
            
            Pol_vec_x_to_plot = []
            Pol_vec_y_to_plot = []
            
            Pol_tick_scale = Scale_factor / sqrt(Intensity)
            
            for x_idx in range(0, Sim_context.X_res, Sim_context.Pixel_skip):
                
                for y_idx in range(0, Sim_context.Y_res, Sim_context.Pixel_skip):
                    
                    X_coords_to_plot.append(Sim_context.X_coords[x_idx][y_idx] - Pol_tick_scale[x_idx][y_idx] * Sim_context.Redshift[x_idx][y_idx]**2 * Sim_context.Pol_vec_x[x_idx][y_idx] / 2)
                    Y_coords_to_plot.append(Sim_context.Y_coords[x_idx][y_idx] - Pol_tick_scale[x_idx][y_idx] * Sim_context.Redshift[x_idx][y_idx]**2 * Sim_context.Pol_vec_y[x_idx][y_idx] / 2)
                    
                    Pol_vec_x_to_plot.append(Pol_tick_scale[x_idx][y_idx] * Sim_context.Redshift[x_idx][y_idx]**2 * Sim_context.Pol_vec_x[x_idx][y_idx])
                    Pol_vec_y_to_plot.append(Pol_tick_scale[x_idx][y_idx] * Sim_context.Redshift[x_idx][y_idx]**2 * Sim_context.Pol_vec_y[x_idx][y_idx])
                    
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
            
        colorbar = Main_Figure.colorbar(Disk_image, ax = [Subplot_A, Subplot_B], orientation = "horizontal", extend = "both", aspect = 50) # type: ignore
        colorbar.set_label(r"Intensity [-]", fontsize = self.fontsize)
        colorbar.ax.tick_params(labelsize = self.fontsize)

        if not Pattern_only:

            """ ============================================================================================================================================================================================================ """

            EVPA_plot = Main_Figure.add_subplot(144)
    
            (Panel_A_Left_Pol_x_slice, Panel_A_Right_Pol_x_slice, 
            Panel_A_Left_x_slice, Panel_A_Right_x_slice) = self.Split_curve_across_brightness_depression(self.Panel_A_context.Pol_vec_x[Panel_A_Slice_idx], 
                                                                                                            self.Panel_A_context.X_coords[Panel_A_Slice_idx])
        
            
            Panel_A_Left_Pol_y_slice, Panel_A_Right_Pol_y_slice, _, _ = self.Split_curve_across_brightness_depression(self.Panel_A_context.Pol_vec_y[Panel_A_Slice_idx],  
                                                                                                                        self.Panel_A_context.X_coords[Panel_A_Slice_idx])
        
            Panel_A_Left_EVPA = arctan(-Panel_A_Left_Pol_x_slice / Panel_A_Left_Pol_y_slice)
            Panel_A_Right_EVPA = arctan(-Panel_A_Right_Pol_x_slice / Panel_A_Right_Pol_y_slice)
            
            (Panel_B_Left_Pol_x_slice, Panel_B_Right_Pol_x_slice, 
            Panel_B_Left_x_slice, Panel_B_Right_x_slice) = self.Split_curve_across_brightness_depression(self.Panel_B_context.Pol_vec_x[Panel_B_Slice_idx], 
                                                                                                                self.Panel_B_context.X_coords[Panel_B_Slice_idx])
        
            
            Left_Lower_Panel_Pol_y_slice, Right_Lower_Panel_Pol_y_slice, _, _ = self.Split_curve_across_brightness_depression(self.Panel_B_context.Pol_vec_y[Panel_B_Slice_idx],  
                                                                                                                            self.Panel_B_context.X_coords[Panel_B_Slice_idx])
        
            Panel_B_Left_EVPA = arctan(-Panel_B_Left_Pol_x_slice / Left_Lower_Panel_Pol_y_slice)
            Panel_B_Right_EVPA = arctan(-Panel_B_Right_Pol_x_slice / Right_Lower_Panel_Pol_y_slice)
            
            """ ====================================================================================================== """

            EVPA_plot.set_title(r"$\text{EVPA at}\,y = 0\,[\text{rad}]$", fontsize = self.fontsize)
            EVPA_plot.set_xlabel(r"$x\,[M_{\text{ADM}}]$", fontsize = self.fontsize)

            EVPA_plot.plot([max(Panel_A_Left_x_slice), max(Panel_A_Left_x_slice)], [-10, 10], "k")
            EVPA_plot.plot([min(Panel_A_Right_x_slice), min(Panel_A_Right_x_slice)], [-10, 10], "k")
            
            EVPA_plot.plot(Panel_A_Left_x_slice, Panel_A_Left_EVPA, color = "m")
            EVPA_plot.plot(Panel_A_Right_x_slice, Panel_A_Right_EVPA, color = "m")

            EVPA_plot.plot([max(Panel_B_Left_x_slice), max(Panel_B_Left_x_slice)], [-10, 10], "k--")
            EVPA_plot.plot([min(Panel_B_Right_x_slice), min(Panel_B_Right_x_slice)], [-10, 10], "k--")
            
            EVPA_plot.plot(Panel_B_Left_x_slice, Panel_B_Left_EVPA, "C0")
            EVPA_plot.plot(Panel_B_Right_x_slice, Panel_B_Right_EVPA, "C0")
            
            Numerical_EVPA_max = max(max(Panel_A_Left_EVPA.flatten()), max(Panel_A_Right_EVPA.flatten()))
            Numerical_EVPA_min = min(min(Panel_A_Left_EVPA.flatten()), min(Panel_A_Right_EVPA.flatten()))
            
            Lower_Panel_EVPA_max = max(max(Panel_B_Left_EVPA.flatten()), max(Panel_B_Right_EVPA.flatten()))
            Lower_Panel_EVPA_min = min(min(Panel_B_Left_EVPA.flatten()), min(Panel_B_Right_EVPA.flatten()))
            
            EVPA_max = max(Lower_Panel_EVPA_max, Numerical_EVPA_max)
            EVPA_min = min(Lower_Panel_EVPA_min, Numerical_EVPA_min)
            
            if EVPA_min < 0:
                EVPA_plot.set_ylim(1.1 * EVPA_min, 1.1 * EVPA_max)
                EVPA_plot.set_aspect((Max_X_range - Min_X_range) / (1.1 * (EVPA_max - EVPA_min)))
                
            else:
                EVPA_plot.set_ylim(0.9 * EVPA_min, 1.1 * EVPA_max)
                EVPA_plot.set_aspect((Max_X_range - Min_X_range) / (1.1 * EVPA_max - 0.9 * EVPA_min))      
                
            EVPA_plot.set_xlim(float(Min_X_range), float(Max_X_range))
                        
            EVPA_plot.minorticks_on()
            EVPA_plot.tick_params(axis = "both", direction = "out")
            EVPA_plot.tick_params(which = 'minor', length = 4, labelsize = self.fontsize)
            EVPA_plot.tick_params(which = 'major', length = 8, labelsize = self.fontsize) 
        
        if not os.path.exists(self.Output_path):
                os.makedirs(self.Output_path)

        Main_Figure.savefig(self.Output_path + "\\Whole_disk_plot_B_{}_{}_{}.pdf".format(self.Panel_A_context.B_field[0],
                                                                                         self.Panel_A_context.B_field[1],
                                                                                         self.Panel_A_context.B_field[2]), 
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
         
        """ ====================================== Lower_Panel image extraction ====================================== """
    
        X_coords = self.Left_Upper_Panel_Parser.X_coords[abs(self.Left_Upper_Panel_Parser.Source_r - r_source) < tolerance]
        Y_coords = self.Left_Upper_Panel_Parser.Y_coords[abs(self.Left_Upper_Panel_Parser.Source_r - r_source) < tolerance]

        Pol_x = self.Left_Upper_Panel_Parser.Polarization_vec_X[abs(self.Left_Upper_Panel_Parser.Source_r - r_source) < tolerance]
        Pol_y = self.Left_Upper_Panel_Parser.Polarization_vec_Y[abs(self.Left_Upper_Panel_Parser.Source_r - r_source) < tolerance]
        
        Redshift = self.Left_Upper_Panel_Parser.Disk_redshift[abs(self.Left_Upper_Panel_Parser.Source_r - r_source) < tolerance]

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
        
        Intensity = self.Polarized_Intensity_data[self.Panel_Enums.Upper.value]
        Intensity[Intensity == 0] = nan
        
        Axes_limits = self.Left_Upper_Panel_Parser.Simulation_metadata["Observation Window Dimentions (-X,+X,-Y,+Y) [M]"].split(",")
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

    def __get_image_at_fixed_orbit_radius(self, source_r_orbit: float, target_r_orbit: float, source_tolerance: float, target_tolerance: float, num_points: int, Curve_visualization: bool, 
                                          Source_center_Y_offset: float, Target_center_Y_offset: float, Bounding_box_vertecies: list[float]):   
        
        Source_Parser = self.Panel_B_parser
        Target_Parser = self.Panel_A_parser
        
        """ Some metric have inner images that I want to filter out. """
        Outside_bounding_box_X = logical_or(Source_Parser.X_coords / self.Panel_B_context.ADM_mass < Bounding_box_vertecies[0], Source_Parser.X_coords / self.Panel_B_context.ADM_mass > Bounding_box_vertecies[1])
        Outside_bounding_box_Y = logical_or(Source_Parser.Y_coords / self.Panel_B_context.ADM_mass < Bounding_box_vertecies[2], Source_Parser.Y_coords / self.Panel_B_context.ADM_mass > Bounding_box_vertecies[3])
                
        Outside_bounding_box = logical_or(Outside_bounding_box_X, Outside_bounding_box_Y)
        Source_radius_condition = abs(Source_Parser.Source_r - source_r_orbit) < source_tolerance * (1 + abs(arctan2((Source_Parser.Y_coords - Source_center_Y_offset * self.Panel_B_context.ADM_mass) , Source_Parser.X_coords))**2)
        
        Filter_condition = logical_and(logical_and(Source_radius_condition, Outside_bounding_box), Source_Parser.Source_r > 0)

        """ ====================================== Source Panel image extraction ====================================== """
    
        Source_X_coords = Source_Parser.X_coords[Filter_condition]
        Source_Y_coords = Source_Parser.Y_coords[Filter_condition]

        Source_Pol_x = Source_Parser.Polarization_vec_X[Filter_condition]
        Source_Pol_y = Source_Parser.Polarization_vec_Y[Filter_condition]
        
        Source_redshift = Source_Parser.Disk_redshift[Filter_condition]

        """ Offset the Y position of the image, so the y = 0 line passes trough the image center - this makes sorting out points based on image azimuth easier. """
        Source_Y_coords = Source_Y_coords - Source_center_Y_offset * self.Panel_B_context.ADM_mass

        Source_Image_azimuth = arctan2(Source_Y_coords, Source_X_coords)
        
        Source_Image_azimuth[Source_Image_azimuth < 0] = Source_Image_azimuth[Source_Image_azimuth < 0] + 2 * pi
        Sorting_idx = argsort(Source_Image_azimuth)
        
        Source_Image_azimuth = Source_Image_azimuth[Sorting_idx]
        Source_X_coords = Source_X_coords[Sorting_idx]
        Source_Y_coords = Source_Y_coords[Sorting_idx]
        Source_redshift = Source_redshift[Sorting_idx]
        Source_Pol_x = Source_Pol_x[Sorting_idx]
        Source_Pol_y = Source_Pol_y[Sorting_idx]

        Total_length = 0
        Current_length = 0
        Final_idx_list = [0]

        for idx, _ in enumerate(Source_Image_azimuth):
            
            if idx > 0:
                Total_length = Total_length + sqrt(Source_X_coords[idx]**2 + Source_Y_coords[idx]**2) * abs(Source_Image_azimuth[idx] - Source_Image_azimuth[idx - 1])
            
        for idx, _ in enumerate(Source_Image_azimuth):
            
            if idx > 0:
                Current_length = Current_length + sqrt(Source_X_coords[idx]**2 + Source_Y_coords[idx]**2) * abs(Source_Image_azimuth[idx] - Source_Image_azimuth[idx - 1])
            
            if Current_length >= Total_length / num_points or Curve_visualization:
                Current_length = 0
                Final_idx_list.append(idx)
        
        Source_X_coords = Source_X_coords[Final_idx_list]  
        Source_Y_coords = Source_Y_coords[Final_idx_list]       
        Source_Pol_x = Source_Pol_x[Final_idx_list]
        Source_Pol_y = Source_Pol_y[Final_idx_list]
        Source_redshift = Source_redshift[Final_idx_list]
        Source_Image_azimuth = Source_Image_azimuth[Final_idx_list]
        
        """ Reverse the Y offset, so I can plot the image and have it make sense. """
        Source_Y_coords = Source_Y_coords + Source_center_Y_offset * self.Panel_B_context.ADM_mass
        
        """ ====================================== Target Panel image extraction ====================================== """

        Source_radius_condition = abs(Target_Parser.Source_r - target_r_orbit) < target_tolerance * (1 + abs(arctan2((Target_Parser.Y_coords - Target_center_Y_offset * self.Panel_A_context.ADM_mass) , Target_Parser.X_coords))**2)
        Filter_condition = logical_and(Source_radius_condition, Target_Parser.Source_r > 0)
        
        Target_X_coords = Target_Parser.X_coords[Filter_condition]
        Target_Y_coords = Target_Parser.Y_coords[Filter_condition]
        
        Target_Pol_x = Target_Parser.Polarization_vec_X[Filter_condition]
        Target_Pol_y = Target_Parser.Polarization_vec_Y[Filter_condition]
        
        Target_redshift = Target_Parser.Disk_redshift[Filter_condition]

        Target_Y_coords = Target_Y_coords - Target_center_Y_offset * self.Panel_A_context.ADM_mass

        Target_Image_azimuth = arctan2(Target_Y_coords, Target_X_coords)
        Target_Image_azimuth[Target_Image_azimuth < 0] = Target_Image_azimuth[Target_Image_azimuth < 0] + 2 * pi
        Sorting_idx = argsort(Target_Image_azimuth)
               
        Target_Image_azimuth = Target_Image_azimuth[Sorting_idx]
        Target_X_coords = Target_X_coords[Sorting_idx]
        Target_Y_coords = Target_Y_coords[Sorting_idx]
        Target_redshift = Target_redshift[Sorting_idx]
        Target_Pol_x = Target_Pol_x[Sorting_idx]
        Target_Pol_y = Target_Pol_y[Sorting_idx]
  
        Return_Target_X_coords = []
        Return_Target_Y_coords = []
        Return_Target_Pol_x = []
        Return_Target_Pol_y = []
        Return_Target_redshift = []
        
        Temp_Target_Azimuth = []

        Total_length = 0
        Current_length = 0
        Final_idx_list = [0]
        
        if Curve_visualization:
        
            for _, Current_azimuth in enumerate(Source_Image_azimuth):
                
                min_idx = argmin(abs(Target_Image_azimuth - Current_azimuth))

                Return_Target_X_coords.append(Target_X_coords[min_idx])
                Return_Target_Y_coords.append(Target_Y_coords[min_idx])
                Return_Target_Pol_x.append(Target_Pol_x[min_idx])
                Return_Target_Pol_y.append(Target_Pol_y[min_idx])
                Return_Target_redshift.append(Target_redshift[min_idx])
                Temp_Target_Azimuth.append(Current_azimuth)
                
            Return_Target_redshift = array(Return_Target_redshift)
            Return_Target_Pol_x = array(Return_Target_Pol_x)
            Return_Target_Pol_y = array(Return_Target_Pol_y)

        else:
            
            Total_length = 0
            Current_length = 0
            Final_idx_list = []
            
            for idx, _ in enumerate(Target_Image_azimuth):
            
                if idx > 0:
                    Total_length = Total_length + sqrt(Target_X_coords[idx]**2 + Target_Y_coords[idx]**2) * abs(Target_Image_azimuth[idx] - Target_Image_azimuth[idx - 1])
                
            for idx, _ in enumerate(Target_Image_azimuth):
                
                if idx > 0:
                    Current_length = Current_length + sqrt(Target_X_coords[idx]**2 + Target_Y_coords[idx]**2) * abs(Target_Image_azimuth[idx] - Target_Image_azimuth[idx - 1])
                
                if Current_length >= Total_length / num_points:
                    Current_length = 0
                    Final_idx_list.append(idx)
            
            Return_Target_X_coords = Target_X_coords[Final_idx_list]  
            Return_Target_Y_coords = Target_Y_coords[Final_idx_list]       
            Return_Target_Pol_x = Target_Pol_x[Final_idx_list]
            Return_Target_Pol_y = Target_Pol_y[Final_idx_list]
            Return_Target_redshift = Target_redshift[Final_idx_list]
                
        Return_Target_Y_coords = array(Return_Target_Y_coords) + Target_center_Y_offset * self.Panel_A_context.ADM_mass
        
        """ Normalize the image coordinates by the ADM mass. """
        Return_Target_X_coords = array(Return_Target_X_coords) / self.Panel_A_context.ADM_mass
        Return_Target_Y_coords = Return_Target_Y_coords / self.Panel_A_context.ADM_mass
        
        Source_X_coords = array(Source_X_coords) / self.Panel_B_context.ADM_mass
        Source_Y_coords = array(Source_Y_coords) / self.Panel_B_context.ADM_mass
        
        return Source_Image_azimuth, Source_Pol_x, Source_Pol_y, Source_redshift, Source_X_coords, Source_Y_coords, Return_Target_Pol_x, Return_Target_Pol_y, Return_Target_redshift, Return_Target_X_coords, Return_Target_Y_coords

    def __get_image_at_fixed_coords(self, r_source: float, tolerance: float, num_points: int, Curve_visualization: bool, Source_Panel: str):

        Source_Parser = self.Lower_Right_Panel_Parser
        Target_Parser = self.Left_Upper_Panel_Parser

        if Source_Panel == "Left/Upper":
            
            Source_Parser = self.Left_Upper_Panel_Parser
            Target_Parser = self.Lower_Right_Panel_Parser

        """ ====================================== Lower_Panel image extraction ====================================== """
    
        Source_X_coords = Source_Parser.X_coords[abs(Source_Parser.Source_r - r_source) < tolerance]
        Source_Y_coords = Source_Parser.Y_coords[abs(Source_Parser.Source_r - r_source) < tolerance]

        Source_Pol_x = Source_Parser.Polarization_vec_X[abs(Source_Parser.Source_r - r_source) < tolerance]
        Source_Pol_y = Source_Parser.Polarization_vec_Y[abs(Source_Parser.Source_r - r_source) < tolerance]
        
        Source_redshift = Source_Parser.Disk_redshift[abs(Source_Parser.Source_r - r_source) < tolerance]

        Image_azimuth = arctan2(Source_Y_coords, Source_X_coords)
        
        Image_azimuth[Image_azimuth < 0] = Image_azimuth[Image_azimuth < 0] + 2 * pi
        Sorting_idx = argsort(Image_azimuth)
        
        Image_azimuth = Image_azimuth[Sorting_idx]
        Source_X_coords = Source_X_coords[Sorting_idx]
        Source_Y_coords = Source_Y_coords[Sorting_idx]
        Source_redshift = Source_redshift[Sorting_idx]
        Source_Pol_x = Source_Pol_x[Sorting_idx]
        Source_Pol_y = Source_Pol_y[Sorting_idx]

        Total_length = 0
        Current_length = 0
        Final_idx_list = [0]

        for idx, _ in enumerate(Image_azimuth):
            
            if idx > 0:
                Total_length = Total_length + sqrt(Source_X_coords[idx]**2 + Source_Y_coords[idx]**2) * abs(Image_azimuth[idx] - Image_azimuth[idx - 1])
            
        for idx, _ in enumerate(Image_azimuth):
            
            if idx > 0:
                Current_length = Current_length + sqrt(Source_X_coords[idx]**2 + Source_Y_coords[idx]**2) * abs(Image_azimuth[idx] - Image_azimuth[idx - 1])
            
            if Current_length >= Total_length / num_points or Curve_visualization:
                Current_length = 0
                Final_idx_list.append(idx)
               
        Source_Pol_x = Source_Pol_x[Final_idx_list]
        Source_Pol_y = Source_Pol_y[Final_idx_list]
        Source_redshift = Source_redshift[Final_idx_list]
        Image_azimuth = Image_azimuth[Final_idx_list]
                
        """ ============================================ Numerical image extraction ============================================ """
    
        Source_X_coords = Source_X_coords[Final_idx_list]
        Source_Y_coords = Source_Y_coords[Final_idx_list]
    
        Target_X_grid = Target_Parser.X_coords[0:self.X_resolution]
        Target_Y_grid = Target_Parser.Y_coords[::self.Y_resolution]

        Target_X_coords = []
        Target_Y_coords = []
        Target_Pol_x = []
        Target_Pol_y = []
        Target_redshift = []
        
        for X_target, Y_target in zip(Source_X_coords, Source_Y_coords):
        
            X_idx = argmin(abs(Target_X_grid - X_target))
            Y_idx = argmin(abs(Target_Y_grid - Y_target))

            Target_X_coords.append(Target_Parser.X_coords[X_idx + self.Y_resolution * Y_idx])
            Target_Y_coords.append(Target_Parser.Y_coords[X_idx + self.Y_resolution * Y_idx])
            Target_Pol_x.append(Target_Parser.Polarization_vec_X[X_idx + self.Y_resolution * Y_idx])
            Target_Pol_y.append(Target_Parser.Polarization_vec_Y[X_idx + self.Y_resolution * Y_idx])
            Target_redshift.append(Target_Parser.Disk_redshift[X_idx + self.Y_resolution * Y_idx])

        Target_X_coords = array(Target_X_coords) / self.Solution_mass
        Target_Y_coords = array(Target_Y_coords) / self.Solution_mass
        
        Source_X_coords = array(Source_X_coords) / self.Solution_mass
        Source_Y_coords = array(Source_Y_coords) / self.Solution_mass
        
        Target_redshift = array(Target_redshift)
        Target_Pol_x = array(Target_Pol_x)
        Target_Pol_y = array(Target_Pol_y)
        
        if Source_Panel == "Left/Upper":
            
            Source_Parser = self.Left_Upper_Panel_Parser
            Target_Parser = self.Lower_Right_Panel_Parser
        
            return Image_azimuth, Source_Pol_x, Source_Pol_y, Source_redshift, Source_X_coords, Source_Y_coords, Target_Pol_x, Target_Pol_y, Target_redshift, Target_X_coords, Target_Y_coords
        
        return Image_azimuth, Target_Pol_x, Target_Pol_y, Target_redshift, Target_X_coords, Target_Y_coords, Source_Pol_x, Source_Pol_y, Source_redshift, Source_X_coords, Source_Y_coords

    def get_direct_image_at_fixed_source_radius(self, r_source: float, source_tolerance: float, target_tolerance: float, num_points: int, Scale_factor: float, axis_limits_upper: list[float], axis_limits_lower: list[float], Tick_visualization: bool, 
                                                Curve_visualization: bool, Subplots: list, Source_coords: str = "Lower_Panel", Plot_disk: bool = True,
                                                Source_center_Y_offset: float = 0, Target_center_Y_offset: float = 0, Bounding_box_vertecies: list[float] = [0, 0, 0, 0]) -> None:
         
        # if Source_coords == "Lower_Panel":
        #     (Image_azimuth, Top_Panel_Pol_x, Top_Panel_Pol_y, Top_Panel_redshift, Top_Panel_X_coords, Top_Panel_Y_coords, 
        #      Bottom_Panel_Pol_x, Bottom_Panel_Pol_y, Bottom_Panel_redshift, Bottom_Panel_X_coords, Bottom_Panel_Y_coords) = self.__get_image_at_fixed_coords(r_source, tolerance_1, num_points, Curve_visualization, Source_Panel = "")
            
        # else:
        #     (Image_azimuth, Top_Panel_Pol_x, Top_Panel_Pol_y, Top_Panel_redshift, Top_Panel_X_coords, Top_Panel_Y_coords, 
        #      Bottom_Panel_Pol_x, Bottom_Panel_Pol_y, Bottom_Panel_redshift, Bottom_Panel_X_coords, Bottom_Panel_Y_coords) = self.__get_image_at_fixed_coords(r_source, tolerance_1, num_points, Curve_visualization, Source_Panel = "Left/Upper")

        (Image_azimuth, Bottom_Panel_Pol_x, Bottom_Panel_Pol_y, Bottom_Panel_redshift, Bottom_Panel_X_coords, Bottom_Panel_Y_coords,
         Top_Panel_Pol_x, Top_Panel_Pol_y, Top_Panel_redshift, Top_Panel_X_coords, Top_Panel_Y_coords) = self.__get_image_at_fixed_orbit_radius(r_source, r_source, source_tolerance, target_tolerance, num_points, Curve_visualization, Source_center_Y_offset, Target_center_Y_offset, Bounding_box_vertecies)

        Image_azimuth = Image_azimuth - Image_azimuth[0]
        Image_azimuth = append(Image_azimuth, Image_azimuth + 2 * pi)

        Panel_A_Intensity = self.Panel_A_context.Redshift**4 * (self.Panel_A_context.Pol_vec_x**2 + self.Panel_A_context.Pol_vec_y**2)
        Panel_B_Intensity = self.Panel_B_context.Redshift**4 * (self.Panel_B_context.Pol_vec_x**2 + self.Panel_B_context.Pol_vec_y**2)
        
        if Tick_visualization:

            """ ==================================================== Lower_Panel Plotting ==================================================== """

            Colormap = colormaps[self.Colormap]
            Colormap.set_bad(self.Colormap_bad_color)
                
            if Subplots[1] != None:

                Panel_B_Intensity[Panel_B_Intensity == 0] = nan
                
                Bottom_Panel_axes_limits = self.Panel_B_parser.Simulation_metadata["Observation Window Dimentions (-X,+X,-Y,+Y) [M]"].split(",")
                Bottom_Panel_axes_limits = array([float(Limit) / self.Panel_B_context.ADM_mass for Limit in Bottom_Panel_axes_limits])
                
                if Plot_disk:
                    Subplots[1].imshow(Panel_B_Intensity, extent = tuple(Bottom_Panel_axes_limits), cmap = Colormap, interpolation = "nearest", vmin = 0)
        
                Final_tick_scale = Scale_factor / sqrt(Bottom_Panel_Pol_x**2 + Bottom_Panel_Pol_y**2 + 1e-10)
                Subplots[1].quiver(Bottom_Panel_X_coords - Final_tick_scale * Bottom_Panel_Pol_x / 2,
                                    Bottom_Panel_Y_coords - Final_tick_scale * Bottom_Panel_Pol_y / 2,
                                    Final_tick_scale * Bottom_Panel_Pol_x,
                                    Final_tick_scale * Bottom_Panel_Pol_y,
                                    headwidth = 0,
                                    headlength = 0,
                                    headaxislength = 0,
                                    angles = 'xy', 
                                    scale_units = 'xy',
                                    scale = 1,
                                    color = "grey",
                                    width = 0.005)
                
                Subplots[1].set_xlim(axis_limits_lower[0], axis_limits_lower[1])
                Subplots[1].set_ylim(axis_limits_lower[2], axis_limits_lower[3])
                
            """ ==================================================== Numerical Plotting ==================================================== """
   
            if Subplots[0] != None:
    
                Panel_A_Intensity[Panel_A_Intensity == 0] = nan
                
                Top_Panel_axes_limits = self.Panel_A_parser.Simulation_metadata["Observation Window Dimentions (-X,+X,-Y,+Y) [M]"].split(",")
                Top_Panel_axes_limits = array([float(Limit) / self.Panel_A_context.ADM_mass for Limit in Top_Panel_axes_limits])
                
                if Plot_disk:
                    Subplots[0].imshow(Panel_A_Intensity, extent = tuple(Top_Panel_axes_limits), cmap = Colormap, interpolation = "nearest", vmin = 0)
                
                Final_tick_scale = Scale_factor / sqrt(Top_Panel_Pol_x**2 + Top_Panel_Pol_y**2 + 1e-10)
                Subplots[0].quiver(Top_Panel_X_coords - Final_tick_scale * Top_Panel_Pol_x / 2,
                                    Top_Panel_Y_coords - Final_tick_scale * Top_Panel_Pol_y / 2,
                                    Final_tick_scale * Top_Panel_Pol_x,
                                    Final_tick_scale * Top_Panel_Pol_y,
                                    headwidth = 0,
                                    headlength = 0,
                                    headaxislength = 0,
                                    angles = 'xy', 
                                    scale_units = 'xy',
                                    scale = 1,
                                    color = "grey",
                                    width = 0.005)
                
                Subplots[0].set_xlim(axis_limits_upper[0], axis_limits_upper[1])
                Subplots[0].set_ylim(axis_limits_upper[2], axis_limits_upper[3])
        
        if Curve_visualization:
            
            if Subplots[1] != None:
    
                Bottom_Panel_X_coords = append(Bottom_Panel_X_coords, Bottom_Panel_X_coords[0])
                Bottom_Panel_Y_coords = append(Bottom_Panel_Y_coords, Bottom_Panel_Y_coords[0])
    
                Subplots[1].plot(Bottom_Panel_X_coords, Bottom_Panel_Y_coords, "r--")
                Subplots[1].set_ylabel(r"y [$M_\text{ADM}$]", fontsize = self.fontsize)
                Subplots[1].set_xlabel(r"x [$M_\text{ADM}$]", fontsize = self.fontsize)
                Subplots[1].set_title(self.Right_panel_sim_name, usetex=True, fontsize = self.fontsize)
                
                Subplots[1].minorticks_on()
                Subplots[1].tick_params(axis = "both", direction = "out")
                Subplots[1].tick_params(which = 'minor', length = 4, labelsize = self.fontsize)
                Subplots[1].tick_params(which = 'major', length = 8, labelsize = self.fontsize) 
                
                Subplots[1].xaxis.set_major_locator(MaxNLocator(nbins = 5))
                Subplots[1].yaxis.set_major_locator(MaxNLocator(nbins = 5))

            if Subplots[0] != None:

                Top_Panel_X_coords = append(Top_Panel_X_coords, Top_Panel_X_coords[0])
                Top_Panel_Y_coords = append(Top_Panel_Y_coords, Top_Panel_Y_coords[0])

                Subplots[0].plot(Top_Panel_X_coords, Top_Panel_Y_coords, "r--")
                Subplots[0].set_ylabel(r"y [$M_\text{ADM}$]", fontsize = self.fontsize)
                Subplots[0].set_xlabel(r"x [$M_\text{ADM}$]", fontsize = self.fontsize)
                Subplots[0].set_title(self.Left_Panel_sim_name, fontsize = self.fontsize)
                
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
                
                Top_Panel_Intensity = Top_Panel_redshift**4 * (Top_Panel_Pol_x**2 + Top_Panel_Pol_y**2)
                Top_Panel_Intensity = append(Top_Panel_Intensity, Top_Panel_Intensity)
                
                Bottom_Panel_Intensity = Bottom_Panel_redshift**4 * (Bottom_Panel_Pol_x**2 + Bottom_Panel_Pol_y**2)
                Bottom_Panel_Intensity = append(Bottom_Panel_Intensity, Bottom_Panel_Intensity)
                
                Subplots[2].plot(Image_azimuth, Top_Panel_Intensity, "m")
                Subplots[2].plot(Image_azimuth, Bottom_Panel_Intensity, "C0")
                
                Intensity_max = max(max(Bottom_Panel_Intensity), max(Top_Panel_Intensity))
                
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
                Subplots[2].legend([self.Left_Panel_sim_name, self.Right_panel_sim_name], fontsize = self.fontsize - 10, loc = 'upper center', ncols = 2)
                
                Subplots[2].set_xticks(x_ticks)
                Subplots[2].set_xticklabels(x_tick_labels)
                
                Subplots[2].ticklabel_format(style = 'sci', axis = 'y', scilimits=(0, 0))
                Subplots[2].yaxis.set_major_locator(MaxNLocator(nbins = 6))
                
            """" ================================================================================================================================================== """
            
            if Subplots[2] != None:
                    
                Bottom_Panel_EVPA = arctan(-Bottom_Panel_Pol_x / Bottom_Panel_Pol_y)
                Bottom_Panel_EVPA = append(Bottom_Panel_EVPA, Bottom_Panel_EVPA)
                
                Top_Panel_EVPA = arctan(-Top_Panel_Pol_x / Top_Panel_Pol_y)
                Top_Panel_EVPA = append(Top_Panel_EVPA, Top_Panel_EVPA)
                
                Bottom_Panel_EVPA_idx_split = self.split_EVPA(EVPA = Bottom_Panel_EVPA)
                Top_Panel_EVPA_idx_split = self.split_EVPA(EVPA = Top_Panel_EVPA)
                
                for idx in range(len(Bottom_Panel_EVPA_idx_split) - 1):
                    Subplots[3].plot(Image_azimuth[Bottom_Panel_EVPA_idx_split[idx] : Bottom_Panel_EVPA_idx_split[idx + 1]] - Image_azimuth[0], 
                                     Bottom_Panel_EVPA[Bottom_Panel_EVPA_idx_split[idx] : Bottom_Panel_EVPA_idx_split[idx + 1]], "C0")    
                       
                    Subplots[3].plot(Image_azimuth[Bottom_Panel_EVPA_idx_split[idx + 1] - 1 : Bottom_Panel_EVPA_idx_split[idx + 1] + 1] - Image_azimuth[0], 
                                     Bottom_Panel_EVPA[Bottom_Panel_EVPA_idx_split[idx + 1] - 1 : Bottom_Panel_EVPA_idx_split[idx + 1] + 1], "C0--")               
                    
                for idx in range(len(Top_Panel_EVPA_idx_split) - 1):        
                    Subplots[3].plot(Image_azimuth[Top_Panel_EVPA_idx_split[idx] : Top_Panel_EVPA_idx_split[idx + 1]] - Image_azimuth[0], 
                                     Top_Panel_EVPA[Top_Panel_EVPA_idx_split[idx] : Top_Panel_EVPA_idx_split[idx + 1]], "m")
                    
                    Subplots[3].plot(Image_azimuth[Top_Panel_EVPA_idx_split[idx + 1] - 1 : Top_Panel_EVPA_idx_split[idx + 1] + 1] - Image_azimuth[0], 
                                     Top_Panel_EVPA[Top_Panel_EVPA_idx_split[idx + 1] - 1 : Top_Panel_EVPA_idx_split[idx + 1] + 1], "m--")
                    
                EVPA_max = max(max(Bottom_Panel_EVPA), max(Top_Panel_EVPA))
                EVPA_min = min(min(Bottom_Panel_EVPA), min(Top_Panel_EVPA))
                
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
                                        
                Top_Panel_Intensity = Top_Panel_redshift**4 * (Top_Panel_Pol_x**2 + Top_Panel_Pol_y**2)
                Top_Panel_Intensity = append(Top_Panel_Intensity, Top_Panel_Intensity)
                
                Bottom_Panel_Intensity = Bottom_Panel_redshift**4 * (Bottom_Panel_Pol_x**2 + Bottom_Panel_Pol_y**2)
                Bottom_Panel_Intensity = append(Bottom_Panel_Intensity, Bottom_Panel_Intensity)
                
                Subplots[4].plot(Image_azimuth, Top_Panel_Intensity - Bottom_Panel_Intensity)
                Delta_I_max = max(Top_Panel_Intensity - Bottom_Panel_Intensity)
                Delta_I_min = min(Top_Panel_Intensity - Bottom_Panel_Intensity)
                
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
                
                Delta_EVPA = arctan((Top_Panel_Pol_x * Bottom_Panel_Pol_y - Top_Panel_Pol_y * Bottom_Panel_Pol_x) / (Top_Panel_Pol_x * Bottom_Panel_Pol_x + Top_Panel_Pol_y * Bottom_Panel_Pol_y))
                Delta_EVPA = append(Delta_EVPA, Delta_EVPA)
                idx = argmax(abs(Delta_EVPA))
                
                print("X = {}".format(Bottom_Panel_X_coords[idx]))
                print("Y = {}".format(Bottom_Panel_Y_coords[idx]))
                
                # Subplots[0].plot(Bottom_Panel_X_coords[idx], Bottom_Panel_Y_coords[idx], "ro")
                
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
    
    #             Sim_path = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\Utilities\\Numerical_metric_runs\\Zero_curvature\\Config_II_70_deg_B_0.87_0.0_0.5_min_order_1_max_more_order_1_zoom".format(Model, inc, B_field[0], B_field[1], B_field[2])
                
    #             Visualizer_instance = Polarization_visuzlier(Sim_path = Sim_path, 
    #                                                         fontsize = 26, 
    #                                                         Fig_output_path = "E:\\Numerical_metric_runs\\Zero_curvature\\Figures\\Config_{}_{}_deg_min_order_0_max_order_0".format(Model, inc),
    #                                                         Numerical_title = r"Model $\mathbf{{{}}}_{{{}}}$".format(Model, horizon_radii[idx]))
                
    #             Visualizer_instance.Plot_polarization_ticks(Pixel_skip_step = 64, Scale_factor = 0.75)
                
    #             plt.show()
                
    Panel_A_path = "E:\\Numerical_metric_runs\\Zero_curvature\\Config_II_70_deg_B_0.0_1.0_0.0_min_order_0_max_order_0_zoom\\Numerical_results\\All_Segments_Results"
    Panel_B_path = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\Utilities\\Numerical_metric_runs\\Zero_curvature\\Config_II_70_deg_B_0.0_1.0_0.0_min_order_1_max_more_order_1_zoom\\Numerical_results\\All_Segments_Results"
                  
    Visualizer_instance = Polarization_visuzlier(Panel_A_path = Panel_A_path, 
                                                 Panel_B_path = Panel_B_path,
                                                 fontsize = 22, 
                                                 Fig_output_path = "E:\\Numerical_metric_runs\\Zero_curvature\\Figures\\Config_V_70_deg_min_order_0_max_order_0",
                                                 Left_Panel_sim_name = r"Direct Image",
                                                 Right_panel_sim_name = r"First Indirect Image")
    
    Visualizer_instance.Plot_polarization_ticks(Scale_factor = 0.5, Pattern_only = True)   
        
    Figure = plt.figure(figsize = (13, 10.5), layout = "compressed")
    
    Figure.suptitle(r"$i = {}$, $B_r = {}$, $B_\theta = {}$, $B_\phi = {}$".format(float(Visualizer_instance.Panel_A_parser.Simulation_metadata["Observer Inclination [Deg]"]),
                                                                                         Visualizer_instance.Panel_A_context.B_field[0],
                                                                                         Visualizer_instance.Panel_A_context.B_field[1],
                                                                                         Visualizer_instance.Panel_A_context.B_field[2]),
                    fontsize = Visualizer_instance.fontsize + 4,
                    bbox = dict(facecolor = 'none', edgecolor = 'black', boxstyle = 'round, pad = 0.3'),
                    y = 0.99)
    
    Numerical_subplot = Figure.add_subplot(231)
    Lower_Panel_subplot = Figure.add_subplot(234)
    Intensity_subplot = Figure.add_subplot(232)
    Delta_Intensity_subplot = Figure.add_subplot(235)
    EVPA_subplot = Figure.add_subplot(233)
    Delta_EVPA_subplot = Figure.add_subplot(236)
    
    Subplots = [Numerical_subplot, Lower_Panel_subplot, Intensity_subplot, EVPA_subplot, Delta_Intensity_subplot, Delta_EVPA_subplot]
    
    for orbit_mult in [1]:
        
        Visualizer_instance.get_direct_image_at_fixed_source_radius(r_source = orbit_mult * 0.054, 
                                                                    source_tolerance = 0.00008, 
                                                                    target_tolerance = 0.00008,
                                                                    num_points = 65, 
                                                                    Scale_factor = 1.25,  
                                                                    axis_limits_upper = [-1.75, 2.5, -0.25, 3.2],
                                                                    axis_limits_lower = [-4, 6, -4, 6],
                                                                    Tick_visualization = True,
                                                                    Curve_visualization = False,
                                                                    Subplots = Subplots,
                                                                    Source_center_Y_offset = -0.94,
                                                                    Target_center_Y_offset = 0.8,
                                                                    Bounding_box_vertecies = [-1.64, 2.07, 0.32, 2.75])
        
        Visualizer_instance.get_direct_image_at_fixed_source_radius(r_source = orbit_mult * 0.054, 
                                                                    source_tolerance = 0.000008, 
                                                                    target_tolerance = 0.00005,
                                                                    num_points = 65, 
                                                                    Scale_factor = 1.25, 
                                                                    axis_limits_upper = [-1.75, 2.5, -0.25, 3.2],
                                                                    axis_limits_lower = [-4, 6, -4, 6],
                                                                    Tick_visualization = False,
                                                                    Curve_visualization = True,
                                                                    Subplots = Subplots,
                                                                    Source_center_Y_offset = -0.94,
                                                                    Target_center_Y_offset = 0.8,
                                                                    Bounding_box_vertecies = [-1.64, 2.07, 0.32, 2.75])
    
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