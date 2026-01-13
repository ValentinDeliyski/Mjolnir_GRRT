import sys
import os

""" Add the parent directory of this file to the search path, 
    so this script can be ran from the "Utilities" folder """
    
parent_directory = os.path.abspath('...')
sys.path.append(parent_directory)

from Mjolnir_Configurator import Simulation_configurator
from Support_functions.Parsers import Units_class, Simulation_Parser
from Support_functions.Spacetimes_new import Kerr

from numpy import pi, tan, sqrt, linspace, arctan, array
from numpy.typing import NDArray

import subprocess

import matplotlib.pyplot as plt

class Simulation:
    
    def __init__(self):
        
        self.Sim_config = Simulation_configurator()
        
        self.Sim_config.simulation_mode = {"Value": 3, "Unit": "[-]"} 
        
        self.Sim_config.metric_parameters.Numerical_metric_spline_path = "C:/Users/Valur/Documents/Repos/Mjolnir_GRRT/Utilities/Galin_zero_curvature_config_I.XML"

        self.Sim_config.min_image_order = {"Value": 0, "Unit": "[-]"}
        self.Sim_config.max_image_order = {"Value": 10, "Unit": "[-]"}

        """ ================================================== The Numerical Metric ================================================== """

        self.ADM_Mass = 0.881990876889021
        self.Black_Hole_Mass = 0.0034981920971278923
        
        self.ADM_Ang_Momentum = 0.7258549181881021
        self.Black_Hole_Ang_Momentum = 9.145481718889068e-5
        
        self.r_ISCO = 0.268499434023205704
        
        self.Sim_config.metric_parameters.Metric_type    = {"Value": "Numerical", "Unit": "[-]"}
        
        self.Sim_config.metric_parameters.Mass           = {"Value": self.ADM_Mass, "Unit": "[M]"}
        self.Sim_config.metric_parameters.Horizon_radius = {"Value": 0.01, "Unit": "[G/c^2]"}
        self.Sim_config.metric_parameters.Spin           = {"Value": self.ADM_Ang_Momentum / self.ADM_Mass, "Unit": "[M]"}
        
        self.Sim_config.metric_parameters.Scattering_radius = {"Value": 300 * self.ADM_Mass, "Unit": "[M]"} 
        self.Sim_config.metric_parameters.Numerical_metric_anzatz_type = {"Value": "Anzatz_1", "Unit": "[-]"} 
        
        self.Sim_config.metric_parameters.Distance_to_singular_point = {"Value": 1e-3, "Unit": "[M]"}
        
        """ ================================================ Observer ================================================ """
   
        self.Sim_config.observer.Distance    = {"Value": 1e4 * self.ADM_Mass, "Unit": "[M]"}
        self.Sim_config.observer.Inclination = {"Value": 80 * pi / 180, "Unit": "[Rad]"}
        self.Sim_config.observer.Obs_frequency = {"Value": 230e9, "Unit": "[Hz]"}
        
        self.Sim_config.observer.Cam_rotation_angle = {"Value": 0, "Unit": "[Hz]"}
        self.Sim_config.observer.Use_angular_coords = {"Value": 1, "Unit": "[M]"}
        
        self.Nominal_Image_x_max = 15 * self.ADM_Mass
        self.Nominal_Image_y_max = 15 * self.ADM_Mass
        
        self.Nominal_Image_x_angle_max = arctan(self.Nominal_Image_x_max / self.Sim_config.observer.Distance["Value"])
        self.Nominal_Image_y_angle_max = arctan(self.Nominal_Image_y_max / self.Sim_config.observer.Distance["Value"])
         
        """ ================================================== Disk ================================================== """
        
        self.Sim_config.disk_model.Disk_Model = {"Value": "Novikov-Thorne", "Unit": "[-]"}
        
        self.Sim_config.disk_model.r_in_NT_disk = {"Value": self.r_ISCO, "Unit": "[M]"} 
        self.Sim_config.disk_model.r_out_NT_disk = {"Value": 5 * self.ADM_Mass, "Unit": "[M]"} 
        
        """ =============================================== Integrator =============================================== """
        
        """ The simulation output file path """
        self.Sim_config.file_manager.Output_file_directory = parent_directory + "Reference_simulations"

    def run_simulation(self):
        
        self.Sim_config.simulation_name = {"Value": "Debug_simulation", "Unit": "[-]"}

        self.Sim_config.generate_simulation_input(Path_to_input_dir = "Reference_simulations\\Debug_simulation",
                                                               Input_file_name = "Debug_simulation.XML")
            
        """ Run the simulation """
        filename = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\Utilities\\Reference_simulations\\Debug_simulation\\Debug_simulation.xml"
        args = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\x64\\Release\\Mjolnir_GRRT.exe -in " + filename + " -print_to_console 1"
        subprocess.call(args, shell = True)          

if __name__ == "__main__":

    Sim_instance = Simulation()
    Sim_instance.run_simulation()
    
    Sim_parser_instance = Simulation_Parser("C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\Utilities\\Reference_simulations\\Debug_simulation\\Debug_log")
    
    plt.plot(array(Sim_parser_instance.NT_Flux_r_coords) / Sim_instance.ADM_Mass, Sim_instance.ADM_Mass * Sim_instance.ADM_Mass * array(Sim_parser_instance.NT_Flux))
    plt.xlim([0, 3])
    plt.show()
