import sys
import os

""" Add the parent directory of this file to the search path, 
    so this script can be ran from the "Utilities" folder """
    
parent_directory = os.path.abspath('...')
sys.path.append(parent_directory)

from Mjolnir_Configurator import Simulation_configurator
from Support_functions.Parsers import Units_class, Simulation_Parser
from Support_functions.Spacetimes_new import Kerr, Gauss_Bonnet

from numpy import pi, tan, sqrt, linspace, arctan, array
from numpy.typing import NDArray

import subprocess

import matplotlib.pyplot as plt

class Simulation:
    
    def __init__(self):
        
        self.Sim_config = Simulation_configurator()
        Kerr_instance = Kerr(mass = 1, spin_param = 0.98)
        Gauss_Bonnet_instance = Gauss_Bonnet(Param = 0.48)
        
        self.Sim_config.simulation_mode = {"Value": 3, "Unit": "[-]"} 
        
        self.Sim_config.metric_parameters.Numerical_metric_spline_path = "C:/Users/Valur/Documents/Repos/Mjolnir_GRRT/Utilities/Galin_zero_curvature_config_VI.XML"

        self.Sim_config.min_image_order = {"Value": 0, "Unit": "[-]"}
        self.Sim_config.max_image_order = {"Value": 10, "Unit": "[-]"}

        """ ================================================== The Numerical Metric ================================================== """
        
        self.Sim_config.metric_parameters.Metric_type = {"Value": "Einstein-Gauss-Bonnet", "Unit": "[-]"}
        self.Sim_config.metric_parameters.Spin        = {"Value": 0.98, "Unit": "[M]"}
        
        """ ================================================ Observer ================================================ """
   
        self.Sim_config.observer.Distance    = {"Value": 1e4 , "Unit": "[M]"}
        self.Sim_config.observer.Inclination = {"Value": 80 * pi / 180, "Unit": "[Rad]"}
        self.Sim_config.observer.Obs_frequency = {"Value": 230e9, "Unit": "[Hz]"}
        
        self.Sim_config.observer.Cam_rotation_angle = {"Value": 0, "Unit": "[Hz]"}
        self.Sim_config.observer.Use_angular_coords = {"Value": 1, "Unit": "[M]"}
        
        """ ================================================== Disk ================================================== """
        
        self.Sim_config.disk_model.Disk_Model = {"Value": "Novikov-Thorne", "Unit": "[-]"}
        
        self.Sim_config.disk_model.r_in_NT_disk = {"Value": Gauss_Bonnet_instance.get_ISCO()[0], "Unit": "[M]"} 
        self.Sim_config.disk_model.r_out_NT_disk = {"Value": 25, "Unit": "[M]"} 
        
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
    
    plt.title(r"Kerr $F(r)$, $a = 0.98M$")
    plt.plot(array(Sim_parser_instance.NT_Flux_r_coords), array(Sim_parser_instance.NT_Flux))
    plt.xlim([0, 10])
    plt.show()
