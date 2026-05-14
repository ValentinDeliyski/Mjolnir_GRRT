import sys
import os

""" Add the parent directory of this file to the search path, 
    so this script can be ran from the "Utilities" folder """
    
parent_directory = os.path.abspath('...')
sys.path.append(parent_directory)

from Mjolnir_Configurator import Simulation_configurator
from Support_functions.Parsers import Units_class, Simulation_Parser
from Support_functions.Spacetimes_new import Kerr

from numpy import pi, tan, sqrt, linspace
from numpy.typing import NDArray

import multiprocessing
import subprocess

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

class Simulation:
    
    def __init__(self, resolution, inclination):
        
        """ These parameters correspond to the ones in table 1 of https://arxiv.org/pdf/2309.10053. """
        
        self.Units = Units_class()
        self.Simulation_configurator = Simulation_configurator()
        
        self.Simulation_configurator.simulation_mode = {"Value": 0, "Unit": "[-]"} 
        self.Simulation_configurator.observer.Cam_rotation_angle = {"Value": 0, "Unit": "[-]"}
        
        """ Central black hole setup"""
        self.Simulation_configurator.metric_parameters.Metric_type = {"Value": "Kerr",  "Unit": "[-]"}
        self.Simulation_configurator.metric_parameters.Spin        = {"Value": 0.5, "Unit": "[M]"}
        
        """ Accretion disk setup """   
        self.Simulation_configurator.disk_model.Disk_Model = {"Value": "Numerical", "Unit": "[-]"}
        self.Simulation_configurator.disk_model.Max_density = {"Value": 1, "Unit": "[-]"}
        self.Simulation_configurator.disk_model.Ensamble_type = {"Value": "Phenomenological", "Unit": "[-]"}
        
        self.Simulation_configurator.emission_models.Absorbtion_coeff = {"Value": 1e2, "Unit": "[-]"}
        
        self.Simulation_configurator.disk_model.Numerical_disk_params.Numerical_XML_path = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\Utilities\\Test.XML"
        
        """ Kill the hotspot """
        self.Simulation_configurator.hotspot_model.Enabled_flag = {"Value": 0, "Unit": "[-]"}
        
        """ Observer setup """
        self.Simulation_configurator.observer.Distance    = {"Value": 1e4,         "Unit":  "[M]" }
        self.Simulation_configurator.observer.Inclination = {"Value": inclination * pi / 180, "Unit": "[Rad]"}
        self.Simulation_configurator.observer.Azimuth     = {"Value": 0,           "Unit": "[Rad]"}
        
        self.Simulation_configurator.observer.Cam_rotation_angle = {"Value": 0, "Unit": "[Hz]"}
        
        self.Observer_FOV = {"Value": 150, "Unit": "[micro-arcsec]"}
    
        self.Simulation_configurator.observer.Image_y_min = {"Value": -15, "Unit": "[M]"}
        self.Simulation_configurator.observer.Image_y_max = {"Value":  15, "Unit": "[M]"}
        self.Simulation_configurator.observer.Image_x_min = {"Value": -15, "Unit": "[M]"}
        self.Simulation_configurator.observer.Image_x_max = {"Value":  15, "Unit": "[M]"}
        
        self.Simulation_configurator.observer.Resolution_x = {"Value": resolution, "Unit": "[-]"}
        self.Simulation_configurator.observer.Resolution_y = {"Value": resolution, "Unit": "[-]"}
        
        """ Configure the integrator """
        
        self.Simulation_configurator.geodesic_integrator.Integrator_type = {"Value": "RK78_Fehlberg", "Unit": "[-]"}
        # self.Simulation_configurator.emission_integrator.Rad_Transfer_Integrator_type = {"Value": "Analytic", "Unit": "[-]"}        
        # # self.Simulation_configurator.geodesic_integrator.min_upper_stepsize = {"Value": 0.15, "Unit": "[-]"}
        # # self.Simulation_configurator.geodesic_integrator.max_step_b_coeff = {"Value": 0.004, "Unit": "[-]"}
        # # self.Simulation_configurator.geodesic_integrator.dist_at_min_upper_stepsize = {"Value": 25, "Unit": "[-]"}
        
        self.Simulation_configurator.geodesic_integrator.max_integration_count = {"Value": 1000000, "Unit": "[-]"}
        self.Simulation_configurator.geodesic_integrator.RK_abs_accuracy = {"Value": 1e-12, "Unit": "[-]"}
        self.Simulation_configurator.geodesic_integrator.RK_rel_accuracy = {"Value": 1e-12, "Unit": "[-]"}
        
        self.Simulation_configurator.geodesic_integrator.max_upper_stepsize = {"Value": 100, "Unit": "[-]"}
        self.Simulation_configurator.geodesic_integrator.Max_rel_step_increase = {"Value": 2, "Unit": "[-]"}

        self.Simulation_configurator.simulation_name = {"Value": "Numerical_disk_Kerr_a_0.5_inc_{}".format(inclination), "Unit": "[-]"}
        
        """ The simulation output file path """
        self.Simulation_configurator.file_manager.Output_file_directory = parent_directory + "Numerical_disks"

    def run_simulation(self):
        
        self.Simulation_configurator.generate_simulation_input(Path_to_input_dir = "Numerical_disks\\" + self.Simulation_configurator.simulation_name["Value"],
                                                               Input_file_name   = self.Simulation_configurator.simulation_name["Value"] + "_input.XML")
            
        """ Run the simulation """
        filename = ("C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\Utilities\\" + 
                    "Numerical_disks\\" + 
                    self.Simulation_configurator.simulation_name["Value"] + 
                    "\\" + 
                    self.Simulation_configurator.simulation_name["Value"] + 
                    "_input.XML")    
        
        args = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\x64\\Release\\Mjolnir_GRRT.exe -in " + filename + " -print_to_console 1"
        
        subprocess.call(args, shell = True)          

if __name__ == "__main__":

    Sim_instance = Simulation(inclination = 60, resolution = 256)
    Sim_instance.run_simulation()
