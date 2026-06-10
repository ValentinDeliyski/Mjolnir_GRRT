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
    
    def __init__(self, resolution, spin, inclination):
 
        self.Units = Units_class()
        self.Simulation_configurator = Simulation_configurator()
        
        self.Simulation_configurator.simulation_mode = {"Value": 0, "Unit": "[-]"} 
        self.Simulation_configurator.observer.Cam_rotation_angle = {"Value": 0, "Unit": "[-]"}
        
        """ Central black hole setup"""
        self.Simulation_configurator.metric_parameters.Metric_type = {"Value": "Kerr",  "Unit": "[-]"}
        self.Simulation_configurator.metric_parameters.Mass        = {"Value": 1, "Unit": "[M]"}
        self.Simulation_configurator.metric_parameters.Spin        = {"Value": spin, "Unit": "[M]"}
        
        """ Accretion disk setup """   
        self.Simulation_configurator.disk_model.Disk_Model = {"Value": "Novikov-Thorne", "Unit": "[-]"}
        
        Kerr_instance = Kerr(self.Simulation_configurator.metric_parameters.Mass["Value"], self.Simulation_configurator.metric_parameters.Spin["Value"])
        
        self.Simulation_configurator.disk_model.NT_disk_params.r_in = {"Value": Kerr_instance.get_ISCO()[0], "Unit": "[M]"} 
        self.Simulation_configurator.disk_model.NT_disk_params.r_out = {"Value": 25 * self.Simulation_configurator.metric_parameters.Mass["Value"], "Unit": "[M]"} 
         
        """ Kill the hotspot """
        self.Simulation_configurator.hotspot_model.Enabled_flag = {"Value": 0, "Unit": "[g/cm^3]"}
        
        """ Observer setup """
        self.Simulation_configurator.observer.Distance    = {"Value": 1e4,         "Unit":  "[M]" }
        self.Simulation_configurator.observer.Inclination = {"Value": inclination * pi / 180, "Unit": "[Rad]"}
        self.Simulation_configurator.observer.Azimuth     = {"Value": 0,           "Unit": "[Rad]"}
        
        self.Simulation_configurator.observer.Cam_rotation_angle = {"Value": 0, "Unit": "[Hz]"}
        
        self.Observer_FOV = {"Value": 150, "Unit": "[micro-arcsec]"}
    
        self.Simulation_configurator.observer.Image_y_min = {"Value": -15, "Unit": "[M]"}
        self.Simulation_configurator.observer.Image_y_max = {"Value":  15, "Unit": "[M]"}
        self.Simulation_configurator.observer.Image_x_min = {"Value": -35, "Unit": "[M]"}
        self.Simulation_configurator.observer.Image_x_max = {"Value":  35, "Unit": "[M]"}
        
        self.Simulation_configurator.observer.Resolution_x = {"Value": resolution, "Unit": "[-]"}
        self.Simulation_configurator.observer.Resolution_y = {"Value": resolution, "Unit": "[-]"}
        
        """ Configure the integrator """
        
        self.Simulation_configurator.geodesic_integrator.Integrator_type = {"Value": "RK78_Fehlberg", "Unit": "[-]"}
        self.Simulation_configurator.geodesic_integrator.max_integration_count = {"Value": 1000000, "Unit": "[-]"}
        self.Simulation_configurator.geodesic_integrator.RK_abs_accuracy = {"Value": 1e-12, "Unit": "[-]"}
        self.Simulation_configurator.geodesic_integrator.RK_rel_accuracy = {"Value": 1e-12, "Unit": "[-]"}
        
        self.Simulation_configurator.geodesic_integrator.max_upper_stepsize = {"Value": 100, "Unit": "[-]"}
        self.Simulation_configurator.geodesic_integrator.Max_rel_step_increase = {"Value": 2, "Unit": "[-]"}

        self.Simulation_configurator.simulation_name = {"Value": "Kerr_a_{}_inc_{}".format(spin, inclination), "Unit": "[-]"}
        
        """ The simulation output file path """
        self.Simulation_configurator.file_manager.Output_file_directory = parent_directory + "Wormhole_sim_paper"

    def run_simulation(self):
        
        self.Simulation_configurator.generate_simulation_input(Path_to_input_dir = "Wormhole_sim_paper\\" + self.Simulation_configurator.simulation_name["Value"],
                                                               Input_file_name   = self.Simulation_configurator.simulation_name["Value"] + "_input.XML")
            
        """ Run the simulation """
        filename = ("C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\Utilities\\" + 
                    "Wormhole_sim_paper\\" + 
                    self.Simulation_configurator.simulation_name["Value"] + 
                    "\\" + 
                    self.Simulation_configurator.simulation_name["Value"] + 
                    "_input.XML")    
        
        args = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\x64\\Release\\Mjolnir_GRRT.exe -in " + filename + " -print_to_console 1"
        
        subprocess.call(args, shell = True)          

if __name__ == "__main__":

    Processses = []

    for spin in [0]:
        
        for inc in [80]:
            
            Sim_instance = Simulation(spin = spin,
                                      inclination = inc,
                                      resolution = 2048)
            
            Process = multiprocessing.Process(target = Sim_instance.run_simulation)  
                
            Process.start()
            Processses.append(Process)     

    for Process in Processses:
        Process.join()
