import sys
import os

""" Add the parent directory of this file to the search path, 
    so this script can be ran from the "Utilities" folder """
    
parent_directory = os.path.abspath('...')
sys.path.append(parent_directory)

from Mjolnir_Configurator import Simulation_configurator
from Support_functions.Parsers import Units_class, Simulation_Parser

from numpy import pi, tan, sqrt, linspace
from numpy.typing import NDArray

from multiprocessing import Pool
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
    
    def __init__(self):
        
        """ These parameters correspond to the ones in table 1 of https://arxiv.org/pdf/2309.10053. """
        
        self.Units = Units_class()
        self.Simulation_configurator = Simulation_configurator()
        self.Simulation_configurator.observer.Include_polarization = {"Value": 0, "Unit": "[-]"}
        self.Simulation_configurator.observer.Cam_rotation_angle = {"Value": 0, "Unit": "[-]"}
        
        """ Central black hole setup"""
        self.Simulation_configurator.metric_parameters.Metric_type = {"Value": "Kerr",  "Unit": "[-]"}
        self.Simulation_configurator.object_mass                   = {"Value": 4.297e6, "Unit": "[M_sun]"}
        self.Simulation_configurator.metric_parameters.Spin        = {"Value": 0.00001, "Unit": "[M]"}
        self.Object_distance                                       = {"Value": 8.277e3, "Unit": "[Pc]"}
        self.Simulation_configurator.max_image_order               = {"Value": 0, "Unit": "[-]"}
        
        """ Accretion disk setup """   
        self.Simulation_configurator.disk_model.Disk_Model = {"Value": "Novikov-Thorne", "Unit": "[-]"}
        
        self.Simulation_configurator.disk_model.r_in_NT_disk = {"Value": 6, "Unit": "[M]"} 
        self.Simulation_configurator.disk_model.r_out_NT_disk = {"Value": 25, "Unit": "[M]"} 
        
        self.Simulation_configurator.disk_model.Mag_field_geometry_r     = {"Value": 0, "Unit": "[-]"}
        self.Simulation_configurator.disk_model.Mag_field_geometry_theta = {"Value": 1, "Unit": "[-]"}
        self.Simulation_configurator.disk_model.Mag_field_geometry_phi   = {"Value": 0, "Unit": "[-]"}
        
        """ Kill the hotspot """
        self.Simulation_configurator.hotspot_model.Density_scale_factor = {"Value": 0, "Unit": "[g/cm^3]"}
        
        """ Observer setup """
        self.Simulation_configurator.observer.Distance    = {"Value": 1e4,            "Unit": "[M]"}
        self.Simulation_configurator.observer.Inclination = {"Value": 150 * pi / 180, "Unit": "[Rad]"}
        self.Simulation_configurator.observer.Azimuth     = {"Value": 0,              "Unit": "[Rad]"}
        
        self.Observer_FOV = {"Value": 400, "Unit": "[micro-arcsec]"}
        
        self.Simulation_configurator.observer.Image_y_min = {"Value": -(self.Object_distance["Value"] * self.Units.PC_TO_METER) / (self.Simulation_configurator.object_mass["Value"] * self.Units.M_SUN_SI * self.Units.GR_MASS_TO_METER) * tan(self.Observer_FOV["Value"] / 2 / self.Units.RAD_TO_MICRO_AS), "Unit": "[M]"}
        self.Simulation_configurator.observer.Image_y_max = {"Value":  (self.Object_distance["Value"] * self.Units.PC_TO_METER) / (self.Simulation_configurator.object_mass["Value"] * self.Units.M_SUN_SI * self.Units.GR_MASS_TO_METER) * tan(self.Observer_FOV["Value"] / 2 / self.Units.RAD_TO_MICRO_AS), "Unit": "[M]"}
        self.Simulation_configurator.observer.Image_x_min = {"Value": -(self.Object_distance["Value"] * self.Units.PC_TO_METER) / (self.Simulation_configurator.object_mass["Value"] * self.Units.M_SUN_SI * self.Units.GR_MASS_TO_METER) * tan(self.Observer_FOV["Value"] / 2 / self.Units.RAD_TO_MICRO_AS), "Unit": "[M]"}
        self.Simulation_configurator.observer.Image_x_max = {"Value":  (self.Object_distance["Value"] * self.Units.PC_TO_METER) / (self.Simulation_configurator.object_mass["Value"] * self.Units.M_SUN_SI * self.Units.GR_MASS_TO_METER) * tan(self.Observer_FOV["Value"] / 2 / self.Units.RAD_TO_MICRO_AS), "Unit": "[M]"}
        
        self.Simulation_configurator.observer.Resolution_x = {"Value": 256, "Unit": "[-]"}
        self.Simulation_configurator.observer.Resolution_y = {"Value": 256, "Unit": "[-]"}
        
        """ Configure the integrator """
        
        self.Simulation_configurator.integrator.RK78_abs_accuracy = {"Value": 1e-13, "Unit": "[-]"}
        self.Simulation_configurator.integrator.RK78_rel_accuracy = {"Value": 1e-13, "Unit": "[-]"}
        self.Simulation_configurator.integrator.Max_rel_step_increase  = {"Value": 5, "Unit": "[-]"}
        self.Simulation_configurator.integrator.max_stepsize = {"Value": 50, "Unit": "[-]"}
        
        """ The simulation output file path """
        self.Simulation_configurator.file_manager.Output_file_directory = parent_directory + "Reference_simulations"

    def run_simulation(self):
        
        self.Simulation_configurator.simulation_name = {"Value": "Thin_Disk_Reference_Simulation", "Unit": "[-]"}

        self.Simulation_configurator.generate_simulation_input(Path_to_input_dir = "Reference_simulations\\Thin_Disk_Reference_Simulation",
                                                               Input_file_name = "Thin_Disk_Reference_Simulation_input.XML")
            
        """ Run the simulation """
        filename = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\Utilities\\Reference_simulations\\Thin_Disk_Reference_Simulation\\Thin_Disk_Reference_Simulation_input.xml"
        args = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\x64\\Release\\Mjolnir_GRRT.exe -in " + filename + " -print_to_console 0"
        subprocess.call(args, shell = True)          

if __name__ == "__main__":

    Sim_instance = Simulation()
    Sim_instance.run_simulation()
