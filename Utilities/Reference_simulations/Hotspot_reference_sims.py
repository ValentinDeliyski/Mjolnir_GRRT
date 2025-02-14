import sys
import os
import threading
import time

""" Add the parent directory of this file to the search path, 
    so this script can be ran from the "Utilities" folder """
    
parent_directory = os.path.abspath('...')
sys.path.append(parent_directory)

from Mjolnir_Configurator import Simulation_configurator
from Support_functions.Parsers import Units_class, Simulation_Parser
from numpy import pi, tan, sqrt
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

class Hotspot_reference_sims:
    
    def __init__(self):
        
        """ These parameters correspond to the ones in table 2 of https://arxiv.org/pdf/2301.11874. """
        
        self.Units = Units_class()
        self.Simulation_configurator = Simulation_configurator()
        self.Simulation_configurator.average_emission_pitch_angle = {"Value": 1, "Unit": "[-]"}
        
        """ Central black hole setup"""
        self.Simulation_configurator.metric_parameters.Metric_type = {"Value": "Kerr",  "Unit": "[-]"}
        self.Simulation_configurator.object_mass                   = {"Value": 4.297e6, "Unit": "[M_sun]"}
        self.Simulation_configurator.metric_parameters.Spin        = {"Value": 0.000, "Unit": "[M]"}
        self.Object_distance                                       = {"Value": 8.277e3, "Unit": "[Pc]"}
        
        """ Kill the accretion disk setup """   
        self.Simulation_configurator.disk_model.Density_scale_factor = {"Value": 0, "Unit": "[g/cm^3]"}

        """ Hotspot setup """
        self.Simulation_configurator.hotspot_model.Density_scale_factor     = {"Value": 1.05e7, "Unit": "[g/cm^3]"}
        self.Simulation_configurator.hotspot_model.Temperature_scale_factor = {"Value": 3.01 * 3e10,  "Unit": "[K]"}
        
        self.Simulation_configurator.hotspot_model.Density_profile     = {"Value": "Sphere", "Unit": "[-]"}
        self.Simulation_configurator.hotspot_model.Temperature_profile = {"Value": "Sphere", "Unit": "[-]"}
        
        self.Simulation_configurator.hotspot_model.Temporal_spread   = {"Value": 850000,      "Unit": "[GM/c^3]"}
        self.Simulation_configurator.hotspot_model.Magnetization     = {"Value": 0.01,        "Unit": "[-]"}
        self.Simulation_configurator.emission_models.Kappa           = {"Value": 5,           "Unit": "[-]"}
        self.Simulation_configurator.hotspot_model.Ensamble_type     = {"Value": "Kappa",     "Unit": "[-]"}
        self.Simulation_configurator.hotspot_model.Distance          = {"Value": 9,           "Unit": "[M]"} 
        self.Simulation_configurator.hotspot_model.Azimuth           = {"Value": -pi/2,       "Unit": "[M]"} 
        self.Simulation_configurator.hotspot_model.Velocity_profile  = {"Value": "Keplarian", "Unit": "[-]"}
        self.Simulation_configurator.hotspot_model.Mag_field_geometry_X = {"Value": 0, "Unit": "[-]"}
        self.Simulation_configurator.hotspot_model.Mag_field_geometry_Y = {"Value": 0, "Unit": "[-]"}
        self.Simulation_configurator.hotspot_model.Mag_field_geometry_Z = {"Value": 1, "Unit": "[-]"}
        
        """ Observer setup """
        self.Simulation_configurator.observer.Distance    = {"Value": 1e4,           "Unit": "[M]"}
        self.Simulation_configurator.observer.Inclination = {"Value": 20 * pi / 180, "Unit": "[Rad]"}
        self.Simulation_configurator.observer.Azimuth     = {"Value": 0,             "Unit": "[Rad]"}
        
        self.Simulation_configurator.observer.Obs_frequency = {"Value": Units_class().C_LIGHT_SI / 2.2e-6, "Unit": "[Hz]"}
        self.Observer_FOV = {"Value": 200, "Unit": "[micro-arcsec]"}
        
        self.Simulation_configurator.observer.Image_y_min = {"Value": -(self.Object_distance["Value"] * self.Units.PC_TO_METER) / (self.Simulation_configurator.object_mass["Value"] * self.Units.M_SUN_SI * self.Units.GR_MASS_TO_METER) * tan(self.Observer_FOV["Value"] / 2 / self.Units.RAD_TO_MICRO_AS), "Unit": "[M]"}
        self.Simulation_configurator.observer.Image_y_max = {"Value":  (self.Object_distance["Value"] * self.Units.PC_TO_METER) / (self.Simulation_configurator.object_mass["Value"] * self.Units.M_SUN_SI * self.Units.GR_MASS_TO_METER) * tan(self.Observer_FOV["Value"] / 2 / self.Units.RAD_TO_MICRO_AS), "Unit": "[M]"}
        self.Simulation_configurator.observer.Image_x_min = {"Value": -(self.Object_distance["Value"] * self.Units.PC_TO_METER) / (self.Simulation_configurator.object_mass["Value"] * self.Units.M_SUN_SI * self.Units.GR_MASS_TO_METER) * tan(self.Observer_FOV["Value"] / 2 / self.Units.RAD_TO_MICRO_AS), "Unit": "[M]"}
        self.Simulation_configurator.observer.Image_x_max = {"Value":  (self.Object_distance["Value"] * self.Units.PC_TO_METER) / (self.Simulation_configurator.object_mass["Value"] * self.Units.M_SUN_SI * self.Units.GR_MASS_TO_METER) * tan(self.Observer_FOV["Value"] / 2 / self.Units.RAD_TO_MICRO_AS), "Unit": "[M]"}
        
        self.Simulation_configurator.observer.Resolution_x = {"Value": 512, "Unit": "[-]"}
        self.Simulation_configurator.observer.Resolution_y = {"Value": 512, "Unit": "[-]"}
       
        """ Kill the Novikov-Thorne disk """
        self.Simulation_configurator.NT_model_params.Evaluate_NT_disk = {"Value": 0, "Unit": "[-]"}
        
        """ Configure the integrator """
        self.Simulation_configurator.integrator.Step_controller_type = {"Value": "Gustafsson", "Unit": "[-]"}
        self.Simulation_configurator.integrator.RK45_accuracy        = {"Value": 1e-12, "Unit": "[-]"}
        
        """ The simulation name and input file path """
        
        self.Simulation_configurator.file_manager.Output_file_directory = parent_directory + "Reference_simulations"

        hotspot_azimuth_position_number = 10

        for time_offset in range(0, hotspot_azimuth_position_number):
            
            self.Simulation_configurator.hotspot_model.Coord_time_at_max = {"Value": 2 * pi / (1 / sqrt(9**3)) * time_offset / hotspot_azimuth_position_number, "Unit": "[GM/c^3]"}
            
            self.Simulation_configurator.simulation_name = {"Value": "Hotspot_Reference_Simulation_{}".format(time_offset), "Unit": "[-]"}

            

            self.Simulation_configurator.generate_simulation_input(Path_to_input_dir = "Reference_simulations\\Hotspot_Reference_Simulation",
                                                                   Input_file_name = "Hotspot_Reference_Simulation_input.XML")
            
            """ Run the simulation """
            filename = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\Utilities\\Reference_simulations\\Hotspot_Reference_Simulation\\Hotspot_Reference_Simulation_input.xml"
            args = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\x64\\Release\\Mjolnir_GRRT.exe -in " + filename + " -print_to_console 1"
            subprocess.call(args, shell = True)
                    
            """ Evaluate the simulataion results """
            Sim_parser_n0 = Simulation_Parser(parent_directory + "Reference_simulations\\Hotspot_Reference_Simulation_{}".format(time_offset) + "\\Kerr_n0")
            Total_flux_n0 = Sim_parser_n0.get_total_flux(self.Units.SGRA_DISTANCE_GEOMETRICAL, unit = "mJy")
            
            Sim_parser_n1 = Simulation_Parser(parent_directory + "Reference_simulations\\Hotspot_Reference_Simulation_{}".format(time_offset) + "\\Kerr_n1")
            Total_flux_n1 = Sim_parser_n1.get_total_flux(self.Units.SGRA_DISTANCE_GEOMETRICAL, unit = "mJy")
            
            Sim_parser_n2 = Simulation_Parser(parent_directory + "Reference_simulations\\Hotspot_Reference_Simulation_{}".format(time_offset) + "\\Kerr_n2")
            Total_flux_n2 = Sim_parser_n2.get_total_flux(self.Units.SGRA_DISTANCE_GEOMETRICAL, unit = "mJy")
            
            Sim_parser_n3 = Simulation_Parser(parent_directory + "Reference_simulations\\Hotspot_Reference_Simulation_{}".format(time_offset) + "\\Kerr_n3")
            Total_flux_n3 = Sim_parser_n3.get_total_flux(self.Units.SGRA_DISTANCE_GEOMETRICAL, unit = "mJy")
            
            Total_flux = Total_flux_n0 + Total_flux_n1 + Total_flux_n2 + Total_flux_n3
            
            print(Total_flux / 15)

Hotspot_reference_sims_instance = Hotspot_reference_sims()