import sys
import os
import threading
import time
from multiprocessing import Pool

""" Add the parent directory of this file to the search path, 
    so this script can be ran from the "Utilities" folder """
    
parent_directory = os.path.abspath('...')
sys.path.append(parent_directory)

from Mjolnir_Configurator import Simulation_configurator
from Support_functions.Parsers import Units_class, Simulation_Parser

from numpy import pi, tan, sqrt, linspace
from numpy.typing import NDArray

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

class Thermal_hotspot_reference_sims:
    
    def __init__(self):
        
        """ These parameters correspond to the ones in table 1 of https://arxiv.org/pdf/2209.09931. """
        
        self.Units = Units_class()
        self.Simulation_configurator = Simulation_configurator()
        
        """ Central black hole setup"""
        self.Simulation_configurator.metric_parameters.Metric_type = {"Value": "Einstein-Gauss-Bonnet", "Unit": "[-]"}
        self.Simulation_configurator.metric_parameters.Spin        = {"Value": 0, "Unit": "[M]"}
        self.Simulation_configurator.object_mass                   = {"Value": 4.297e6,  "Unit": "[M_sun]"}
        
        self.Object_distance = {"Value": 8.277e3, "Unit": "[Pc]"}
        
        """ Accretion disk setup """
        
        self.Simulation_configurator.disk_model.Ensamble_type = {"Value": "Thermal",   "Unit": "[-]"}
        self.Simulation_configurator.disk_model.Disk_Model    = {"Value": "Phenom_RIAF_2", "Unit": "[-]"}
        
        self.Simulation_configurator.average_emission_pitch_angle    = {"Value": 0, "Unit": "[-]"}
        
        self.Simulation_configurator.disk_model.Density_scale_factor = {"Value": 6e6, "Unit": "[g/cm^3]"}
        self.Simulation_configurator.disk_model.Opening_angle        = {"Value": 0.3, "Unit": "[cos(angle)]"}
        
        Temp_coversion_factor = Units_class.BOLTZMANN_SI / Units_class.M_ELECTRON_SI / Units_class.C_LIGHT_SI**2
        self.Simulation_configurator.disk_model.Temperature_scale_factor = {"Value": 200 / Temp_coversion_factor, "Unit": "[K]"}
        
        self.Simulation_configurator.disk_model.Density_cutoff_radius     = {"Value": 6, "Unit": "[-]"}
        self.Simulation_configurator.disk_model.Temperature_cutoff_radius = {"Value": 6, "Unit": "[-]"}
        
        self.Simulation_configurator.disk_model.Density_cutoff_scale     = {"Value": 0.001, "Unit": "[-]"}
        self.Simulation_configurator.disk_model.Temperature_cutoff_scale = {"Value": 0.001, "Unit": "[-]"}
        
        self.Simulation_configurator.disk_model.Density_power_law_power     = {"Value": 3 / 2, "Unit": "[-]"}
        self.Simulation_configurator.disk_model.Density_power_law_scale     = {"Value": 1, "Unit": "[-]"}
        self.Simulation_configurator.disk_model.Temperature_power_law_power = {"Value": 0.84, "Unit": "[-]"}
        self.Simulation_configurator.disk_model.Temperature_power_law_scale = {"Value": 1, "Unit": "[-]"}
        
        self.Simulation_configurator.disk_model.Mag_field_magnitude_profile = {"Value": "Power_law_based", "Unit": "[-]"}
        self.Simulation_configurator.disk_model.Mag_field_magnitude_scale   = {"Value": 100, "Unit": "[G]"}
        self.Simulation_configurator.disk_model.Mag_field_power             = {"Value": 1, "Unit": "[-]"}
        self.Simulation_configurator.disk_model.Mag_field_radial_scale      = {"Value": 1, "Unit": "[-]"}
        
        self.Simulation_configurator.disk_model.Mag_field_geometry       = {"Value": "Constant", "Unit": "[-]"}
        self.Simulation_configurator.disk_model.Mag_field_geometry_r     = {"Value": 0, "Unit": "[-]"}
        self.Simulation_configurator.disk_model.Mag_field_geometry_theta = {"Value": 0, "Unit": "[-]"}
        self.Simulation_configurator.disk_model.Mag_field_geometry_phi   = {"Value": 1, "Unit": "[-]"}
        
        self.Simulation_configurator.disk_model.Velocity_profile = {"Value": "Theta Dependant", "Unit": "[-]"}
        
        """ Observer setup """
        self.Simulation_configurator.observer.Distance    = {"Value": 1e4,            "Unit": "[M]"}
        self.Simulation_configurator.observer.Inclination = {"Value": 20 * pi / 180, "Unit": "[Rad]"}
        self.Simulation_configurator.observer.Azimuth     = {"Value": 0,              "Unit": "[Rad]"}
        
        self.Simulation_configurator.observer.Obs_frequency = {"Value": 230e9, "Unit": "[Hz]"}
        
        self.Simulation_configurator.observer.Image_y_min = {"Value": -20, "Unit": "[M]"}
        self.Simulation_configurator.observer.Image_y_max = {"Value":  20, "Unit": "[M]"}
        self.Simulation_configurator.observer.Image_x_min = {"Value": -20, "Unit": "[M]"}
        self.Simulation_configurator.observer.Image_x_max = {"Value":  20, "Unit": "[M]"}
        
        self.Simulation_configurator.observer.Resolution_x = {"Value": 256, "Unit": "[-]"}
        self.Simulation_configurator.observer.Resolution_y = {"Value": 256, "Unit": "[-]"}
    
        self.Simulation_configurator.integrator.RK78_accuracy         = {"Value": 1e-13, "Unit": "[-]"}
        self.Simulation_configurator.integrator.Max_rel_step_increase = {"Value": 5, "Unit": "[-]"}
        self.Simulation_configurator.observer.Include_polarization    = {"Value": 1, "Unit": "[-]"}

    def Run_background_sim(self):
              
        """ Kill the hotspot """
        self.Simulation_configurator.hotspot_model.Density_scale_factor = {"Value": 0, "Unit": "[g/cm^3]"}
        
        """ The simulation name and input file path """
        self.Simulation_configurator.simulation_name = {"Value": "Thermal_Hotspot", "Unit": "[-]"}
        
        """ The simulation output file path """
        self.Simulation_configurator.file_manager.Output_file_directory = parent_directory + "Reference_simulations"

        self.Simulation_configurator.generate_simulation_input(Path_to_input_dir = "Reference_simulations\\Thermal_Hotspot",
                                                               Input_file_name = "Background_only_input.XML")
        
        """ Run the simulation """
        filename = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\Utilities\\Reference_simulations\\Thermal_Hotspot\\Background_only_input.xml"
        args = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\x64\\Release\\Mjolnir_GRRT.exe -in " + filename + " -print_to_console 0"
        subprocess.call(args, shell = True)

    def Run_sim_w_hotspot(self, obs_time: float, idx: int):
              
        self.Simulation_configurator.observer.Init_time = {"Value": obs_time, "Unit": "[GM/c^3]"}
        self.Simulation_configurator.thermalize_emission_medium = {"Value": 1, "Unit": "[GM/c^3]"}
            
        """ Setup the hotspot """
        self.Simulation_configurator.hotspot_model.Ensamble_type        = {"Value": "Thermal", "Unit": "[-]"}
        
        self.Simulation_configurator.hotspot_model.Density_scale_factor = {"Value": 2e6, "Unit": "[g/cm^3]"}
        self.Simulation_configurator.hotspot_model.Density_profile      = {"Value": "Gaussian", "Unit": "[-]"}
        self.Simulation_configurator.hotspot_model.Density_spread       = {"Value": 1.5, "Unit": "[M]"}
        
        Temp_coversion_factor = Units_class.BOLTZMANN_SI / Units_class.M_ELECTRON_SI / Units_class.C_LIGHT_SI**2
        self.Simulation_configurator.hotspot_model.Temperature_scale_factor = {"Value": 200 / Temp_coversion_factor, "Unit": "[g/cm^3]"}
        
        self.Simulation_configurator.hotspot_model.Temperature_profile      = {"Value": "Hybrid_power_law_gaussian", "Unit": "[-]"}
        self.Simulation_configurator.hotspot_model.Temperature_spread       = {"Value": 1.5, "Unit": "[M]"}
        self.Simulation_configurator.hotspot_model.Temperature_power_law_power = {"Value": 0.84, "Unit": "[M]"}
        self.Simulation_configurator.hotspot_model.Density_power_law_scale     = {"Value": 1.0, "Unit": "[M]"}
        
        self.Simulation_configurator.hotspot_model.Temporal_spread   = {"Value": 850000,      "Unit": "[GM/c^3]"}
        self.Simulation_configurator.hotspot_model.Ensamble_type     = {"Value": "Thermal",   "Unit": "[-]"}
        self.Simulation_configurator.hotspot_model.Distance          = {"Value": 11,          "Unit": "[M]"} 
        self.Simulation_configurator.hotspot_model.Azimuth           = {"Value": 0,          "Unit": "[M]"} 
        self.Simulation_configurator.hotspot_model.Velocity_profile  = {"Value": "Keplarian", "Unit": "[-]"}
        
        self.Simulation_configurator.hotspot_model.Mag_field_magnitude_profile = {"Value": "Power_law_based", "Unit": "[-]"}
        self.Simulation_configurator.hotspot_model.Mag_field_magnitude_scale   = {"Value": 100, "Unit": "[G]"}
        self.Simulation_configurator.hotspot_model.Mag_field_power             = {"Value": 1, "Unit": "[-]"}
        self.Simulation_configurator.hotspot_model.Mag_field_radial_scale      = {"Value": 1, "Unit": "[-]"}
        
        self.Simulation_configurator.hotspot_model.Mag_field_geometry       = {"Value": "Constant", "Unit": "[-]"}
        self.Simulation_configurator.hotspot_model.Mag_field_geometry_r     = {"Value": 0, "Unit": "[-]"}
        self.Simulation_configurator.hotspot_model.Mag_field_geometry_theta = {"Value": 0, "Unit": "[-]"}
        self.Simulation_configurator.hotspot_model.Mag_field_geometry_phi   = {"Value": 1, "Unit": "[-]"}
        
        """ The simulation name and input file path """
        self.Simulation_configurator.simulation_name = {"Value": "Thermal_Hotspot_{}".format(idx), "Unit": "[-]"}
        
        """ The simulation output file path """
        self.Simulation_configurator.file_manager.Output_file_directory = parent_directory + "Reference_simulations"

        self.Simulation_configurator.generate_simulation_input(Path_to_input_dir = "Reference_simulations\\Thermal_Hotspot_{}".format(idx),
                                                               Input_file_name = "Background_only_input.XML")
        
        """ Run the simulation """
        filename = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\Utilities\\Reference_simulations\\Thermal_Hotspot_{}\\Background_only_input.xml".format(idx)
        args = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\x64\\Release\\Mjolnir_GRRT.exe -in " + filename + " -print_to_console 0"
        subprocess.call(args, shell = True)

if __name__ == "__main__":

    Thermal_hotspot_reference_sims_instance = Thermal_hotspot_reference_sims()

    Spot_period: float = 2 * pi * Thermal_hotspot_reference_sims_instance.Simulation_configurator.hotspot_model.Distance["Value"]**(3 / 2)
     
    Hotspot_number: int = 40
    
    Obs_times: NDArray = linspace(0, Spot_period, Hotspot_number)
    File_idx: NDArray  = linspace(0, Hotspot_number - 1, Hotspot_number)

    Sim_args: list[list[float | int]] = []

    for t_obs, idx in zip(Obs_times, File_idx):
        Sim_args.append([t_obs, int(idx)])
        
    # with Pool(10) as pool:
    #     pool.starmap(Thermal_hotspot_reference_sims_instance.Run_sim_w_hotspot, Sim_args)
        
    Thermal_hotspot_reference_sims_instance.Run_sim_w_hotspot(25, int(File_idx[0]))
