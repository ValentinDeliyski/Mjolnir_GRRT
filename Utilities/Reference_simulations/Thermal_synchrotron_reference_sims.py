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
from numpy import pi, tan, sqrt, arctan
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

class Thermal_syhnchrotron_reference_sims:
    
    def __init__(self):
        
        """ These parameters correspond to the ones in table 1 of https://arxiv.org/pdf/2206.12066. """
        
        self.Units = Units_class()
        self.Simulation_configurator = Simulation_configurator()
        
        """ Central black hole setup"""
        self.Simulation_configurator.metric_parameters.Metric_type = {"Value": "Kerr", "Unit": "[-]"}
        self.Simulation_configurator.object_mass                   = {"Value": 6.2e9,  "Unit": "[M_sun]"}
        
        self.Object_distance = {"Value": 16.9e6, "Unit": "[Pc]"}
        
        """ Accretion disk setup """
        self.Simulation_configurator.disk_model.Common_RIAF_disk_params.Temperature_scale_factor = {"Value": 1e11, "Unit": "[K]"}
        
        self.Simulation_configurator.disk_model.Common_RIAF_disk_params.Magnetization = {"Value": 0.01, "Unit": "[-]"}
        
        self.Simulation_configurator.disk_model.Common_RIAF_disk_params.Density_cutoff_radius     = {"Value": 0, "Unit": "[-]"}
        self.Simulation_configurator.disk_model.Common_RIAF_disk_params.Temperature_cutoff_radius = {"Value": 0, "Unit": "[-]"}
        
        self.Simulation_configurator.disk_model.Ensamble_type = {"Value": "Thermal",   "Unit": "[-]"}
        self.Simulation_configurator.disk_model.Disk_Model    = {"Value": "Phenom_RIAF_1", "Unit": "[-]"}
        
        self.Simulation_configurator.disk_model.Common_RIAF_disk_params.Density_power_law_power     = {"Value": 2.0, "Unit": "[-]"}
        self.Simulation_configurator.disk_model.Common_RIAF_disk_params.Temperature_power_law_power = {"Value": 1.0, "Unit": "[-]"}
        
        self.Simulation_configurator.disk_model.Velocity_profile = {"Value": "Theta Dependant", "Unit": "[-]"}
        
        """ Observer setup """
        self.Simulation_configurator.observer.Distance    = {"Value": 1e4,            "Unit": "[M]"}
        self.Simulation_configurator.observer.Inclination = {"Value": 160 * pi / 180, "Unit": "[Rad]"}
        self.Simulation_configurator.observer.Azimuth     = {"Value": 0,              "Unit": "[Rad]"}
        
        self.Simulation_configurator.observer.Obs_frequency = {"Value": 230e9, "Unit": "[Hz]"}
        self.Observer_FOV = {"Value": 100, "Unit": "[micro-arcsec]"}
        
        self.Simulation_configurator.observer.Image_y_min = {"Value": -(self.Object_distance["Value"] * self.Units.PC_TO_METER) / (self.Simulation_configurator.object_mass["Value"] * self.Units.M_SUN_SI * self.Units.GR_MASS_TO_METER) * tan(self.Observer_FOV["Value"] / 2 / self.Units.RAD_TO_MICRO_AS), "Unit": "[M]"}
        self.Simulation_configurator.observer.Image_y_max = {"Value":  (self.Object_distance["Value"] * self.Units.PC_TO_METER) / (self.Simulation_configurator.object_mass["Value"] * self.Units.M_SUN_SI * self.Units.GR_MASS_TO_METER) * tan(self.Observer_FOV["Value"] / 2 / self.Units.RAD_TO_MICRO_AS), "Unit": "[M]"}
        self.Simulation_configurator.observer.Image_x_min = {"Value": -(self.Object_distance["Value"] * self.Units.PC_TO_METER) / (self.Simulation_configurator.object_mass["Value"] * self.Units.M_SUN_SI * self.Units.GR_MASS_TO_METER) * tan(self.Observer_FOV["Value"] / 2 / self.Units.RAD_TO_MICRO_AS), "Unit": "[M]"}
        self.Simulation_configurator.observer.Image_x_max = {"Value":  (self.Object_distance["Value"] * self.Units.PC_TO_METER) / (self.Simulation_configurator.object_mass["Value"] * self.Units.M_SUN_SI * self.Units.GR_MASS_TO_METER) * tan(self.Observer_FOV["Value"] / 2 / self.Units.RAD_TO_MICRO_AS), "Unit": "[M]"}

        self.Simulation_configurator.observer.Resolution_x = {"Value": 256, "Unit": "[-]"}
        self.Simulation_configurator.observer.Resolution_y = {"Value": 256, "Unit": "[-]"}
        
        """ Kill the hotspot """
        self.Simulation_configurator.hotspot_model.Density_scale_factor = {"Value": 0, "Unit": "[g/cm^3]"}
    
        self.Simulation_configurator.geodesic_integrator.Integrator_type = {"Value": "RK78_Fehlberg", "Unit": "[-]"}
        self.Simulation_configurator.geodesic_integrator.RK_abs_accuracy = {"Value": 1e-12, "Unit": "[-]"}
        self.Simulation_configurator.geodesic_integrator.RK_rel_accuracy = {"Value": 1e-12, "Unit": "[-]"}
        
        self.Simulation_configurator.geodesic_integrator.max_upper_stepsize = {"Value": 100, "Unit": "[-]"}
        self.Simulation_configurator.geodesic_integrator.Max_rel_step_increase = {"Value": 10, "Unit": "[-]"}
        
        self.Simulation_configurator.emission_integrator.Rad_Transfer_Integrator_type = {"Value": "RK78_Fehlberg", "Unit": "[-]"}
        
        self.Simulation_configurator.observer.Include_polarization = {"Value": 0, "Unit": "[-]"}
        
    def Run_and_eval_sim_1(self):
        
        """ This simulation corresponds to the "low spin circ/thin" one in https://arxiv.org/pdf/2206.12066. Their results are presented on the top 4 panels
            of figure 3. They find a total flux of 546 mJy """
        
        """ Central black hole setup """
        self.Simulation_configurator.metric_parameters.Spin = {"Value": 0.01, "Unit": "[M]"}

        """ Accretion disk setup """
        self.Simulation_configurator.disk_model.Common_RIAF_disk_params.Opening_angle = {"Value": 0.1, "Unit": "[tan(angle)]"}
        self.Simulation_configurator.disk_model.Common_RIAF_disk_params.Density_power_law_scale     = {"Value": 1 + sqrt(1 - self.Simulation_configurator.metric_parameters.Spin["Value"]**2), "Unit": "[M]"}
        self.Simulation_configurator.disk_model.Common_RIAF_disk_params.Temperature_power_law_scale = {"Value": 1 + sqrt(1 - self.Simulation_configurator.metric_parameters.Spin["Value"]**2), "Unit": "[M]"}
        
        self.Simulation_configurator.disk_model.Common_RIAF_disk_params.Density_scale_factor = {"Value": 1.5e6, "Unit": "[g/cm^3]"}
        
        """ The simulation name and input file path """
        self.Simulation_configurator.simulation_name = {"Value": "Reference_Simulation_1", "Unit": "[-]"}
        
        """ The simulation output file path """
        self.Simulation_configurator.file_manager.Output_file_directory = parent_directory + "Reference_simulations"

        self.Simulation_configurator.generate_simulation_input(Path_to_input_dir = "Reference_simulations\\Reference_Simulation_1",
                                                               Input_file_name = "Reference_Simulation_1_input.XML")
        
        """ Run the simulation """
        filename = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\Utilities\\Reference_simulations\\Reference_Simulation_1\\Reference_Simulation_1_input.xml"
        args = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\x64\\Release\\Mjolnir_GRRT.exe -in " + filename + " -print_to_console 0"
        subprocess.call(args, shell = True)
                
        """ Evaluate the simulataion results """
        Sim_parser = Simulation_Parser(parent_directory + "Reference_simulations\\Reference_Simulation_1" + "\\Kerr")
        Total_flux = Sim_parser.get_total_flux(self.Units.M87_DISTANCE_GEOMETRICAL, unit = "mJy")
        
        print("Total Flux = {} [Jy]".format(round(Total_flux, 4)))
        
        Relative_error = abs(Total_flux - 546) / 546
        
        try:
            assert(Relative_error < 0.01)
            print(f"{bcolors.OKGREEN}Simulation 1 pass with a relative error of {{}} %.{bcolors.ENDC}".format(round(Relative_error * 100, 2)))
        except:
            print(f"{bcolors.FAIL}Simulation 1 fail with a relative error of {{}} %.{bcolors.ENDC}".format(round(Relative_error * 100, 2)))
        
    def Run_and_eval_sim_2(self):
        
        """ This simulation corresponds to the "low spin circ/thick" one in https://arxiv.org/pdf/2206.12066. Their results are presented on the third 4 panels
            of figure 3. They find a total flux of 651 mJy """
        
        """ Central black hole setup """
        self.Simulation_configurator.metric_parameters.Spin = {"Value": 0.01, "Unit": "[M]"}

        """ Accretion disk setup """
        self.Simulation_configurator.disk_model.Common_RIAF_disk_params.Opening_angle = {"Value": 1, "Unit": "[tan(angle)]"}
        self.Simulation_configurator.disk_model.Common_RIAF_disk_params.Density_power_law_scale  = {"Value": 1 + sqrt(1 - self.Simulation_configurator.metric_parameters.Spin["Value"]**2), "Unit": "[M]"}
        self.Simulation_configurator.disk_model.Common_RIAF_disk_params.Temperature_cutoff_scale = {"Value": 1 + sqrt(1 - self.Simulation_configurator.metric_parameters.Spin["Value"]**2), "Unit": "[M]"}
        
        self.Simulation_configurator.disk_model.Common_RIAF_disk_params.Density_scale_factor = {"Value": 0.7e6, "Unit": "[g/cm^3]"}
        
        """ The simulation name and input file path """
        self.Simulation_configurator.simulation_name = {"Value": "Reference_Simulation_2", "Unit": "[-]"}

        """ The simulation output file path """
        self.Simulation_configurator.file_manager.Output_file_directory = parent_directory + "Reference_simulations"

        self.Simulation_configurator.generate_simulation_input(Path_to_input_dir = "Reference_simulations\\Reference_Simulation_2",
                                                               Input_file_name = "Reference_Simulation_2_input.XML")
        
        """ Run the simulation """
        filename = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\Utilities\\Reference_simulations\\Reference_Simulation_2\\Reference_Simulation_2_input.xml"
        args = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\x64\\Release\\Mjolnir_GRRT.exe -in " + filename + " -print_to_console 0"
        subprocess.call(args, shell = True)
                
        """ Evaluate the simulataion results """
        Sim_parser = Simulation_Parser(parent_directory + "Reference_simulations\\Reference_Simulation_2" + "\\Kerr")
        Total_flux = Sim_parser.get_total_flux(self.Units.M87_DISTANCE_GEOMETRICAL, unit = "mJy")
        
        print("Total Flux = {} [Jy]".format(round(Total_flux, 4)))
        
        Relative_error = abs(Total_flux - 651) / 651
        
        try:
            assert(Relative_error < 0.01)
            print(f"{bcolors.OKGREEN}Simulation 2 pass with a relative error of {{}} %.{bcolors.ENDC}".format(round(Relative_error * 100, 2)))
        except:
            print(f"{bcolors.FAIL}Simulation 2 fail with a relative error of {{}} %.{bcolors.ENDC}".format(round(Relative_error * 100, 2)))
        
Thermal_syhnchrotron_reference_sims_instance = Thermal_syhnchrotron_reference_sims()

Sim_1_thread = threading.Thread(target = Thermal_syhnchrotron_reference_sims_instance.Run_and_eval_sim_1, args = [])
Sim_2_thread = threading.Thread(target = Thermal_syhnchrotron_reference_sims_instance.Run_and_eval_sim_2, args = [])

Sim_1_thread.start()
time.sleep(1)

Sim_2_thread.start()
