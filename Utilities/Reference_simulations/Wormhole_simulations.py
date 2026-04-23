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
from Support_functions.Spacetimes_new import Kerr
from numpy import pi, tan, sqrt, array, where, all
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

class Wormhole_simulation_cofigurator:
    
    def __init__(self, spin = 0.0, redshift_f = 0.0, inclination = pi / 4, central_object = "M87", obs_freq = 230e9, resolution = 1024):
        
        """ These parameters correspond to the ones in table 1 of https://arxiv.org/pdf/2206.12066. """
        
        self.Units = Units_class()
        self.Kerr_instance = Kerr(mass = 1, spin_param = spin)
        
        self.Simulation_configurator = Simulation_configurator()
        self.central_object = central_object
        
        """ Central Wormhole setup"""
        
        self.Simulation_configurator.metric_parameters.Metric_type = {"Value": "Wormhole", "Unit": "[-]"}
        self.Simulation_configurator.metric_parameters.WH_redshift = {"Value": redshift_f, "Unit": "[M]"}
        self.Simulation_configurator.metric_parameters.Spin        = {"Value":  spin,    "Unit": "[M]"}
        self.Simulation_configurator.metric_parameters.WH_r_throat = {"Value":    1,    "Unit": "[M]"}
        
        self.Simulation_configurator.metric_parameters.WH_stop_at_throat = {"Value":    0,    "Unit": "[M]"}
        
        match self.central_object:
            case "M87":
                self.Simulation_configurator.object_mass = {"Value": self.Units.M_M87_BH_SI / self.Units.M_SUN_SI,  "Unit": "[M_sun]"}
                self.Object_distance = {"Value": self.Units.M87_DISTANCE_PC, "Unit": "[Pc]"}
                
                r_0 = self.Simulation_configurator.metric_parameters.WH_r_throat["Value"]
                T_max = 0.63e11
                n_max = 5e5
                
        
            case "Sgr_A":
                self.Simulation_configurator.object_mass = {"Value": self.Units.M_SGRA_BH_SI / self.Units.M_SUN_SI,  "Unit": "[M_sun]"}
                self.Object_distance = {"Value": self.Units.SGRA_DISTANCE_PC, "Unit": "[Pc]"}
                
                r_0 = self.Simulation_configurator.metric_parameters.WH_r_throat["Value"]
                T_max = 1e11
                n_max = 5e6
                
                self.Simulation_configurator.disk_model.Density_cutoff_scale = {"Value": 0.4 * 9.75684 / 13.8033 , "Unit": "[-]"}
                
            case _:
                
                print("Wrong central object!")
                
                exit(1)
        
        """ =============== Accretion disk setup =============== """
        
        self.Simulation_configurator.disk_model.Ensamble_type = {"Value": "Thermal",   "Unit": "[-]"}
        self.Simulation_configurator.disk_model.Disk_Model    = {"Value": "Phenom_RIAF_3", "Unit": "[-]"}
        
        """ ----------- Density ----------- """
        
        self.Simulation_configurator.disk_model.Density_scale_factor   = {"Value": n_max, "Unit":   "[g/cm^3]"  }
        self.Simulation_configurator.disk_model.Density_cutoff_radius  = {"Value":  r_0,  "Unit":     "[-]"     }
        self.Simulation_configurator.disk_model.Density_cutoff_scale   = {"Value":   5,   "Unit":     "[M]"     }
        self.Simulation_configurator.disk_model.Opening_angle          = {"Value":  0.1,  "Unit": "[tan(angle)]"}
        
        """ --------- Temperature --------- """
        
        self.Simulation_configurator.disk_model.Temperature_scale_factor   = {"Value": T_max, "Unit": "[K]"}
        self.Simulation_configurator.disk_model.Temperature_cutoff_radius  = {"Value":  r_0,  "Unit": "[M]"}
        self.Simulation_configurator.disk_model.Temperature_cutoff_scale   = {"Value":  5,    "Unit": "[M]"}
        
        """ --------- Misc --------- """
        
        self.Simulation_configurator.disk_model.Magnetization    = {"Value":      0.01,         "Unit": "[-]"}
        self.Simulation_configurator.disk_model.Ensamble_type    = {"Value":    "Thermal",      "Unit": "[-]"}
        self.Simulation_configurator.disk_model.Velocity_profile = {"Value": "Theta Dependant", "Unit": "[-]"}

        """ =============== Observer setup =============== """
        
        self.Simulation_configurator.observer.Distance    = {"Value": 1e4,         "Unit":  "[M]" }
        self.Simulation_configurator.observer.Inclination = {"Value": inclination, "Unit": "[Rad]"}
        self.Simulation_configurator.observer.Azimuth     = {"Value": 0,           "Unit": "[Rad]"}
        
        self.Simulation_configurator.observer.Obs_frequency = {"Value": obs_freq, "Unit": "[Hz]"}
        self.Simulation_configurator.observer.Cam_rotation_angle = {"Value": -70 * pi / 180 - pi / 4, "Unit": "[Hz]"}
       
        self.Observer_FOV = {"Value": 100, "Unit": "[micro-arcsec]"}
    
        self.Simulation_configurator.observer.Image_y_min = {"Value": -(self.Object_distance["Value"] * self.Units.PC_TO_METER) / (self.Simulation_configurator.object_mass["Value"] * self.Units.M_SUN_SI * self.Units.GR_MASS_TO_METER) * tan(self.Observer_FOV["Value"] / 2 / self.Units.RAD_TO_MICRO_AS), "Unit": "[M]"}
        self.Simulation_configurator.observer.Image_y_max = {"Value":  (self.Object_distance["Value"] * self.Units.PC_TO_METER) / (self.Simulation_configurator.object_mass["Value"] * self.Units.M_SUN_SI * self.Units.GR_MASS_TO_METER) * tan(self.Observer_FOV["Value"] / 2 / self.Units.RAD_TO_MICRO_AS), "Unit": "[M]"}
        self.Simulation_configurator.observer.Image_x_min = {"Value": -(self.Object_distance["Value"] * self.Units.PC_TO_METER) / (self.Simulation_configurator.object_mass["Value"] * self.Units.M_SUN_SI * self.Units.GR_MASS_TO_METER) * tan(self.Observer_FOV["Value"] / 2 / self.Units.RAD_TO_MICRO_AS), "Unit": "[M]"}
        self.Simulation_configurator.observer.Image_x_max = {"Value":  (self.Object_distance["Value"] * self.Units.PC_TO_METER) / (self.Simulation_configurator.object_mass["Value"] * self.Units.M_SUN_SI * self.Units.GR_MASS_TO_METER) * tan(self.Observer_FOV["Value"] / 2 / self.Units.RAD_TO_MICRO_AS), "Unit": "[M]"}
        
        self.Simulation_configurator.observer.Resolution_x = {"Value": resolution, "Unit": "[-]"}
        self.Simulation_configurator.observer.Resolution_y = {"Value": resolution, "Unit": "[-]"}
        
        """ =============== Kill the hotspot =============== """
        self.Simulation_configurator.hotspot_model.Density_scale_factor = {"Value": 0, "Unit": "[g/cm^3]"}
        
        self.Simulation_configurator.geodesic_integrator.Integrator_type = {"Value": "RK78_Fehlberg", "Unit": "[-]"}
        self.Simulation_configurator.geodesic_integrator.RK_abs_accuracy = {"Value": 1e-12, "Unit": "[-]"}
        self.Simulation_configurator.geodesic_integrator.RK_rel_accuracy = {"Value": 1e-12, "Unit": "[-]"}
        
        self.Simulation_configurator.rad_transfer_integrator.Integrator_type = {"Value": "RK78_Fehlberg", "Unit": "[-]"}
        self.Simulation_configurator.geodesic_integrator.max_stepsize = {"Value": 100, "Unit": "[-]"}
        
        self.Simulation_configurator.geodesic_integrator.Max_rel_step_increase = {"Value": 5, "Unit": "[-]"}
        
        """ =============== The simulation name and input file path =============== """
        self.Simulation_configurator.simulation_name = {"Value": central_object + "_Wormhole_a_{}_redshift_{}_obsf_{}".format(spin, redshift_f, int(obs_freq / 1e9)), "Unit": "[-]"}
        
        """ =============== The simulation output file path =============== """
        self.Simulation_configurator.file_manager.Output_file_directory = parent_directory + "Wormhole_sim_paper"
        
    def Run_and_evaluate_simulation(self):
        
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

for freq in [230e9, 345e9]:

    Sim_instance = Wormhole_simulation_cofigurator(spin = 0.9, 
                                                    redshift_f = 2,
                                                    inclination = 160 * pi / 180,
                                                    resolution = 1024,
                                                    central_object = "M87", 
                                                    obs_freq = freq)

    Sim_instance.Run_and_evaluate_simulation()