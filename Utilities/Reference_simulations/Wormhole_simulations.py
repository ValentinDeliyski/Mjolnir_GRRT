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
        self.Simulation_configurator = Simulation_configurator()
        self.central_object = central_object
        
        """ Central Wormhole setup"""
        
        self.Simulation_configurator.metric_parameters.Metric_type = {"Value": "Wormhole", "Unit": "[-]"}
        self.Simulation_configurator.metric_parameters.WH_redshift = {"Value": redshift_f, "Unit": "[M]"}
        self.Simulation_configurator.metric_parameters.Spin        = {"Value":    spin,    "Unit": "[M]"}
        
        match self.central_object:
            case "M87":
                self.Simulation_configurator.object_mass = {"Value": self.Units.M_M87_BH_SI / self.Units.M_SUN_SI,  "Unit": "[M_sun]"}
                self.Object_distance = {"Value": self.Units.M87_DISTANCE_PC, "Unit": "[Pc]"}
                
                r_0 = 5
                T_max = 5.1e10
                n_max = 5e5
                
        
            case "Sgr_A":
                self.Simulation_configurator.object_mass = {"Value": self.Units.M_SGRA_BH_SI / self.Units.M_SUN_SI,  "Unit": "[M_sun]"}
                self.Object_distance = {"Value": self.Units.SGRA_DISTANCE_PC, "Unit": "[Pc]"}
                
                r_0 = 2
                T_max = 1e11
                n_max = 5e6
                
                self.Simulation_configurator.disk_model.Density_cutoff_scale = {"Value": 0.4 * 9.75684 / 13.8033 , "Unit": "[-]"}
                
            case _:
                
                print("Wrong central object!")
                
                exit(1)
        
        """ =============== Accretion disk setup =============== """
        
        """ ----------- Density ----------- """
        
        self.Simulation_configurator.disk_model.Density_radial_power_law = {"Value":     2.0,     "Unit":     "[-]"     }
        self.Simulation_configurator.disk_model.Density_scale_factor     = {"Value":    n_max,    "Unit":   "[g/cm^3]"  }
        self.Simulation_configurator.disk_model.Density_r_cutoff         = {"Value":     r_0,     "Unit":     "[-]"     }
        self.Simulation_configurator.disk_model.Density_profile          = {"Value": "Power Law", "Unit":     "[-]"     }
        self.Simulation_configurator.disk_model.Opening_angle            = {"Value":     0.1,     "Unit": "[tan(angle)]"}
        self.Simulation_configurator.disk_model.Density_r_0              = {"Value":     r_0,     "Unit":     "[M]"     }
        
        """ --------- Temperature --------- """
        
        self.Simulation_configurator.disk_model.Temperature_radial_power_law = {"Value":     1.0,     "Unit": "[-]"}
        self.Simulation_configurator.disk_model.Temperature_scale_factor     = {"Value":    T_max,    "Unit": "[K]"}
        self.Simulation_configurator.disk_model.Temperature_r_cutoff         = {"Value":     r_0,     "Unit": "[-]"}
        self.Simulation_configurator.disk_model.Temperature_profile          = {"Value": "Power Law", "Unit": "[-]"}
        self.Simulation_configurator.disk_model.Temperature_r_0              = {"Value":     r_0,     "Unit": "[M]"}
        
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
        
        """ =============== Kill the Novikov-Thorne disk =============== """
        self.Simulation_configurator.NT_model_params.Evaluate_NT_disk = {"Value": 0, "Unit": "[-]"}
    
        self.Simulation_configurator.integrator.RK45_accuracy      = {"Value": 1e-13, "Unit": "[-]"}
        self.Simulation_configurator.observer.Include_polarization = {"Value": 0, "Unit": "[-]"}
        
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
        
        args     = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\x64\\Release\\Mjolnir_GRRT.exe -in " + filename + " -print_to_console 0"
        subprocess.call(args, shell = True)
                
        if self.central_object == "M87":
            Distance = self.Units.M87_DISTANCE_GEOMETRICAL
        else:
            Distance = self.Units.SGRA_DISTANCE_GEOMETRICAL
                
        """ Evaluate the simulataion results """
        Sim_parser_n0 = Simulation_Parser(parent_directory + "Wormhole_sim_paper\\"+ self.Simulation_configurator.simulation_name["Value"] + "\\Wormhole_n0")
        Total_flux_n0 = Sim_parser_n0.get_total_flux(obs_pos = Distance, unit = "mJy")
        
        Sim_parser_n1 = Simulation_Parser(parent_directory + "Wormhole_sim_paper\\"+ self.Simulation_configurator.simulation_name["Value"] + "\\Wormhole_n1")
        Total_flux_n1 = Sim_parser_n1.get_total_flux(obs_pos = Distance, unit = "mJy")
        
        Sim_parser_n2 = Simulation_Parser(parent_directory + "Wormhole_sim_paper\\"+ self.Simulation_configurator.simulation_name["Value"] + "\\Wormhole_n2")
        Total_flux_n2 = Sim_parser_n2.get_total_flux(obs_pos = Distance, unit = "mJy")
        
        Sim_parser_n3 = Simulation_Parser(parent_directory + "Wormhole_sim_paper\\"+ self.Simulation_configurator.simulation_name["Value"] + "\\Wormhole_n3")
        Total_flux_n3 = Sim_parser_n3.get_total_flux(obs_pos = Distance, unit = "mJy")
        
        Total_flux = Total_flux_n0 + Total_flux_n1 + Total_flux_n2 + Total_flux_n3
        
        print("Finished simulation {}!".format(self.Simulation_configurator.simulation_name["Value"]))
        print("Total Flux [mJy] = {}".format(round(Total_flux, 2)))
        
Sim_threads = []
        
for central_object in ["Sgr_A"]:
        
    for obs_freq in [230e9]:
        
        for redshift_values in [0]:
            
            for spin_values in [0.99]:

                Sim_threads.append(threading.Thread(target = Wormhole_simulation_cofigurator(spin = spin_values,
                                                                                             redshift_f = redshift_values,
                                                                                             inclination = 160 * pi / 180,
                                                                                             resolution = 2048,
                                                                                             central_object = central_object, 
                                                                                             obs_freq = obs_freq).Run_and_evaluate_simulation, args = []))

# Thread legend
# 0  = not started
# 1  = running
# -1 = finished

Running_threads_mask = array([0 for _ in Sim_threads])
Max_thread_number = len(Sim_threads)

for Thred_idx, Thread in enumerate(Sim_threads):
    
        while sum(Running_threads_mask == 1) >= 9:
            
            time.sleep(5)
            
            running_threads_idx_list, = where(Running_threads_mask == 1)
            
            for idx in running_threads_idx_list:
                
                if not Sim_threads[idx].is_alive():
                    
                    Running_threads_mask[idx] = -1
                    
        Thread.start()
                        
        time.sleep(1)
            
        Running_threads_mask[Thred_idx] = 1