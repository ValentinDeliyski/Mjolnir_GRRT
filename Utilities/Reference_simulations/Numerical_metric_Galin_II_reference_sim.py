import sys
import os
import threading
from time import sleep

""" Add the parent directory of this file to the search path, 
    so this script can be ran from the "Utilities" folder """
    
parent_directory = os.path.abspath('...')
sys.path.append(parent_directory)

from Mjolnir_Configurator import Simulation_configurator
from numpy import pi
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
        
        self.Sim_config = Simulation_configurator()
        
        self.Sim_config.metric_parameters.Numerical_metric_spline_path = "C:/Users/Valur/Documents/Repos/Mjolnir_GRRT/Utilities/Galin_numerical_config_II.XML"

        self.Sim_config.simulation_mode = {"Value": 1, "Unit": "[-]"}

        self.Sim_config.object_mass = {"Value": 6.2e9, "Unit": "[M_sun]"}

        # ================================================== Metric ================================================== #

        self.Sim_config.metric_parameters.Metric_type    = {"Value": "Numerical", "Unit": "[-]"}
        self.Sim_config.metric_parameters.Mass           = {"Value": 0.881990876889021, "Unit": "[M]"}
        self.Sim_config.metric_parameters.Horizon_radius = {"Value": 0.01, "Unit": "[G/c^2]"}
        self.Sim_config.metric_parameters.Spin           = {"Value": 0.7258549181881021 / 0.881990876889021, "Unit": "[M]"}
        self.Sim_config.metric_parameters.Numerical_metric_anzatz_type = {"Value": "Anzatz_1", "Unit": "[M]"} 
        
        self.Sim_config.metric_parameters.Scattering_radius = {"Value": 30, "Unit": "[M]"} 
        
        # ================================================== Observer ================================================== #

        self.Sim_config.observer.Resolution_x = {"Value": 1024, "Unit": "[-]"}
        self.Sim_config.observer.Resolution_y = {"Value": 1024, "Unit": "[-]"}
        
        self.Sim_config.observer.Distance    = {"Value": 1e4, "Unit": "[M]"}
        self.Sim_config.observer.Inclination = {"Value": 90 * pi / 180, "Unit": "[Rad]"}
        self.Sim_config.observer.Obs_frequency = {"Value": 230e9, "Unit": "[Hz]"}
        self.Sim_config.observer.Cam_rotation_angle = {"Value": 0, "Unit": "[Hz]"}

        # ================================================== Disk ================================================== #
        
        self.Sim_config.disk_model.Ensamble_type = {"Value": "Thermal",   "Unit": "[-]"}
        self.Sim_config.disk_model.Disk_Model    = {"Value": "Phenom_RIAF_1", "Unit": "[-]"}
        self.Sim_config.disk_model.Mag_field_geometry = {"Value": "Constant", "Unit": "[-]"}
        
        self.Sim_config.disk_model.Density_scale_factor = {"Value": 500000, "Unit": "[g/cm^3]"}
        self.Sim_config.disk_model.Temperature_scale_factor = {"Value": 4.1e+10, "Unit": "[K]"}
                
        self.Sim_config.disk_model.Density_cutoff_radius = {"Value": 5, "Unit": "[M]"}
        self.Sim_config.disk_model.Temperature_cutoff_radius = {"Value": 5, "Unit": "[M]"}

        self.Sim_config.disk_model.Density_power_law_scale     = {"Value": 5, "Unit": "[M]"}
        self.Sim_config.disk_model.Temperature_power_law_scale = {"Value": 5, "Unit": "[M]"}

        self.Sim_config.disk_model.Opening_angle = {"Value": 0.4, "Unit": "[tan(angle)]"}
        
        self.Sim_config.disk_model.Density_power_law_power     = {"Value": 2.0, "Unit": "[-]"}
        self.Sim_config.disk_model.Temperature_power_law_power = {"Value": 1.0, "Unit": "[-]"}
        
        self.Sim_config.disk_model.Velocity_profile = {"Value": "Theta Dependant", "Unit": "[-]"}
        self.Sim_config.integrator.RK78_abs_accuracy  = {"Value": 1e-10, "Unit": "[-]"}
        self.Sim_config.integrator.RK78_rel_accuracy  = {"Value": 1e-10, "Unit": "[-]"}
        self.Sim_config.observer.Include_polarization = {"Value": 0, "Unit": "[-]"}
        
        # self.Sim_config.integrator.Step_controller_type  = {"Value": "PID", "Unit": "[-]"}
        self.Sim_config.integrator.max_integration_count = {"Value": 1000000, "Unit": "[-]"}
        self.Sim_config.integrator.max_affine_parameter  = {"Value": 100000, "Unit": "[-]"}
        
        self.Sim_config.metric_parameters.Distance_to_singular_point = {"Value": 1e-3, "Unit": "[M]"}
        
        self.Sim_config.integrator.Max_rel_step_increase = {"Value": 1.01, "Unit": "[-]"}
        self.Sim_config.integrator.Min_rel_step_increase = {"Value": 0, "Unit": "[-]"}
        self.Sim_config.integrator.init_stepsize         = {"Value": 1, "Unit": "[-]"}
        
        self.Sim_config.integrator.max_stepsize         = {"Value": 100, "Unit": "[-]"}
        
        self.Sim_config.integrator.radiative_transfer_integrator_type = {"Value": "Implicit Trapezoid", "Unit": "[-]"}
        # ================================================== Hotspot ================================================== #

        self.Sim_config.hotspot_model.Density_scale_factor = {"Value": 0, "Unit": "[g / cm^3]"}

    def run_simulation(self):
            
        self.Sim_config.observer.Image_y_min = {"Value": -2.0, "Unit": "[M]"}
        self.Sim_config.observer.Image_y_max = {"Value":  10.0, "Unit": "[M]"}
        self.Sim_config.observer.Image_x_min = {"Value": -10.0, "Unit": "[M]"}
        self.Sim_config.observer.Image_x_max = {"Value":  10.0, "Unit": "[M]"}
            
        """ The simulation output file path """
        self.Sim_config.file_manager.Output_file_directory = parent_directory + "Reference_simulations"
        self.Sim_config.simulation_name = {"Value": "Galin_Numerical_config_II_4", "Unit": "[-]"}
        
        self.Sim_config.generate_simulation_input(Path_to_input_dir = "Reference_simulations\\Galin_Numerical_config_II_4",
                                                  Input_file_name = "Galin_Numerical_config_II_4.XML")
        
        filename = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\Utilities\\Reference_simulations\\Galin_Numerical_config_II_4\\Galin_Numerical_config_II_4.xml"
        args = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\x64\\Release\\Mjolnir_GRRT.exe -in " + filename + " -print_to_console 1"
           
        subprocess.call(args, shell = True)
        
if __name__ == "__main__":

    Simulation_instance = Simulation()

    Sim_1_thread = threading.Thread(target = Simulation_instance.run_simulation, args = [])
    Sim_1_thread.start()


    