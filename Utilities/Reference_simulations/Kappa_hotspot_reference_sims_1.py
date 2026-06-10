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

class Hotspot_reference_sims:
    
    def __init__(self):
        
        """ These parameters correspond to the ones in table 1 of https://arxiv.org/pdf/2309.10053. """
        
        self.Units = Units_class()
        self.Simulation_configurator = Simulation_configurator()
        self.Simulation_configurator.average_emission_pitch_angle =  {"Value": 0, "Unit": "[-]"}
        self.Simulation_configurator.observer.Include_polarization = {"Value": 1, "Unit": "[-]"}
        self.Simulation_configurator.observer.Cam_rotation_angle =   {"Value": 0, "Unit": "[-]"}
        
        self.Simulation_configurator.min_image_order               = {"Value": 0, "Unit": "[-]"}
        self.Simulation_configurator.max_image_order               = {"Value": 0, "Unit": "[-]"}
        
        """ Central black hole setup"""
        self.Simulation_configurator.metric_parameters.Metric_type = {"Value": "Kerr",  "Unit": "[-]"}
        self.Simulation_configurator.object_mass                   = {"Value": 4.297e6, "Unit": "[M_sun]"}
        self.Simulation_configurator.metric_parameters.Spin        = {"Value": 0.000, "Unit": "[M]"}
        self.Object_distance                                       = {"Value": 8.277e3, "Unit": "[Pc]"}
        
        """ Kill the accretion disk """   
        self.Simulation_configurator.disk_model.Enabled_flag = {"Value": 0, "Unit": "[-]"}

        """ Hotspot setup """
        self.Simulation_configurator.hotspot_model.Density_scale_factor     = {"Value": 2e6, "Unit": "[g/cm^3]"}
        self.Simulation_configurator.hotspot_model.Temperature_scale_factor = {"Value": 1e11,  "Unit": "[K]"}
        
        self.Simulation_configurator.hotspot_model.Density_profile     = {"Value": "Gaussian", "Unit": "[-]"}
        self.Simulation_configurator.hotspot_model.Temperature_profile = {"Value": "Gaussian", "Unit": "[-]"}
        
        self.Simulation_configurator.hotspot_model.Density_spread = {"Value": 1, "Unit": "[-]"}
        self.Simulation_configurator.hotspot_model.Temperature_spread = {"Value": 1, "Unit": "[-]"}
        
        self.Simulation_configurator.hotspot_model.Temporal_spread   = {"Value": 60000000000, "Unit": "[GM/c^3]"}
        self.Simulation_configurator.hotspot_model.Magnetization     = {"Value": 1,       "Unit": "[-]"}
        self.Simulation_configurator.emission_models.Kappa           = {"Value": 4,       "Unit": "[-]"}
        self.Simulation_configurator.hotspot_model.Ensamble_type     = {"Value": "Kappa", "Unit": "[-]"}
        self.Simulation_configurator.hotspot_model.Distance          = {"Value": 8,       "Unit": "[M]"} 
        
        self.Simulation_configurator.hotspot_model.Mag_field_geometry       = {"Value": "Vertical", "Unit": "[-]"}
        self.Simulation_configurator.hotspot_model.Threshold_relative_density   = {"Value": 1e-5, "Unit": "[-]"}
        
        """ This value for the initial hotspot azimuth makes it appear on the anti-beaming size at t_obs = 0. This makes the light curve look nicer. """
        self.Simulation_configurator.hotspot_model.Azimuth           = {"Value": pi * 0.50,  "Unit": "[M]"} 
        
        """ Observer setup """
        self.Simulation_configurator.observer.Distance    = {"Value": 1e4,            "Unit": "[M]"}
        self.Simulation_configurator.observer.Inclination = {"Value": 150 * pi / 180, "Unit": "[Rad]"}
        self.Simulation_configurator.observer.Azimuth     = {"Value": 0,              "Unit": "[Rad]"}
        
        self.Simulation_configurator.observer.Obs_frequency = {"Value": self.Units.C_LIGHT_SI / 2.2e-6, "Unit": "[Hz]"}
        self.Observer_FOV = {"Value": 200, "Unit": "[micro-arcsec]"}
        
        self.Simulation_configurator.observer.Image_y_min = {"Value": -(self.Object_distance["Value"] * self.Units.PC_TO_METER) / (self.Simulation_configurator.object_mass["Value"] * self.Units.M_SUN_SI * self.Units.GR_MASS_TO_METER) * tan(self.Observer_FOV["Value"] / 2 / self.Units.RAD_TO_MICRO_AS), "Unit": "[M]"}
        self.Simulation_configurator.observer.Image_y_max = {"Value":  (self.Object_distance["Value"] * self.Units.PC_TO_METER) / (self.Simulation_configurator.object_mass["Value"] * self.Units.M_SUN_SI * self.Units.GR_MASS_TO_METER) * tan(self.Observer_FOV["Value"] / 2 / self.Units.RAD_TO_MICRO_AS), "Unit": "[M]"}
        self.Simulation_configurator.observer.Image_x_min = {"Value": -(self.Object_distance["Value"] * self.Units.PC_TO_METER) / (self.Simulation_configurator.object_mass["Value"] * self.Units.M_SUN_SI * self.Units.GR_MASS_TO_METER) * tan(self.Observer_FOV["Value"] / 2 / self.Units.RAD_TO_MICRO_AS), "Unit": "[M]"}
        self.Simulation_configurator.observer.Image_x_max = {"Value":  (self.Object_distance["Value"] * self.Units.PC_TO_METER) / (self.Simulation_configurator.object_mass["Value"] * self.Units.M_SUN_SI * self.Units.GR_MASS_TO_METER) * tan(self.Observer_FOV["Value"] / 2 / self.Units.RAD_TO_MICRO_AS), "Unit": "[M]"}
        
        self.Simulation_configurator.observer.Resolution_x = {"Value": 256, "Unit": "[-]"}
        self.Simulation_configurator.observer.Resolution_y = {"Value": 256, "Unit": "[-]"}
        
        """ Configure the integrator """
        
        # self.Simulation_configurator.geodesic_integrator.Integrator_type = {"Value": "RK54", "Unit": "[-]"}
        # self.Simulation_configurator.emission_integrator.Rad_Transfer_Integrator_type = {"Value": "RK54", "Unit": "[-]"}
        # self.Simulation_configurator.emission_integrator.Parallel_Transport_Integrator_type = {"Value": "RK54", "Unit": "[-]"}
        self.Simulation_configurator.geodesic_integrator.RK_abs_accuracy = {"Value": 1e-12, "Unit": "[-]"}
        self.Simulation_configurator.geodesic_integrator.RK_rel_accuracy = {"Value": 1e-12, "Unit": "[-]"}
            
        self.Simulation_configurator.geodesic_integrator.max_upper_stepsize = {"Value": 100, "Unit": "[-]"}
        self.Simulation_configurator.geodesic_integrator.Max_rel_step_increase = {"Value": 10, "Unit": "[-]"}
        self.Simulation_configurator.geodesic_integrator.Max_step_in_emission_medium = {"Value": 0.01, "Unit": "[-]"}
        # self.Simulation_configurator.geodesic_integrator.Max_step_in_emission_medium = {"Value": 0.05, "Unit": "[-]"}
        
        """ The simulation output file path """
        self.Simulation_configurator.file_manager.Output_file_directory = parent_directory + "Reference_simulations"

    def run_simulation(self, obs_time: float, idx: int):
        
        self.Simulation_configurator.simulation_name = {"Value": "Hotspot_Reference_Simulation_{}".format(idx), "Unit": "[-]"}
        
        """ This observation time offset is to synch the hotspot temporal profile with its azimuth coordinate and get the maximum emission right at
            the beaming point (on the left side of the image). """
        self.Simulation_configurator.observer.Init_time = {"Value": obs_time - 80, "Unit": "[GM/c^3]"}
            
        self.Simulation_configurator.generate_simulation_input(Path_to_input_dir = "Reference_simulations\\Hotspot_Reference_Simulation_{}".format(idx),
                                                               Input_file_name = "Hotspot_Reference_Simulation_input.XML")
            
        """ Run the simulation """
        filename = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\Utilities\\Reference_simulations\\Hotspot_Reference_Simulation_{}\\Hotspot_Reference_Simulation_input.xml".format(idx)
        args = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\x64\\Release\\Mjolnir_GRRT.exe -in " + filename + " -print_to_console 0"
        subprocess.call(args, shell = True)          

if __name__ == "__main__":

    Hotspot_reference_sims_instance = Hotspot_reference_sims()

    Spot_period: float = 2 * pi * Hotspot_reference_sims_instance.Simulation_configurator.hotspot_model.Distance["Value"]**(3 / 2)
    # Spot_period: float = 2 * pi / 0.043801516696240203
     
    Hotspot_number: int = 20
    
    Obs_times: NDArray = linspace(0, Spot_period, Hotspot_number)
    File_idx: NDArray  = linspace(0, Hotspot_number - 1, Hotspot_number)

    Sim_args: list[list[float | int]] = []

    for t_obs, idx in zip(Obs_times, File_idx):
        Sim_args.append([t_obs, int(idx)])
        
    with Pool(10) as pool:
        pool.starmap(Hotspot_reference_sims_instance.run_simulation, Sim_args)