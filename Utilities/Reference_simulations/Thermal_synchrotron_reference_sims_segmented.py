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
from numpy import pi, tan, sqrt, array, flip
import matplotlib.pyplot as plt
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
        self.Simulation_configurator.disk_model.Temperature_scale_factor = {"Value": 1e11, "Unit": "[K]"}
        
        self.Simulation_configurator.disk_model.Magnetization = {"Value": 0.01, "Unit": "[-]"}
        
        self.Simulation_configurator.disk_model.Density_cutoff_radius     = {"Value": 0, "Unit": "[-]"}
        self.Simulation_configurator.disk_model.Temperature_cutoff_radius = {"Value": 0, "Unit": "[-]"}
        
        self.Simulation_configurator.disk_model.Ensamble_type = {"Value": "Thermal",   "Unit": "[-]"}
        self.Simulation_configurator.disk_model.Disk_Model    = {"Value": "Phenom_RIAF_1", "Unit": "[-]"}
        
        self.Simulation_configurator.disk_model.Density_power_law_power     = {"Value": 2.0, "Unit": "[-]"}
        self.Simulation_configurator.disk_model.Temperature_power_law_power = {"Value": 1.0, "Unit": "[-]"}
        
        self.Simulation_configurator.disk_model.Velocity_profile = {"Value": "Theta Dependant", "Unit": "[-]"}
        
        """ Observer setup """
        self.Simulation_configurator.observer.Distance    = {"Value": 1e4,            "Unit": "[M]"}
        self.Simulation_configurator.observer.Inclination = {"Value": 160 * pi / 180, "Unit": "[Rad]"}
        self.Simulation_configurator.observer.Azimuth     = {"Value": 0,              "Unit": "[Rad]"}
        
        self.Simulation_configurator.observer.Obs_frequency = {"Value": 230e9, "Unit": "[Hz]"}
        self.Observer_FOV = {"Value": 100, "Unit": "[micro-arcsec]"}

        self.Nominal_Image_x_max = (self.Object_distance["Value"] * self.Units.PC_TO_METER) / (self.Simulation_configurator.object_mass["Value"] * self.Units.M_SUN_SI * self.Units.GR_MASS_TO_METER) * tan(self.Observer_FOV["Value"] / 2 / self.Units.RAD_TO_MICRO_AS) 
        self.Nominal_Image_y_max = (self.Object_distance["Value"] * self.Units.PC_TO_METER) / (self.Simulation_configurator.object_mass["Value"] * self.Units.M_SUN_SI * self.Units.GR_MASS_TO_METER) * tan(self.Observer_FOV["Value"] / 2 / self.Units.RAD_TO_MICRO_AS) 

        self.Nominal_Resolution_x = 512
        self.Nominal_Resolution_y = 512

        self.Segment_number = 8
        
        """ Kill the hotspot """
        self.Simulation_configurator.hotspot_model.Density_scale_factor = {"Value": 0, "Unit": "[g/cm^3]"}
    
        self.Simulation_configurator.integrator.RK78_abs_accuracy  = {"Value": 1e-13, "Unit": "[-]"}
        self.Simulation_configurator.integrator.RK78_rel_accuracy  = {"Value": 1e-13, "Unit": "[-]"}
        self.Simulation_configurator.observer.Include_polarization = {"Value": 0, "Unit": "[-]"}

        """ The simulation output file path """
        self.Simulation_configurator.file_manager.Output_file_directory = parent_directory + "Reference_simulations"

    def Run_simulation(self, Segment_idx: int) -> None:
        
        """ This simulation corresponds to the "low spin circ/thin" one in https://arxiv.org/pdf/2206.12066. Their results are presented on the top 4 panels
            of figure 3. They find a total flux of 546 mJy """
        
        """ Central black hole setup """
        self.Simulation_configurator.metric_parameters.Spin = {"Value": 0.01, "Unit": "[M]"}

        """ Accretion disk setup """
        self.Simulation_configurator.disk_model.Opening_angle = {"Value": 0.1, "Unit": "[tan(angle)]"}
        self.Simulation_configurator.disk_model.Density_power_law_scale     = {"Value": 1 + sqrt(1 - self.Simulation_configurator.metric_parameters.Spin["Value"]**2), "Unit": "[M]"}
        self.Simulation_configurator.disk_model.Temperature_power_law_scale = {"Value": 1 + sqrt(1 - self.Simulation_configurator.metric_parameters.Spin["Value"]**2), "Unit": "[M]"}
        
        self.Simulation_configurator.disk_model.Density_scale_factor = {"Value": 1.5e6, "Unit": "[g/cm^3]"}

        """ The segmented observation window """
        Nominal_y_scan_step = 2 * self.Nominal_Image_y_max / (self.Nominal_Resolution_y - 1)

        self.Simulation_configurator.observer.Image_y_max = {"Value": -self.Nominal_Image_y_max + (Segment_idx + 1) * (2 * self.Nominal_Image_y_max / self.Segment_number), "Unit": "[M]"}
        self.Simulation_configurator.observer.Image_y_min = {"Value": -self.Nominal_Image_y_max + (Segment_idx + 0) * (2 * self.Nominal_Image_y_max / self.Segment_number), "Unit": "[M]"}

        if Segment_idx > 0 and Segment_idx < self.Segment_number - 1:
            self.Simulation_configurator.observer.Image_y_max["Value"] -= Nominal_y_scan_step / 2
            self.Simulation_configurator.observer.Image_y_min["Value"] += Nominal_y_scan_step / 2

        elif Segment_idx == self.Segment_number - 1:
            self.Simulation_configurator.observer.Image_y_min["Value"] += Nominal_y_scan_step / 2

        elif Segment_idx == 0:
            self.Simulation_configurator.observer.Image_y_max["Value"] -= Nominal_y_scan_step / 2

        self.Simulation_configurator.observer.Image_x_max = {"Value":  self.Nominal_Image_x_max, "Unit": "[M]"}
        self.Simulation_configurator.observer.Image_x_min = {"Value": -self.Nominal_Image_x_max, "Unit": "[M]"}

        self.Simulation_configurator.observer.Resolution_x = {"Value": self.Nominal_Resolution_x, "Unit": "[-]"}
        self.Simulation_configurator.observer.Resolution_y = {"Value": int(self.Nominal_Resolution_y / self.Segment_number), "Unit": "[-]"}
        
        """ The simulation name and input file path """
        self.Simulation_configurator.simulation_name = {"Value": "Reference_Simulation_1_Segment_{}".format(Segment_idx), "Unit": "[-]"}
        
        self.Simulation_configurator.generate_simulation_input(Path_to_input_dir = "Reference_simulations\\Reference_Simulation_1_Segment_{}".format(Segment_idx),
                                                               Input_file_name = "Reference_Simulation_1_Segment_{}_input.XML".format(Segment_idx))
        
        """ Run the simulation """
        filename = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\Utilities\\Reference_simulations\\Reference_Simulation_1_Segment_{}\\Reference_Simulation_1_Segment_{}_input.xml".format(Segment_idx, Segment_idx)
        args = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\x64\\Release\\Mjolnir_GRRT.exe -in " + filename + " -print_to_console 0"
        subprocess.call(args, shell = True)

    def Combine_segmented_output_files(self):

        Total_image_intensity = []

        for Segment_idx in range(self.Segment_number):

            Output_file_path = self.Simulation_configurator.file_manager.Output_file_directory + "\\Reference_Simulation_1_Segment_{}\\Kerr".format(Segment_idx) 
            Sim_parser = Simulation_Parser(Output_file_path)
            Total_image_intensity.append(Sim_parser.I_Intensity)

        Total_image_intensity = array(Total_image_intensity)
        Total_image_intensity = Total_image_intensity.flatten().reshape(self.Nominal_Resolution_x, self.Nominal_Resolution_y)

        Total_image_intensity = flip(Total_image_intensity, axis = 0)

        Fig = plt.figure().add_subplot(111)
        Fig.imshow(Total_image_intensity, cmap="hot")

        plt.show()

if __name__ == "__main__":

    Thermal_syhnchrotron_reference_sims_instance = Thermal_syhnchrotron_reference_sims()

    Thermal_syhnchrotron_reference_sims_instance.Segment_number = 8

    Segment_list = []

    for Segment_idx in range(Thermal_syhnchrotron_reference_sims_instance.Segment_number):
        Segment_list.append(threading.Thread(target = Thermal_syhnchrotron_reference_sims_instance.Run_simulation, args = [Segment_idx]))

    for Segment_idx in range(Thermal_syhnchrotron_reference_sims_instance.Segment_number):
        Segment_list[Segment_idx].start()

    for Segment_idx in range(Thermal_syhnchrotron_reference_sims_instance.Segment_number):
        Segment_list[Segment_idx].join()

    Thermal_syhnchrotron_reference_sims_instance.Combine_segmented_output_files()
