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
from numpy import pi, tan, sqrt, array, flip, append, float64
from numpy.typing import NDArray
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

        self.Nominal_Resolution_x = 2048
        self.Nominal_Resolution_y = 2048

        self.Segment_number = 8
        
        """ Kill the hotspot """
        self.Simulation_configurator.hotspot_model.Density_scale_factor = {"Value": 0, "Unit": "[g/cm^3]"}
        
        self.Simulation_configurator.geodesic_integrator.Integrator_type = {"Value": "RK78_Fehlberg", "Unit": "[-]"}
        self.Simulation_configurator.geodesic_integrator.RK_abs_accuracy = {"Value": 1e-12, "Unit": "[-]"}
        self.Simulation_configurator.geodesic_integrator.RK_rel_accuracy = {"Value": 1e-12, "Unit": "[-]"}
        
        self.Simulation_configurator.rad_transfer_integrator.Integrator_type = {"Value": "RK78_Fehlberg", "Unit": "[-]"}
        self.Simulation_configurator.rad_transfer_integrator.RK_abs_accuracy = {"Value": 1e-10, "Unit": "[-]"}
        self.Simulation_configurator.rad_transfer_integrator.RK_rel_accuracy = {"Value": 1e-10, "Unit": "[-]"}
            
        self.Simulation_configurator.geodesic_integrator.max_stepsize = {"Value": 100, "Unit": "[-]"}
        
        self.Simulation_configurator.geodesic_integrator.Max_rel_step_increase = {"Value": 2, "Unit": "[-]"}
        self.Simulation_configurator.rad_transfer_integrator.Max_rel_step_increase = {"Value": 2, "Unit": "[-]"}
        
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
        
        X_coords: NDArray[float64]      = array([])
        Y_coords: NDArray[float64]      = array([])
        I_Intensity: NDArray[float64]   = array([])
        Q_Intensity: NDArray[float64]   = array([])
        U_Intensity: NDArray[float64]   = array([])
        V_Intensity: NDArray[float64]   = array([])
        Disk_redshift: NDArray[float64] = array([])
        Disk_flux: NDArray[float64]     = array([])
        Celestial_theta: NDArray[float64] = array([])
        Celestial_phi: NDArray[float64]   = array([])

        Raw_simulation_header: str = ""
        PT_model_active: bool = False
        
        for Segment_idx in range(self.Segment_number):

            Output_file_path = self.Simulation_configurator.file_manager.Output_file_directory + "\\Reference_Simulation_1_Segment_{}\\Kerr".format(Segment_idx) 
            Sim_parser = Simulation_Parser(Output_file_path)
        
            PT_model_active = Sim_parser.Simulation_metadata["Active disk model"] == "Page-Thorne"
            
            X_coords = append(X_coords, Sim_parser.X_coords)
            Y_coords = append(Y_coords, Sim_parser.Y_coords)
            I_Intensity = append(I_Intensity, Sim_parser.I_Intensity)
            Q_Intensity = append(Q_Intensity, Sim_parser.Q_Intensity)
            U_Intensity = append(U_Intensity, Sim_parser.U_Intensity)
            V_Intensity = append(V_Intensity, Sim_parser.V_Intensity)
            Disk_redshift = append(Disk_redshift, Sim_parser.Disk_redshift)
            Disk_flux = append(Disk_flux, Sim_parser.Disk_flux)
            Celestial_theta = append(Celestial_theta, Sim_parser.Celestial_theta)
            Celestial_phi = append(Celestial_phi, Sim_parser.Celestial_phi)
            
            Raw_simulation_header = Sim_parser.Raw_simulation_header
            
        """ ============================ Flatten all arrays, so I can write in the combined file easily ============================ """
            
        X_coords = X_coords.flatten()
        Y_coords = Y_coords.flatten()
        I_Intensity = I_Intensity.flatten()
        Q_Intensity = Q_Intensity.flatten()
        U_Intensity = U_Intensity.flatten()
        V_Intensity = V_Intensity.flatten()
        Disk_redshift = Disk_redshift.flatten()
        Disk_flux = Disk_flux.flatten()
        Celestial_theta = Celestial_theta.flatten()
        Celestial_phi = Celestial_phi.flatten()
            
        if not os.path.isdir(self.Simulation_configurator.file_manager.Output_file_directory + "\\Reference_Simulation_1"):
            os.mkdir(self.Simulation_configurator.file_manager.Output_file_directory + "\\Reference_Simulation_1")
        
        with open(self.Simulation_configurator.file_manager.Output_file_directory + "\\Reference_Simulation_1\\" + self.Simulation_configurator.simulation_name["Value"].split("_Segment")[0] + ".txt", "w") as file:
            file.write(Raw_simulation_header)
            
            if not PT_model_active:
            
                file.write("Image X Coord [M],Image Y Coord [M],Synchotron Intensity I [Jy/sRad],Synchotron Intensity Q [Jy/sRad],Synchotron Intensity U [Jy/sRad],Synchotron Intensity V [Jy/sRad],Celestial Sphere Crossing Theta [Rad],Celestial Sphere Crossing Phi [Rad]\n")
            
                for X_coord, Y_coord, I, Q, U, V, Theta, Phi in zip(X_coords, Y_coords, I_Intensity, Q_Intensity, U_Intensity, V_Intensity, Celestial_theta, Celestial_phi):
                    file.write("{},{},{},{},{},{},{}, {}\n".format(X_coord, Y_coord, I, Q, U, V, Theta, Phi))

        print("Segments exported to single file in " + self.Simulation_configurator.file_manager.Output_file_directory + "\\Reference_Simulation_1")

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
