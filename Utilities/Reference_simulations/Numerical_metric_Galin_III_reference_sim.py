import sys
import os
import threading
from time import sleep

""" Add the parent directory of this file to the search path, 
    so this script can be ran from the "Utilities" folder """
    
parent_directory = os.path.abspath('...')
sys.path.append(parent_directory)

from Mjolnir_Configurator import Simulation_configurator
from Support_functions.Parsers import Units_class, Simulation_Parser
from numpy import pi, tan, sqrt, array, flip, append, float64, arctan
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


class Simulation:
    
    def __init__(self):
        
        self.Sim_config = Simulation_configurator()
        
        self.Sim_config.metric_parameters.Numerical_metric_spline_path = "C:/Users/Valur/Documents/Repos/Mjolnir_GRRT/Utilities/Galin_numerical_config_III.XML"

        self.Sim_config.simulation_mode = {"Value": 1, "Unit": "[-]"}

        self.Sim_config.object_mass = {"Value": 6.2e9, "Unit": "[M_sun]"}

        # ================================================== Metric ================================================== #

        self.Sim_config.metric_parameters.Metric_type    = {"Value": "Numerical", "Unit": "[-]"}
        self.Sim_config.metric_parameters.Mass           = {"Value": 0.8904892552349474, "Unit": "[M]"}
        self.Sim_config.metric_parameters.Horizon_radius = {"Value": 0.01, "Unit": "[G/c^2]"}
        self.Sim_config.metric_parameters.Spin           = {"Value": 0.7813066738190858 / 0.8904892552349474, "Unit": "[M]"}
        self.Sim_config.metric_parameters.Numerical_metric_anzatz_type = {"Value": "Anzatz_1", "Unit": "[M]"} 
        
        self.Sim_config.metric_parameters.Scattering_radius = {"Value": 25.81413133692742, "Unit": "[M]"} 
        
        # ================================================== Observer ================================================== #
        
        self.Sim_config.observer.Distance    = {"Value": 12.439795771152403, "Unit": "[M]"}
        self.Sim_config.observer.Inclination = {"Value": 90 * pi / 180, "Unit": "[Rad]"}
        self.Sim_config.observer.Obs_frequency = {"Value": 230e9, "Unit": "[Hz]"}
        self.Sim_config.observer.Cam_rotation_angle = {"Value": 0, "Unit": "[Hz]"}
        
        self.Sim_config.observer.Use_angular_coords = {"Value": 1, "Unit": "[M]"}
        
        self.Nominal_Image_x_max = 10
        self.Nominal_Image_y_max = 10
        
        self.Nominal_Image_x_angle_max = arctan(self.Nominal_Image_x_max / self.Sim_config.observer.Distance["Value"])
        self.Nominal_Image_y_angle_max = arctan(self.Nominal_Image_y_max / self.Sim_config.observer.Distance["Value"])
        
        self.Nominal_Resolution_x = 1024
        self.Nominal_Resolution_y = 1024

        self.Segment_number = 8
        
        # ================================================== Disk ================================================== #
        
        self.Sim_config.disk_model.Ensamble_type = {"Value": "Thermal",   "Unit": "[-]"}
        self.Sim_config.disk_model.Disk_Model    = {"Value": "Phenom_RIAF_1", "Unit": "[-]"}
        self.Sim_config.disk_model.Mag_field_geometry = {"Value": "Constant", "Unit": "[-]"}
        
        self.Sim_config.disk_model.Density_scale_factor = {"Value": 500000, "Unit": "[g/cm^3]"}
        self.Sim_config.disk_model.Temperature_scale_factor = {"Value": 4.8e+10, "Unit": "[K]"}
                
        self.Sim_config.disk_model.Density_cutoff_radius = {"Value": 5, "Unit": "[M]"}
        self.Sim_config.disk_model.Temperature_cutoff_radius = {"Value": 5, "Unit": "[M]"}

        self.Sim_config.disk_model.Density_power_law_scale     = {"Value": 5, "Unit": "[M]"}
        self.Sim_config.disk_model.Temperature_power_law_scale = {"Value": 5, "Unit": "[M]"}

        self.Sim_config.disk_model.Opening_angle = {"Value": 0.1, "Unit": "[tan(angle)]"}
        
        self.Sim_config.disk_model.Density_power_law_power     = {"Value": 2.0, "Unit": "[-]"}
        self.Sim_config.disk_model.Temperature_power_law_power = {"Value": 1.0, "Unit": "[-]"}
        
        self.Sim_config.disk_model.Velocity_profile = {"Value": "Theta Dependant", "Unit": "[-]"}

        self.Sim_config.integrator.RK78_abs_accuracy = {"Value": 1e-14, "Unit": "[-]"}
        self.Sim_config.integrator.RK78_rel_accuracy = {"Value": 1e-14, "Unit": "[-]"}
        
        self.Sim_config.integrator.ESDIRK54_abs_accuracy = {"Value": 5e-9, "Unit": "[-]"}
        self.Sim_config.integrator.ESDIRK54_rel_accuracy = {"Value": 5e-9, "Unit": "[-]"}
        self.Sim_config.observer.Include_polarization = {"Value": 0, "Unit": "[-]"}
        
        self.Sim_config.integrator.max_integration_count = {"Value": 1000000, "Unit": "[-]"}
        self.Sim_config.integrator.max_affine_parameter  = {"Value": 10000, "Unit": "[-]"}
        
        self.Sim_config.metric_parameters.Distance_to_singular_point = {"Value": 1e-3, "Unit": "[M]"}
        
        self.Sim_config.integrator.Max_rel_step_increase = {"Value": 5, "Unit": "[-]"}
        self.Sim_config.integrator.Min_rel_step_increase = {"Value": 0, "Unit": "[-]"}
        self.Sim_config.integrator.init_stepsize         = {"Value": 1e-5, "Unit": "[-]"}
        
        self.Sim_config.integrator.max_stepsize = {"Value": 10, "Unit": "[-]"}
        self.Sim_config.integrator.default_geodesic_integrator_type = {"Value": "RK78_DP", "Unit": "[-]"}
        
        self.Sim_config.integrator.radiative_transfer_integrator_type = {"Value": "Implicit Trapezoid", "Unit": "[-]"}
        
        # ================================================== Hotspot ================================================== #

        self.Sim_config.hotspot_model.Density_scale_factor = {"Value": 0, "Unit": "[g / cm^3]"}

    def Run_simulation(self, Segment_idx: int):
            
        """ The segmented observation window """
        Nominal_y_scan_step = 2 * self.Nominal_Image_y_angle_max / (self.Nominal_Resolution_y - 1)
        
        self.Sim_config.observer.Image_y_angle_max = {"Value": -self.Nominal_Image_y_angle_max + (Segment_idx + 1) * (2 * self.Nominal_Image_y_angle_max / self.Segment_number), "Unit": "[M]"}
        self.Sim_config.observer.Image_y_angle_min = {"Value": -self.Nominal_Image_y_angle_max + (Segment_idx + 0) * (2 * self.Nominal_Image_y_angle_max / self.Segment_number), "Unit": "[M]"}

        if Segment_idx > 0 and Segment_idx < self.Segment_number - 1:
            self.Sim_config.observer.Image_y_angle_max["Value"] -= Nominal_y_scan_step / 2
            self.Sim_config.observer.Image_y_angle_min["Value"] += Nominal_y_scan_step / 2

        elif Segment_idx == self.Segment_number - 1:
            self.Sim_config.observer.Image_y_angle_min["Value"] += Nominal_y_scan_step / 2

        elif Segment_idx == 0:
            self.Sim_config.observer.Image_y_angle_max["Value"] -= Nominal_y_scan_step / 2
            

        self.Sim_config.observer.Image_x_angle_max = {"Value":  self.Nominal_Image_x_angle_max, "Unit": "[M]"}
        self.Sim_config.observer.Image_x_angle_min = {"Value": -self.Nominal_Image_x_angle_max, "Unit": "[M]"}

        self.Sim_config.observer.Resolution_x = {"Value": self.Nominal_Resolution_x, "Unit": "[-]"}
        self.Sim_config.observer.Resolution_y = {"Value": int(self.Nominal_Resolution_y / self.Segment_number), "Unit": "[-]"}
                 
        """ The simulation output file path """
        self.Sim_config.file_manager.Output_file_directory = parent_directory + "Reference_simulations"
        self.Sim_config.simulation_name = {"Value": "Galin_Numerical_config_III_5_345_Segment_{}".format(Segment_idx), "Unit": "[-]"}
        
        self.Sim_config.generate_simulation_input(Path_to_input_dir = "Reference_simulations\\Galin_Numerical_config_III_5_345_Segment_{}".format(Segment_idx),
                                                  Input_file_name = "Galin_Numerical_config_III_5_345_Segment_{}.XML".format(Segment_idx))
        
        filename = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\Utilities\\Reference_simulations\\Galin_Numerical_config_III_5_345_Segment_{}\\Galin_Numerical_config_III_5_345_Segment_{}.xml".format(Segment_idx, Segment_idx)
        args = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\x64\\Release\\Mjolnir_GRRT.exe -in " + filename + " -print_to_console 0"
           
        # subprocess.call(args, shell = True)
        
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

            Output_file_path = self.Sim_config.file_manager.Output_file_directory + "\\Galin_Numerical_config_III_5_345_Segment_{}\\Numerical".format(Segment_idx) 
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
            
        if not os.path.isdir(self.Sim_config.file_manager.Output_file_directory + "\\Galin_Numerical_config_III_5_345"):
            os.mkdir(self.Sim_config.file_manager.Output_file_directory + "\\Galin_Numerical_config_III_5_345")
        
        with open(self.Sim_config.file_manager.Output_file_directory + "\\Galin_Numerical_config_III_5_345\\" + self.Sim_config.simulation_name["Value"].split("_Segment")[0] + ".txt", "w") as file:
            file.write(Raw_simulation_header)
            
            if not PT_model_active:
            
                file.write("Image X Coord [M],Image Y Coord [M],Synchotron Intensity I [Jy/sRad],Synchotron Intensity Q [Jy/sRad],Synchotron Intensity U [Jy/sRad],Synchotron Intensity V [Jy/sRad],Celestial Sphere Crossing Theta [Rad],Celestial Sphere Crossing Phi [Rad]\n")
            
                for X_coord, Y_coord, I, Q, U, V, Theta, Phi in zip(X_coords, Y_coords, I_Intensity, Q_Intensity, U_Intensity, V_Intensity, Celestial_theta, Celestial_phi):
                    file.write("{},{},{},{},{},{},{}, {}\n".format(X_coord, Y_coord, I, Q, U, V, Theta, Phi))

        print("Segments exported to single file in " + self.Sim_config.file_manager.Output_file_directory + "\\Galin_Numerical_config_III_5_345")

        
if __name__ == "__main__":

    Simulation_instance = Simulation()
    
    Simulation_instance.Segment_number = 16
    
    Segment_list = []

    for Segment_idx in range(Simulation_instance.Segment_number):
        Segment_list.append(threading.Thread(target = Simulation_instance.Run_simulation, args = [Segment_idx]))

    for Segment_idx in range(Simulation_instance.Segment_number):
        Segment_list[Segment_idx].start()

    for Segment_idx in range(Simulation_instance.Segment_number):
        Segment_list[Segment_idx].join()

    Simulation_instance.Combine_segmented_output_files()


    