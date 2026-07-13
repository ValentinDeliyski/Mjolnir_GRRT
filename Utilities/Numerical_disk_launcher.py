import sys
import os

""" Add the parent directory of this file to the search path, 
    so this script can be ran from the "Utilities" folder """
    
parent_directory = os.path.abspath('...')
sys.path.append(parent_directory)

from Mjolnir_Configurator import Simulation_configurator
from Support_functions.Parsers import Units_class, Simulation_Parser
from Support_functions.Spacetimes_new import Kerr

from numpy import pi, tan, sqrt, linspace, array, float64, append, arctan
from numpy.typing import NDArray

from time import sleep

import matplotlib.pyplot as plt
import subprocess
import multiprocessing
import shutil

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
    
    def __init__(self, resolution, inclination):
        
        self.Output_folder = "C:/Users/Valur/Documents/Repos/Mjolnir_GRRT/Utilities/Numerical_disk_runs/M87_mass_{}_deg/".format(inclination)
        
        self.Segment_number = 16
        
        self.Nominal_Resolution_x = resolution
        self.Nominal_Resolution_y = resolution
        
        self.Units = Units_class()
        self.Simulation_configurator = Simulation_configurator()
        
        self.Simulation_configurator.simulation_mode = {"Value": 0, "Unit": "[-]"} 
        self.Simulation_configurator.object_mass = {"Value": 6.2e9, "Unit": "[M_sun]"} 
        self.Simulation_configurator.average_emission_pitch_angle = {"Value": 0, "Unit": "[-]"}
        
        self.Simulation_configurator.observer.Cam_rotation_angle = {"Value": 0, "Unit": "[-]"}
        self.Simulation_configurator.observer.Use_angular_coords = {"Value": 1, "Unit": "[-]"}
        
        self.Simulation_configurator.disk_model.Threshold_relative_density = {"Value": 1e-5, "Unit": "[-]"}
        
        self.Simulation_configurator.observer.Include_polarization = {"Value": 1, "Unit": "[-]"}
        
        """ Central black hole setup"""
        self.Simulation_configurator.metric_parameters.Metric_type = {"Value": "Kerr",  "Unit": "[-]"}
        self.Simulation_configurator.metric_parameters.Spin        = {"Value": 0.5, "Unit": "[M]"}
        
        """ Accretion disk setup """   
        self.Simulation_configurator.disk_model.Disk_Model = {"Value": "Numerical", "Unit": "[-]"}
        self.Simulation_configurator.disk_model.Ensamble_type = {"Value": "Thermal", "Unit": "[-]"}
        
        self.Simulation_configurator.disk_model.Numerical_disk_params.Spline_type = {"Value": "GSL_linear", "Unit": "[-]"}
        
        self.Simulation_configurator.disk_model.r_ISCO = {"Value": 4.233, "Unit": "[M]"}
        self.Simulation_configurator.disk_model.Ang_momentum_below_ISCO = {"Value": 3.414213760169089, "Unit": "[-]"}
        self.Simulation_configurator.disk_model.Ang_momentum_exponent = {"Value": 0.5, "Unit": "[-]"}
        
        self.Simulation_configurator.disk_model.Velocity_profile = {"Value": "von_Zeipel_cylinder", "Unit": "[-]"}
        
        self.Simulation_configurator.disk_model.Numerical_disk_params.Numerical_XML_path = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\Utilities\\Test.XML"
        
        """ Kill the hotspot """
        self.Simulation_configurator.hotspot_model.Enabled_flag = {"Value": 0, "Unit": "[-]"}
        
        """ Observer setup """
        self.Simulation_configurator.observer.Distance    = {"Value": 1e4,         "Unit":  "[M]" }
        self.Simulation_configurator.observer.Inclination = {"Value": inclination * pi / 180, "Unit": "[Rad]"}
        self.Simulation_configurator.observer.Azimuth     = {"Value": 0,           "Unit": "[Rad]"}
        
        self.Simulation_configurator.observer.Cam_rotation_angle = {"Value": 0, "Unit": "[Hz]"}
        
        self.Observer_FOV = {"Value": 150, "Unit": "[micro-arcsec]"}
    
        self.Nominal_Image_x_max = 15
        self.Nominal_Image_y_max = 15
    
        self.Simulation_configurator.observer.Image_y_min = {"Value": -15, "Unit": "[M]"}
        self.Simulation_configurator.observer.Image_y_max = {"Value":  15, "Unit": "[M]"}
        self.Simulation_configurator.observer.Image_x_min = {"Value": -15, "Unit": "[M]"}
        self.Simulation_configurator.observer.Image_x_max = {"Value":  15, "Unit": "[M]"}
        
        self.Nominal_Image_x_angle_max = arctan(self.Nominal_Image_x_max / self.Simulation_configurator.observer.Distance["Value"])
        self.Nominal_Image_y_angle_max = arctan(self.Nominal_Image_y_max / self.Simulation_configurator.observer.Distance["Value"])
         
        self.Simulation_configurator.observer.Resolution_x = {"Value": resolution, "Unit": "[-]"}
        self.Simulation_configurator.observer.Resolution_y = {"Value": resolution, "Unit": "[-]"}
        
        """ Configure the integrator """
        
        self.Simulation_configurator.geodesic_integrator.Integrator_type = {"Value": "RK78_Fehlberg", "Unit": "[-]"}
        self.Simulation_configurator.emission_integrator.Rad_Transfer_Integrator_type = {"Value": "Analytic", "Unit": "[-]"}       
         
        self.Simulation_configurator.geodesic_integrator.min_upper_stepsize = {"Value": 1, "Unit": "[-]"}
        self.Simulation_configurator.geodesic_integrator.max_step_b_coeff = {"Value": 0.0001, "Unit": "[-]"}
        self.Simulation_configurator.geodesic_integrator.Max_step_in_emission_medium = {"Value": 0.1, "Unit": "[-]"}
        self.Simulation_configurator.geodesic_integrator.dist_at_min_upper_stepsize = {"Value": 25, "Unit": "[-]"}
        
        self.Simulation_configurator.geodesic_integrator.max_integration_count = {"Value": 10000000, "Unit": "[-]"}
        self.Simulation_configurator.geodesic_integrator.RK_abs_accuracy = {"Value": 1e-12, "Unit": "[-]"}
        self.Simulation_configurator.geodesic_integrator.RK_rel_accuracy = {"Value": 1e-12, "Unit": "[-]"}
        
        self.Simulation_configurator.geodesic_integrator.max_upper_stepsize = {"Value": 100, "Unit": "[-]"}
        self.Simulation_configurator.geodesic_integrator.Max_rel_step_increase = {"Value": 2, "Unit": "[-]"}

    def run_simulation(self):
        
        self.Simulation_configurator.generate_simulation_input(Path_to_input_dir = "Numerical_disks\\" + self.Simulation_configurator.simulation_name["Value"],
                                                               Input_file_name   = self.Simulation_configurator.simulation_name["Value"] + "_input.XML")
            
        """ Run the simulation """
        filename = ("C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\Utilities\\" + 
                    "Numerical_disks\\" + 
                    self.Simulation_configurator.simulation_name["Value"] + 
                    "\\" + 
                    self.Simulation_configurator.simulation_name["Value"] + 
                    "_input.XML")    
        
        args = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\x64\\Release\\Mjolnir_GRRT.exe -in " + filename + " -print_to_console 0"
        
        subprocess.call(args, shell = True)          

    def run_numerical_metric_segment(self, Segment_idx):

        glock.acquire()

        """ The segmented observation window """
        Nominal_y_scan_step = 2 * self.Nominal_Image_y_angle_max / (self.Nominal_Resolution_y - 1)
        
        self.Simulation_configurator.observer.Image_y_angle_max = {"Value": -self.Nominal_Image_y_angle_max + (Segment_idx + 1) * (2 * self.Nominal_Image_y_angle_max / self.Segment_number), "Unit": "[M]"}
        self.Simulation_configurator.observer.Image_y_angle_min = {"Value": -self.Nominal_Image_y_angle_max + (Segment_idx + 0) * (2 * self.Nominal_Image_y_angle_max / self.Segment_number), "Unit": "[M]"}

        if Segment_idx > 0 and Segment_idx < self.Segment_number - 1:
            self.Simulation_configurator.observer.Image_y_angle_max["Value"] -= Nominal_y_scan_step / 2
            self.Simulation_configurator.observer.Image_y_angle_min["Value"] += Nominal_y_scan_step / 2

        elif Segment_idx == self.Segment_number - 1:
            self.Simulation_configurator.observer.Image_y_angle_min["Value"] += Nominal_y_scan_step / 2

        elif Segment_idx == 0:
            self.Simulation_configurator.observer.Image_y_angle_max["Value"] -= Nominal_y_scan_step / 2

        self.Simulation_configurator.observer.Image_x_angle_max = {"Value":  self.Nominal_Image_x_angle_max, "Unit": "[M]"}
        self.Simulation_configurator.observer.Image_x_angle_min = {"Value": -self.Nominal_Image_x_angle_max, "Unit": "[M]"}

        self.Simulation_configurator.observer.Resolution_x = {"Value": self.Nominal_Resolution_x, "Unit": "[-]"}
        self.Simulation_configurator.observer.Resolution_y = {"Value": int(self.Nominal_Resolution_y / self.Segment_number), "Unit": "[-]"}
        
        """ The simulation output file path """
        self.Simulation_configurator.file_manager.Output_file_directory = self.Output_folder
        self.Simulation_configurator.simulation_name = {"Value": "Numerical_Segment_{}".format(Segment_idx), "Unit": "[-]"}
        
        self.Simulation_configurator.generate_simulation_input(Path_to_input_dir = self.Output_folder + self.Simulation_configurator.simulation_name["Value"],
                                                            Input_file_name = "Numerical_Segment_{}.XML".format(Segment_idx))
        
        filename =  self.Output_folder + self.Simulation_configurator.simulation_name["Value"] + "/Numerical_Segment_{}.XML".format(Segment_idx)
        args = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\x64\\Release\\Mjolnir_GRRT.exe -in " + filename + " -print_to_console 0"
        
        sleep(0.1)
        glock.release()
           
        subprocess.run(args, shell = False)
        
    def Combine_segmented_output_files(self):
            
            X_coords: NDArray[float64]      = array([])
            Y_coords: NDArray[float64]      = array([])
            I_Intensity: NDArray[float64]   = array([])
            Q_Intensity: NDArray[float64]   = array([])
            U_Intensity: NDArray[float64]   = array([])
            V_Intensity: NDArray[float64]   = array([])
            Optical_Depth: NDArray[float64]   = array([])
            Faraday_Q_Depth: NDArray[float64]   = array([])
            Faraday_V_Depth: NDArray[float64]   = array([])
            Disk_redshift: NDArray[float64] = array([])
            Disk_flux: NDArray[float64]     = array([])
            Final_t_coord : NDArray[float64]     = array([])
            Celestial_theta: NDArray[float64] = array([])
            Celestial_phi: NDArray[float64]   = array([])

            Source_t: NDArray[float64] = array([])
            Source_r: NDArray[float64] = array([])
            Source_phi: NDArray[float64] = array([])
            Source_p_r: NDArray[float64] = array([])
            Source_p_theta: NDArray[float64] = array([])
            Source_p_phi: NDArray[float64] = array([])
            
            Polarization_vec_X: NDArray[float64] = array([])
            Polarization_vec_Y: NDArray[float64] = array([])

            Raw_simulation_header: str = ""
            NT_model_active: bool = False
            
            for Segment_idx in range(self.Segment_number):

                Output_file_path = self.Output_folder + "\\Numerical_Segment_{}\\Kerr".format(Segment_idx) 
                Sim_parser = Simulation_Parser(Output_file_path)
            
                NT_model_active = Sim_parser.Simulation_metadata["Active disk model"] == "Novikov-Thorne"
                
                X_coords = append(X_coords, Sim_parser.X_coords)
                Y_coords = append(Y_coords, Sim_parser.Y_coords)
                I_Intensity = append(I_Intensity, Sim_parser.I_Intensity)
                Q_Intensity = append(Q_Intensity, Sim_parser.Q_Intensity)
                U_Intensity = append(U_Intensity, Sim_parser.U_Intensity)
                V_Intensity = append(V_Intensity, Sim_parser.V_Intensity)
                Faraday_Q_Depth = append(Faraday_Q_Depth, Sim_parser.Faraday_Q_Depth)
                Faraday_V_Depth = append(Faraday_V_Depth, Sim_parser.Faraday_V_Depth)
                Optical_Depth = append(Optical_Depth, Sim_parser.Optical_Depth)
                Disk_redshift = append(Disk_redshift, Sim_parser.Disk_redshift)
                Disk_flux = append(Disk_flux, Sim_parser.Disk_flux)
                Final_t_coord = append(Final_t_coord, Sim_parser.Final_t_coord)
                
                Source_t = append(Source_t, Sim_parser.Source_t)
                Source_r = append(Source_r, Sim_parser.Source_r)
                Source_phi = append(Source_phi, Sim_parser.Source_phi)
                Source_p_r = append(Source_p_r, Sim_parser.Source_p_r)
                Source_p_theta = append(Source_p_theta, Sim_parser.Source_p_theta)
                Source_p_phi = append(Source_p_phi, Sim_parser.Source_p_phi)

                Polarization_vec_X = append(Polarization_vec_X, Sim_parser.Polarization_vec_X)
                Polarization_vec_Y = append(Polarization_vec_Y, Sim_parser.Polarization_vec_Y)
                
                Celestial_theta = append(Celestial_theta, Sim_parser.Celestial_theta)
                Celestial_phi = append(Celestial_phi, Sim_parser.Celestial_phi)
                
                Raw_simulation_header_split = Sim_parser.Raw_simulation_header.split("\n")
                Raw_simulation_header_split[13] = Raw_simulation_header_split[13].replace("{} x {}".format(self.Nominal_Resolution_x, int(self.Nominal_Resolution_y / self.Segment_number)), 
                                                                                          "{} x {}".format(self.Nominal_Resolution_x, self.Nominal_Resolution_y))
                
                Raw_simulation_header_split[12] = Raw_simulation_header_split[12].replace(Sim_parser.Simulation_metadata["Observation Window Dimentions (-X,+X,-Y,+Y) [M]"], 
                                                                                         "{},{},{},{}".format(-self.Nominal_Image_x_max, self.Nominal_Image_x_max, -self.Nominal_Image_y_max, self.Nominal_Image_y_max))
                
                Raw_simulation_header = '\n'.join([x for x in Raw_simulation_header_split if not x == ''])
                
            """ ============================ Flatten all arrays, so I can write in the combined file easily ============================ """
                
            X_coords = X_coords.flatten()
            Y_coords = Y_coords.flatten()
            I_Intensity = I_Intensity.flatten()
            Q_Intensity = Q_Intensity.flatten()
            U_Intensity = U_Intensity.flatten()
            V_Intensity = V_Intensity.flatten()
            Optical_Depth = Optical_Depth.flatten()
            Faraday_Q_Depth = Faraday_Q_Depth.flatten()
            Faraday_V_Depth = Faraday_V_Depth.flatten()
            Disk_redshift = Disk_redshift.flatten()
            Disk_flux = Disk_flux.flatten()
            Celestial_theta = Celestial_theta.flatten()
            Celestial_phi = Celestial_phi.flatten()
            
            Source_t = Source_t.flatten()
            Source_r = Source_r.flatten()
            Source_phi = Source_phi.flatten()
            Source_p_r = Source_p_r.flatten()
            Source_p_theta = Source_p_theta.flatten()
            Source_p_phi = Source_p_phi.flatten()

            Polarization_vec_X = Polarization_vec_X.flatten()
            Polarization_vec_Y = Polarization_vec_Y.flatten()
                         
            if not os.path.isdir(self.Output_folder + "\\Numerical_results"):
                os.mkdir(self.Output_folder + "\\Numerical_results")
            
            with open(self.Output_folder + "\\Numerical_results\\All_Segments_Results.txt", "w") as file:
                file.write(Raw_simulation_header)
                
                if not NT_model_active:
                
                    file.write("\nImage X Coord [M],Image Y Coord [M],Synchotron Intensity I [Jy/sRad],Synchotron Intensity Q [Jy/sRad],Synchotron Intensity U [Jy/sRad],Synchotron Intensity V [Jy/sRad],Final t Coordinate [M],Total Optical Depth [-],Total Faraday Q Depth [-],Total Faraday V Depth [-],Celestial Sphere Crossing Theta [Rad],Celestial Sphere Crossing Phi [Rad],\n")
                
                    for Tuple_to_write in zip(X_coords, Y_coords, I_Intensity, Q_Intensity, U_Intensity, V_Intensity, Final_t_coord, Optical_Depth, Faraday_Q_Depth, Faraday_V_Depth, Celestial_theta, Celestial_phi):
                        file.write("{},{},{},{},{},{},{}, {}, {}, {}, {}, {}\n".format(*Tuple_to_write))

                else:
                    
                    file.write("\nImage X Coord [M],Image Y Coord [M],Disk Redshift [-],Disk Flux [M_dot/M^2],Polarization vector X [-],Polarization vector Y [-],Source t Coord [M],Source r Coord [M],Source phi Coord [Rad],Radial Momentum (covariant),Theta Momentum (covariant),Phi Momentum (covariant),Celestial Sphere Crossing Theta [Rad],Celestial Sphere Crossing Phi [Rad],\n")
                
                    for Tuple_to_write in zip(X_coords, Y_coords, Disk_redshift, Disk_flux, Polarization_vec_X, Polarization_vec_Y, Source_t, 
                                              Source_r, Source_phi, Source_p_r, Source_p_theta, Source_p_phi, Celestial_theta, Celestial_phi):
                        file.write("{},{},{},{},{},{},{},{},{},{},{},{},{},{},\n".format(*Tuple_to_write))

            print("Segments exported to single file in " + self.Output_folder + "\\Numerical_results\\All_Segments_Results.txt")
          
            # for Segment_idx in range(self.Segment_number):
                
            #     shutil.rmtree(self.Output_folder + "\\Numerical_Segment_{}".format(Segment_idx))
                
            # print("Temporary folders and results files deleted.")
            
def init_global_lock(lock):
    global glock
    glock = lock
    
if __name__ == "__main__":
    
    for inc in [20, 70]:
        
        Sim_instance = Simulation(inclination = inc, resolution = 256)
        
        # Sim_instance.run_numerical_metric_segment(6)
        
        if 1: #not os.path.isfile(Sim_instance.Output_folder + "Numerical_results\\All_Segments_Results.txt"):
                    
            Process_lock = multiprocessing.Lock()
            Process_pool = multiprocessing.Pool(processes = Sim_instance.Segment_number, initializer = init_global_lock, initargs = (Process_lock,))
            Process_pool.map(Sim_instance.run_numerical_metric_segment, range(Sim_instance.Segment_number))
            Process_pool.close()
            Process_pool.join()

            Sim_instance.Combine_segmented_output_files()
    