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
from numpy import pi, std, array, average
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
    
class Target_fluxes:
    
    """ Below are the average results from the 7 codes in the reference 
        https://iopscience.iop.org/article/10.3847/1538-4357/ab96c6 in units of [Jy]"""
    
    Test_1 = (1.6466 + 1.6869 + 1.6604 + 1.6466 + 1.6617 + 1.6609 + 1.6694) / 7
    Test_2 = (1.4361 + 1.4727 + 1.4486 + 1.4361 + 1.4710 + 1.4486 + 1.4568) / 7
    Test_3 = (0.4419 + 0.4527 + 0.4454 + 0.4419 + 0.4508 + 0.4456 + 0.4480) / 7
    Test_4 = (0.2711 + 0.2771 + 0.2729 + 0.2709 + 0.2763 + 0.2729 + 0.2749) / 7
    Test_5 = array([0.0256, 0.0261, 0.0258, 0.0254, 0.0260, 0.0258, 0.0259])
    Test_5_avg = average(Test_5)
    Test_5_std = std(Test_5)
    
class Phenomenological_syhnchrotron_reference_sims:
    
    def __init__(self) -> None:
        
        """ These parameters correspond to the ones in table 1 of the reference 
            https://iopscience.iop.org/article/10.3847/1538-4357/ab96c6. """
        
        self.Units = Units_class()
        self.Simulation_configurator = Simulation_configurator()
        
        """ Central black hole setup"""
        self.Simulation_configurator.metric_parameters.Metric_type = {"Value": "Kerr", "Unit": "[-]"}
        self.Simulation_configurator.object_mass                   = {"Value": 6e11 / 100 / self.Units.GR_MASS_TO_METER / self.Units.M_SUN_SI,  "Unit": "[M_sun]"}
        
        self.Object_distance = {"Value": 2.4e22 / 100 / self.Units.PC_TO_METER, "Unit": "[Pc]"}
        self.Object_Geometrical_distance = self.Object_distance["Value"] * self.Units.PC_TO_METER / self.Units.GR_MASS_TO_METER / (self.Simulation_configurator.object_mass["Value"] * self.Units.M_SUN_SI)
        
        """ Accretion disk setup """
            
        self.Simulation_configurator.disk_model.Ensamble_type = {"Value": "Phenomenological", "Unit": "[-]"}
        self.Simulation_configurator.disk_model.Disk_Model    = {"Value": "Colab_test_1",  "Unit": "[-]"}
        self.Simulation_configurator.emission_models.Emission_coeff = {"Value": 3e-18, "Unit": "[erg / (cm^3 s sr Hz)]"}
        
        self.Simulation_configurator.disk_model.Velocity_profile = {"Value": "Theta Dependant", "Unit": "[-]"}
        
        """ Observer setup """
        self.Simulation_configurator.observer.Distance    = {"Value": 1e3,           "Unit": "[M]"}
        self.Simulation_configurator.observer.Inclination = {"Value": 60 * pi / 180, "Unit": "[Rad]"}
        self.Simulation_configurator.observer.Azimuth     = {"Value": 0,             "Unit": "[Rad]"}
        
        self.Simulation_configurator.observer.Obs_frequency = {"Value": 230e9, "Unit": "[Hz]"}

        self.Simulation_configurator.observer.Image_y_min = {"Value": -15, "Unit": "[M]"}
        self.Simulation_configurator.observer.Image_y_max = {"Value":  15, "Unit": "[M]"}
        self.Simulation_configurator.observer.Image_x_min = {"Value": -15, "Unit": "[M]"}
        self.Simulation_configurator.observer.Image_x_max = {"Value":  15, "Unit": "[M]"}
        
        self.Simulation_configurator.observer.Resolution_x = {"Value": 256, "Unit": "[-]"}
        self.Simulation_configurator.observer.Resolution_y = {"Value": 256, "Unit": "[-]"}
        
        """ Kill the hotspot """
        self.Simulation_configurator.hotspot_model.Density_scale_factor = {"Value": 0, "Unit": "[g/cm^3]"}
    
        self.Simulation_configurator.observer.Include_polarization = {"Value": 0, "Unit": "[-]"}
        
        self.Simulation_configurator.geodesic_integrator.Integrator_type = {"Value": "RK78_Fehlberg", "Unit": "[-]"}
        self.Simulation_configurator.emission_integrator.Rad_Transfer_Integrator_type = {"Value": "RK78_Fehlberg", "Unit": "[-]"}
        self.Simulation_configurator.geodesic_integrator.RK_abs_accuracy = {"Value": 1e-12, "Unit": "[-]"}
        self.Simulation_configurator.geodesic_integrator.RK_rel_accuracy = {"Value": 1e-12, "Unit": "[-]"}
            
        self.Simulation_configurator.geodesic_integrator.max_upper_stepsize = {"Value": 100, "Unit": "[-]"}
        self.Simulation_configurator.geodesic_integrator.Max_rel_step_increase = {"Value": 2, "Unit": "[-]"}

    def get_total_flux(self, sim_path: str) -> float:
        
        return Simulation_Parser(sim_path + "\\Kerr").get_total_flux(self.Object_Geometrical_distance, unit = "Jy")

    def Run_and_eval_test_sim_2(self) -> None:
        
        """ Central black hole setup """
        self.Simulation_configurator.metric_parameters.Spin = {"Value": 0, "Unit": "[M]"}

        """ Accretion disk setup """
        self.Simulation_configurator.disk_model.Colab_test_1_disk_params.Vertical_scale = {"Value": 0, "Unit": "[cos(angle)]"}
        self.Simulation_configurator.disk_model.Colab_test_1_disk_params.Radial_scale = {"Value": 10, "Unit": "[M]"}
        
        self.Simulation_configurator.emission_models.Absorbtion_coeff   = {"Value": 0, "Unit": "[?]"}
        self.Simulation_configurator.emission_models.Emission_power_law = {"Value": -2, "Unit": "[-]"}

        """ The simulation name """
        self.Simulation_configurator.simulation_name = {"Value": "Phenomenological_Reference_Simulation_2", "Unit": "[-]"}
        
        """ Generate the simulation input/output file paths """
        # This is the "common" output directory, and each simulation will be in its own sub-folder
        self.Simulation_configurator.file_manager.Output_file_directory = parent_directory + "Reference_simulations"

        self.Simulation_configurator.generate_simulation_input(Path_to_input_dir = self.Simulation_configurator.file_manager.Output_file_directory + "\\" + self.Simulation_configurator.simulation_name["Value"],
                                                               Input_file_name = self.Simulation_configurator.simulation_name["Value"] + "_input.XML")
        
        # This is the path to the simulation sub-folder with the output files from the ray tracer and the input XML
        Sim_subfolder = self.Simulation_configurator.file_manager.Output_file_directory + "\\" + self.Simulation_configurator.simulation_name["Value"]
        
        # This is the full path to the input XML
        Input_file_path = Sim_subfolder + "\\" + self.Simulation_configurator.simulation_name["Value"] + "_input.XML"
        
        """ Run the simulation """
        Comman_line_args = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\x64\\Release\\Mjolnir_GRRT.exe -in " + Input_file_path + " -print_to_console 0"
        subprocess.call(Comman_line_args, shell = True)
                
        """ Evaluate the simulataion results """
        Total_flux = self.get_total_flux(Sim_subfolder)
        Relative_error = (Total_flux - Target_fluxes.Test_2) / Target_fluxes.Test_2
        
        try:
            assert(abs(Relative_error) < 0.01)
            print(f"{bcolors.BOLD}" + f"{bcolors.OKGREEN}Reference Simulation 2 pass with a relative error of {{}} %.{bcolors.ENDC}".format(round(Relative_error * 100, 3)))
        except:
            print(f"{bcolors.BOLD}" + f"{bcolors.FAIL}Reference Simulation 2 fail with a relative error of {{}} %.{bcolors.ENDC}".format(round(Relative_error * 100, 3)))

    def Run_and_eval_test_sim_3(self) -> None:
        
        """ Central black hole setup """
        self.Simulation_configurator.metric_parameters.Spin = {"Value": 0.9, "Unit": "[M]"}

        """ Accretion disk setup """
        self.Simulation_configurator.disk_model.Colab_test_1_disk_params.Vertical_scale  = {"Value": 10 / 3, "Unit": "[tan(angle)]"}
        self.Simulation_configurator.disk_model.Colab_test_1_disk_params.Radial_scale  = {"Value": 10, "Unit": "[M]"}
        
        self.Simulation_configurator.emission_models.Absorbtion_coeff   = {"Value": 0, "Unit": "[?]"}
        self.Simulation_configurator.emission_models.Emission_power_law = {"Value": 0, "Unit": "[-]"}

        """ The simulation name """
        self.Simulation_configurator.simulation_name = {"Value": "Phenomenological_Reference_Simulation_3", "Unit": "[-]"}
        
        """ Generate the simulation input/output file paths """
        # This is the "common" output directory, and each simulation will be in its own sub-folder
        self.Simulation_configurator.file_manager.Output_file_directory = parent_directory + "Reference_simulations"

        self.Simulation_configurator.generate_simulation_input(Path_to_input_dir = self.Simulation_configurator.file_manager.Output_file_directory + "\\" + self.Simulation_configurator.simulation_name["Value"],
                                                               Input_file_name = self.Simulation_configurator.simulation_name["Value"] + "_input.XML")
        
        # This is the path to the simulation sub-folder with the output files from the ray tracer and the input XML
        Sim_subfolder = self.Simulation_configurator.file_manager.Output_file_directory + "\\" + self.Simulation_configurator.simulation_name["Value"]
        
        # This is the full path to the input XML
        Input_file_path = Sim_subfolder + "\\" + self.Simulation_configurator.simulation_name["Value"] + "_input.XML"
        
        """ Run the simulation """
        Comman_line_args = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\x64\\Release\\Mjolnir_GRRT.exe -in " + Input_file_path + " -print_to_console 0"
        subprocess.call(Comman_line_args, shell = True)
                
        """ Evaluate the simulataion results """
        Total_flux = self.get_total_flux(Sim_subfolder)
        Relative_error = (Total_flux - Target_fluxes.Test_3) / Target_fluxes.Test_3
        
        try:
            assert(abs(Relative_error) < 0.01)
            print(f"{bcolors.BOLD}" + f"{bcolors.OKGREEN}Reference Simulation 3 pass with a relative error of {{}} %.{bcolors.ENDC}".format(round(Relative_error * 100, 3)))
        except:
            print(f"{bcolors.BOLD}" + f"{bcolors.FAIL}Reference Simulation 3 fail with a relative error of {{}} %.{bcolors.ENDC}".format(round(Relative_error * 100, 3)))
        
    def Run_and_eval_test_sim_4(self) -> None:
        
        """ Central black hole setup """
        self.Simulation_configurator.metric_parameters.Spin = {"Value": 0.9, "Unit": "[M]"}

        """ Accretion disk setup """
        self.Simulation_configurator.disk_model.Colab_test_1_disk_params.Vertical_scale  = {"Value": 10 / 3, "Unit": "[tan(angle)]"}
        self.Simulation_configurator.disk_model.Colab_test_1_disk_params.Radial_scale  = {"Value": 10, "Unit": "[M]"}
        
        self.Simulation_configurator.emission_models.Absorbtion_coeff   = {"Value": 1e5, "Unit": "[?]"}
        self.Simulation_configurator.emission_models.Emission_power_law = {"Value": 0, "Unit": "[-]"}

        """ The simulation name """
        self.Simulation_configurator.simulation_name = {"Value": "Phenomenological_Reference_Simulation_4", "Unit": "[-]"}
        
        """ Generate the simulation input/output file paths """
        # This is the "common" output directory, and each simulation will be in its own sub-folder
        self.Simulation_configurator.file_manager.Output_file_directory = parent_directory + "Reference_simulations"

        self.Simulation_configurator.generate_simulation_input(Path_to_input_dir = self.Simulation_configurator.file_manager.Output_file_directory + "\\" + self.Simulation_configurator.simulation_name["Value"],
                                                               Input_file_name = self.Simulation_configurator.simulation_name["Value"] + "_input.XML")
        
        # This is the path to the simulation sub-folder with the output files from the ray tracer and the input XML
        Sim_subfolder = self.Simulation_configurator.file_manager.Output_file_directory + "\\" + self.Simulation_configurator.simulation_name["Value"]
        
        # This is the full path to the input XML
        Input_file_path = Sim_subfolder + "\\" + self.Simulation_configurator.simulation_name["Value"] + "_input.XML"

        """ Run the simulation """
        Comman_line_args = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\x64\\Release\\Mjolnir_GRRT.exe -in " + Input_file_path + " -print_to_console 0"
        subprocess.call(Comman_line_args, shell = True)
                
        """ Evaluate the simulataion results """
        Total_flux = self.get_total_flux(Sim_subfolder)
        Relative_error = (Total_flux - Target_fluxes.Test_4) / Target_fluxes.Test_4
        
        try:
            assert(abs(Relative_error) < 0.01)
            print(f"{bcolors.BOLD}" + f"{bcolors.OKGREEN}Reference Simulation 4 pass with a relative error of {{}} %.{bcolors.ENDC}".format(round(Relative_error * 100, 3)))
        except:
            print(f"{bcolors.BOLD}" + f"{bcolors.FAIL}Reference Simulation 4 fail with a relative error of {{}} %.{bcolors.ENDC}".format(round(Relative_error * 100, 3)))
            
    def Run_and_eval_test_sim_5(self) -> None:
        
        self.Simulation_configurator.geodesic_integrator.min_upper_stepsize = {"Value": 0.15, "Unit": "[-]"}
        self.Simulation_configurator.geodesic_integrator.max_step_b_coeff = {"Value": 0.004, "Unit": "[-]"}
        self.Simulation_configurator.geodesic_integrator.dist_at_min_upper_stepsize = {"Value": 25, "Unit": "[-]"}
        
        """ Central black hole setup """
        self.Simulation_configurator.metric_parameters.Spin = {"Value": 0.0, "Unit": "[M]"}

        """ Accretion disk setup """
        self.Simulation_configurator.disk_model.Colab_test_1_disk_params.Vertical_scale  = {"Value": 100 / 3, "Unit": "[tan(angle)]"}
        self.Simulation_configurator.disk_model.Colab_test_1_disk_params.Radial_scale  = {"Value": 10, "Unit": "[M]"}
        
        self.Simulation_configurator.emission_models.Absorbtion_coeff   = {"Value": 1e6, "Unit": "[?]"}
        self.Simulation_configurator.emission_models.Emission_power_law = {"Value": 0, "Unit": "[-]"}

        """ The simulation name """
        self.Simulation_configurator.simulation_name = {"Value": "Phenomenological_Reference_Simulation_5", "Unit": "[-]"}
        
        """ Generate the simulation input/output file paths """
        # This is the "common" output directory, and each simulation will be in its own sub-folder
        self.Simulation_configurator.file_manager.Output_file_directory = parent_directory + "Reference_simulations"

        self.Simulation_configurator.generate_simulation_input(Path_to_input_dir = self.Simulation_configurator.file_manager.Output_file_directory + "\\" + self.Simulation_configurator.simulation_name["Value"],
                                                               Input_file_name = self.Simulation_configurator.simulation_name["Value"] + "_input.XML")
        
        # This is the path to the simulation sub-folder with the output files from the ray tracer and the input XML
        Sim_subfolder = self.Simulation_configurator.file_manager.Output_file_directory + "\\" + self.Simulation_configurator.simulation_name["Value"]
        
        # This is the full path to the input XML
        Input_file_path = Sim_subfolder + "\\" + self.Simulation_configurator.simulation_name["Value"] + "_input.XML"

        """ Run the simulation """
        Comman_line_args = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\x64\\Release\\Mjolnir_GRRT.exe -in " + Input_file_path + " -print_to_console 0"
        subprocess.call(Comman_line_args, shell = True)
                
        """ Evaluate the simulataion results """
        Total_flux = self.get_total_flux(Sim_subfolder)
        Error = (Total_flux - Target_fluxes.Test_5_avg)
        
        try:
            assert(Target_fluxes.Test_5.min() < Total_flux < Target_fluxes.Test_5.max())
            print(f"{bcolors.BOLD}" + f"{bcolors.OKGREEN}Reference Simulation 5 pass with a total flux of {{}} [Jy] and relative error of {{}} standard deviations from the average.{bcolors.ENDC}".format(round(Total_flux, 5), 
                                                                                                                                                                                                           round(abs(Error / Target_fluxes.Test_5_std), 3)))
        except:
            print(f"{bcolors.BOLD}" + f"{bcolors.FAIL}Reference Simulation 5 fail with a total flux of {{}} [Jy] and relative error of {{}} standard deviations from the average.{bcolors.ENDC}".format(round(Total_flux, 5), 
                                                                                                                                                                                                        round(abs(Error / Target_fluxes.Test_5_std), 3)))

if __name__ == "__main__": 
  
    Phenomenological_syhnchrotron_reference_sims_instance = Phenomenological_syhnchrotron_reference_sims()

    Sim_2_thread = threading.Thread(target = Phenomenological_syhnchrotron_reference_sims_instance.Run_and_eval_test_sim_2, args = [])
    Sim_3_thread = threading.Thread(target = Phenomenological_syhnchrotron_reference_sims_instance.Run_and_eval_test_sim_3, args = [])
    Sim_4_thread = threading.Thread(target = Phenomenological_syhnchrotron_reference_sims_instance.Run_and_eval_test_sim_4, args = [])
    Sim_5_thread = threading.Thread(target = Phenomenological_syhnchrotron_reference_sims_instance.Run_and_eval_test_sim_5, args = [])

    """ Run the simulation threads - the sleep calls inbetween are so the input file has time to actually generate. """

    # Sim_2_thread.start()
    time.sleep(1)
    # Sim_3_thread.start()
    time.sleep(1)
    # Sim_4_thread.start()
    time.sleep(1)
    Sim_5_thread.start()