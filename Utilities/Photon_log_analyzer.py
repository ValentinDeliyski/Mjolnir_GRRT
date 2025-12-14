from Support_functions.Parsers import Units_class, Simulation_Parser
from Mjolnir_Configurator import Simulation_configurator
import matplotlib.pyplot as plt

from numpy import pi, array, flip, cos, sin, sqrt, tan, linspace, meshgrid
from numpy.linalg import norm
from numpy.typing import NDArray

import subprocess
import os 

class Simulation_runner():
    
    def __init__(self):
        
        """ These parameters correspond to the ones in table 1 of https://arxiv.org/pdf/2206.12066. """
        
        self.Sim_config = Simulation_configurator()
        
        self.Sim_config.metric_parameters.Numerical_metric_spline_path = "C:/Users/Valur/Documents/Repos/Mjolnir_GRRT/Utilities/Galin_numerical_config_III.XML"

        self.Sim_config.simulation_mode = {"Value": 3, "Unit": "[-]"}

        self.Sim_config.object_mass = {"Value": 6.2e9, "Unit": "[M_sun]"}

        # ================================================== Metric ================================================== #

        self.Sim_config.metric_parameters.Metric_type    = {"Value": "Numerical", "Unit": "[-]"}
        self.Sim_config.metric_parameters.Mass           = {"Value": 0.8904892552349474, "Unit": "[M]"}
        self.Sim_config.metric_parameters.Horizon_radius = {"Value": 0.01, "Unit": "[G/c^2]"}
        self.Sim_config.metric_parameters.Spin           = {"Value": 0.7813066738190858 / 0.8904892552349474, "Unit": "[M]"}
        self.Sim_config.metric_parameters.Numerical_metric_anzatz_type = {"Value": "Anzatz_1", "Unit": "[M]"} 
        
        self.Sim_config.metric_parameters.Scattering_radius = {"Value": 300, "Unit": "[M]"} 
        
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
        self.Sim_config.disk_model.Temperature_scale_factor = {"Value": 4.8e+10, "Unit": "[K]"}
                
        self.Sim_config.disk_model.Density_cutoff_radius = {"Value": 5, "Unit": "[M]"}
        self.Sim_config.disk_model.Temperature_cutoff_radius = {"Value": 5, "Unit": "[M]"}

        self.Sim_config.disk_model.Density_power_law_scale     = {"Value": 5, "Unit": "[M]"}
        self.Sim_config.disk_model.Temperature_power_law_scale = {"Value": 5, "Unit": "[M]"}

        self.Sim_config.disk_model.Opening_angle = {"Value": 0.1, "Unit": "[tan(angle)]"}
        
        self.Sim_config.disk_model.Density_power_law_power     = {"Value": 2.0, "Unit": "[-]"}
        self.Sim_config.disk_model.Temperature_power_law_power = {"Value": 1.0, "Unit": "[-]"}
        
        self.Sim_config.disk_model.Velocity_profile = {"Value": "Theta Dependant", "Unit": "[-]"}

        self.Sim_config.integrator.RK78_abs_accuracy  = {"Value": 1e-12, "Unit": "[-]"}
        self.Sim_config.integrator.RK78_rel_accuracy  = {"Value": 1e-12, "Unit": "[-]"}
        self.Sim_config.integrator.ESDIRK54_rel_accuracy  = {"Value": 1e-8, "Unit": "[-]"}
        self.Sim_config.integrator.ESDIRK54_abs_accuracy  = {"Value": 1e-8, "Unit": "[-]"}
        self.Sim_config.observer.Include_polarization = {"Value": 0, "Unit": "[-]"}
        
        # self.Sim_config.integrator.Step_controller_type  = {"Value": "PID", "Unit": "[-]"}
        self.Sim_config.integrator.max_integration_count = {"Value": 1000000, "Unit": "[-]"}
        self.Sim_config.integrator.max_affine_parameter  = {"Value": 10000000, "Unit": "[-]"}
        
        self.Sim_config.metric_parameters.Distance_to_singular_point = {"Value": 1e-3, "Unit": "[M]"}
        
        self.Sim_config.integrator.Max_rel_step_increase = {"Value": 10, "Unit": "[-]"}
        self.Sim_config.integrator.Min_rel_step_increase = {"Value": 0.1, "Unit": "[-]"}
        self.Sim_config.integrator.init_stepsize         = {"Value": 1, "Unit": "[-]"}
        
        self.Sim_config.integrator.max_stepsize         = {"Value": 10, "Unit": "[-]"}
        
        self.Sim_config.integrator.radiative_transfer_integrator_type = {"Value": "Implicit Trapezoid", "Unit": "[-]"}
        # ================================================== Hotspot ================================================== #

        self.Sim_config.hotspot_model.Density_scale_factor = {"Value": 0, "Unit": "[g / cm^3]"}
        
    def Run_simulation(self, X_init: float = 1, Y_init: float = 1):
        
        """ This simulation corresponds to the "low spin circ/thin" one in https://arxiv.org/pdf/2206.12066. Their results are presented on the top 4 panels
            of figure 3. They find a total flux of 546 mJy """
        
        self.Sim_config.sim_mode_3_X_init = {"Value": X_init, "Unit": "[M]"}
        self.Sim_config.sim_mode_3_Y_init = {"Value": Y_init, "Unit": "[M]"}
        
        """ Central black hole setup """
        self.Sim_config.metric_parameters.Spin = {"Value": 0.01, "Unit": "[M]"}

        """ Accretion disk setup """
        self.Sim_config.disk_model.Opening_angle = {"Value": 0.1, "Unit": "[tan(angle)]"}
        # self.Sim_config.disk_model.Density_power_law_scale     = {"Value": 1 + sqrt(1 - self.Sim_config.metric_parameters.Spin["Value"]**2), "Unit": "[M]"}
        # self.Sim_config.disk_model.Temperature_power_law_scale = {"Value": 1 + sqrt(1 - self.Sim_config.metric_parameters.Spin["Value"]**2), "Unit": "[M]"}
        
        self.Sim_config.disk_model.Density_scale_factor = {"Value": 1.5e6, "Unit": "[g/cm^3]"}
        
        """ The simulation name and input file path """
        self.Sim_config.simulation_name = {"Value": "Reference_Simulation_1", "Unit": "[-]"}
        
        """ The simulation output file path """
        self.Sim_config.file_manager.Output_file_directory = parent_directory + "Reference_simulations"

        self.Sim_config.generate_simulation_input(Path_to_input_dir = "Reference_simulations\\Reference_Simulation_1",
                                                               Input_file_name = "Reference_Simulation_1_input.XML")
        
        """ Run the simulation """
        filename = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\Utilities\\Reference_simulations\\Reference_Simulation_1\\Reference_Simulation_1_input.xml"
        args = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\x64\\Release\\Mjolnir_GRRT.exe -in " + filename + " -print_to_console 1"
        subprocess.call(args, shell = True)
        
        self.Sim_parser = Simulation_Parser(parent_directory + "Reference_simulations\\Reference_Simulation_1\\Numerical_photon_log")

    def Convert_spherical_to_cartesian(self, Position_tuple: tuple, BH_grid: bool = False) -> tuple[NDArray, NDArray, NDArray]:
        
        x: list[float] = []
        y: list[float] = []
        z: list[float] = []
        
        if not BH_grid:
            
            _, r_list, theta_list, phi_list = Position_tuple
        
            for r, theta, phi in zip(r_list, theta_list, phi_list):
                
                x.append(r * sin(theta) * cos(phi))
                y.append(r * sin(theta) * sin(phi))
                z.append(r * cos(theta))
                
        else:
            
            theta, phi = Position_tuple
            Event_horizon_radius = 0.01
        
            x = Event_horizon_radius * sin(phi) * cos(theta)
            y = Event_horizon_radius * sin(phi) * sin(theta)
            z = Event_horizon_radius * cos(phi)
            
        return array(x), array(y), array(z)
        
    def Plot_photon_trajectory(self, R_cutoff: float = 50) -> None:
        
        Position_tuple, _, _, _, _, _ = self.Sim_parser.get_photon_log()
        
        """ ------------ The ray trajectory ------------ """
        _, r_list, _, _ = Position_tuple
        x_coord, y_coord, z_coord = self.Convert_spherical_to_cartesian(Position_tuple = Position_tuple)

        """ ------------ Central black hole visualization ------------ """

        BH_theta_grid, BH_phi_grid = meshgrid(linspace(0, 2 * pi, 100), linspace(0, pi, 50))
        x_coord_BH, y_coord_BH, z_coord_BH = self.Convert_spherical_to_cartesian(Position_tuple = (BH_theta_grid, BH_phi_grid), BH_grid = True)

        Position_figure = plt.figure()
        Position_subplot = Position_figure.add_subplot(projection = '3d')
        Position_subplot.plot(x_coord[array(r_list) < R_cutoff], y_coord[array(r_list) < R_cutoff], z_coord[array(r_list) < R_cutoff])
        
        Position_subplot.plot_surface(x_coord_BH, y_coord_BH, z_coord_BH, color = "k", alpha = 1) # type: ignore
        
        Position_subplot.set_xlim((-R_cutoff, R_cutoff))
        Position_subplot.set_ylim((-R_cutoff, R_cutoff))
        Position_subplot.set_zlim((-R_cutoff, R_cutoff)) # type: ignore
        
        Position_subplot.set_xlabel("X [M]")
        Position_subplot.set_ylabel("Y [M]")
        Position_subplot.set_zlabel("Z [M]") # type: ignore
        Position_subplot.set_aspect('equal', 'box')
        
    def Plot_debug_data(self):
        
        _, _, _, integration_step, affine_param, Debug_tuple = self.Sim_parser.get_photon_log()
        
        Error_state, Rejected_Steps = Debug_tuple
        
        Debug_figure = plt.figure()
        
        Error_state_subplot = Debug_figure.add_subplot(131)
        Error_state_subplot.plot(affine_param, Error_state)
        
        Error_state_subplot.set_xlabel(r"$\lambda$ [M]")
        Error_state_subplot.set_ylabel(r"Normed Error State [-]")
        
        Rejected_Steps_subplot = Debug_figure.add_subplot(132)
        Rejected_Steps_subplot.plot(affine_param, Rejected_Steps)
        
        Rejected_Steps_subplot.set_xlabel(r"$\lambda$ [M]")
        Rejected_Steps_subplot.set_ylabel(r"Number of rejected steps [-]")
        
        
        Stepsize_subplot = Debug_figure.add_subplot(133)
        Stepsize_subplot.plot(affine_param, integration_step)
        
        Stepsize_subplot.set_xlabel(r"$\lambda$ [M]")
        Stepsize_subplot.set_ylabel(r"Stepsize [M]")
        
        
        
    
if __name__ == "__main__":
    
    """ ======================================= Setup ======================================= """
    
    params = {"ytick.color" : "black",
              "xtick.color" : "black",
              "axes.labelcolor" : "black",
              "axes.edgecolor" : "black",
              "text.usetex" : True,
              "font.family" : "serif",
              "font.serif" : ["Computer Modern Serif"]}
    
    plt.rcParams.update(params)
    parent_directory = os.path.abspath('...')
    
    Runner_instance = Simulation_runner()
    
    """ ================================= Run the simulation ================================= """
    
    Runner_instance.Run_simulation(X_init = 2.3196, Y_init = 4.24745)
 
    """ ====================================== Plotting ====================================== """

    Runner_instance.Plot_photon_trajectory(R_cutoff = 2)
    Runner_instance.Plot_debug_data()
    
    plt.show()