from Support_functions.Parsers import Units_class, Simulation_Parser
from Mjolnir_Configurator import Simulation_configurator
import matplotlib.pyplot as plt

from numpy import pi, dot, float64, cosh, sinh, exp, array, flip, cos, sin, ones, arctan
from numpy.linalg import norm
from numpy.typing import NDArray

import os 

def Configure_and_make_Mjolnir_photon_log():
    
    Sim_config = Simulation_configurator()

    Sim_config.simulation_mode = {"Value": 2, "Unit": "[-]"}
    
    """ Magic number equal to 1 / the conversion factor from geometric to CGS stepsizes. """
    Sim_config.object_mass = {"Value": 1 / 147706.32775277024, "Unit": "[M_sun]"}

    # ================================================== Metric ================================================== #

    Sim_config.metric_parameters.Metric_type = {"Value": "Kerr", "Unit": "[-]"}
    Sim_config.metric_parameters.Spin        = {"Value": 0.999, "Unit": "[M]"}
    
    # ================================================== Observer ================================================== #

    Sim_config.observer.Distance    = {"Value": 1e4, "Unit": "[M]"}
    Sim_config.observer.Inclination = {"Value": 80 * pi / 180, "Unit": "[Rad]"}
    Sim_config.observer.Obs_frequency = {"Value": 230e9, "Unit": "[Hz]"}
    Sim_config.observer.Include_polarization = {"Value": 0, "Unit": "[-]"}
    Sim_config.metric_parameters.Scattering_radius = {"Value": 1000, "Unit": "[-]"}

    # ================================================== Disk ================================================== #
    
    Sim_config.disk_model.Disk_Model = {"Value": "Novikov-Thorne", "Unit": "[-]"}   
    Sim_config.disk_model.r_in_NT_disk = {"Value": 1.2, "Unit": "[M]"} 
    Sim_config.disk_model.r_out_NT_disk = {"Value": 500, "Unit": "[M]"} 
    
    Sim_config.min_image_order = {"Value": 1, "Unit": "[M]"} 
    Sim_config.sim_mode_3_X_init = {"Value": 0.01, "Unit": "[M]"} 
    Sim_config.sim_mode_3_Y_init = {"Value": -8, "Unit": "[M]"} 
        
    # Sim_config.disk_model.Ensamble_type = {"Value": "Debug_constant_functions", "Unit": "[-]"}
    # Sim_config.disk_model.Disk_Model    = {"Value": "Debug_constant_density", "Unit": "[-]"}
    
    # ================================================== Emission models ================================================== #
    
    Sim_config.emission_models.Debug_j_I_value = {"Value": 0, "Unit": "[-]"}
    Sim_config.emission_models.Debug_j_Q_value = {"Value": 0.1, "Unit": "[-]"}
    Sim_config.emission_models.Debug_j_U_value = {"Value": 0.1, "Unit": "[-]"}
    Sim_config.emission_models.Debug_j_V_value = {"Value": 0.1, "Unit": "[-]"}
    
    Sim_config.emission_models.Debug_alpha_I_value = {"Value": 0, "Unit": "[-]"}
    Sim_config.emission_models.Debug_alpha_Q_value = {"Value": 0, "Unit": "[-]"}
    Sim_config.emission_models.Debug_alpha_U_value = {"Value": 0, "Unit": "[-]"}
    Sim_config.emission_models.Debug_alpha_V_value = {"Value": 0, "Unit": "[-]"}
    
    Sim_config.emission_models.Debug_rho_I_value = {"Value": 0, "Unit": "[-]"}
    Sim_config.emission_models.Debug_rho_Q_value = {"Value": 10, "Unit": "[-]"}
    Sim_config.emission_models.Debug_rho_U_value = {"Value": 0, "Unit": "[-]"}
    Sim_config.emission_models.Debug_rho_V_value = {"Value": -4, "Unit": "[-]"}
    
    Sim_config.disk_model.Mag_field_geometry_r     = {"Value": 0, "Unit": "[-]"}
    Sim_config.disk_model.Mag_field_geometry_theta = {"Value": 1, "Unit": "[-]"}
    Sim_config.disk_model.Mag_field_geometry_phi   = {"Value": 1, "Unit": "[-]"}
        
    # ================================================== Integrator ================================================== #
    
    Sim_config.geodesic_integrator.RK_abs_accuracy = {"Value": 1e-14, "Unit": "[-]"}
    Sim_config.geodesic_integrator.RK_rel_accuracy = {"Value": 1e-14, "Unit": "[-]"}
    Sim_config.geodesic_integrator.max_affine_parameter = {"Value": 5e4, "Unit": "[M]"}
    Sim_config.geodesic_integrator.use_adaptive_step    = {"Value": 1, "Unit": "[-]"}
    Sim_config.geodesic_integrator.max_stepsize  = {"Value": 100, "Unit": "[-]"}
     
    # ================================================== Hotspot ================================================== #

    Sim_config.hotspot_model.Density_scale_factor = {"Value": 0, "Unit": "[g / cm^3]"}
    
    """ The simulation name and input file path """
    Sim_config.simulation_name = {"Value": "Plasma_integration_test_1", "Unit": "[-]"}

    """ The simulation output file path """
    Sim_config.file_manager.Output_file_directory = parent_directory + "Reference_simulations"
    Sim_config.simulation_name = {"Value": "Integration_tests", "Unit": "[-]"}

    Sim_config.generate_simulation_input(Path_to_input_dir = "Reference_simulations\\Integration_tests",
                                         Input_file_name = "Plasma_integration_test_1.XML")

    import subprocess

    filename = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\Utilities\\Reference_simulations\\Integration_tests\\Plasma_integration_test_1.xml"
    args = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\x64\\Release\\Mjolnir_GRRT.exe -in " + filename + " -print_to_console 1"
    
    subprocess.call(args, shell = True)

def Compute_analytic_radiative_transfer(Affine_parameter_log: NDArray, Emission_functions: list[float], Absorbtions_functions: list[float], Faradey_functions: list[float], 
                                        I_init: float = 0, Q_init: float = 0, U_init: float = 0, V_init: float = 0):
    
    Aff_param = Affine_parameter_log
    
    Init_Polarization_State: list[float] = [I_init, Q_init, U_init, V_init]
    
    j_I: float = Emission_functions[0]
    j_Q: float = Emission_functions[1]
    j_U: float = Emission_functions[2]
    j_V: float = Emission_functions[3]
    
    alpha_I: float = Absorbtions_functions[0]
    alpha_Q: float = Absorbtions_functions[1]
    alpha_U: float = Absorbtions_functions[2]
    alpha_V: float = Absorbtions_functions[3]
    alpha_p: float64 = norm(Absorbtions_functions[1:])
    
    rho_I: float = Faradey_functions[0]
    rho_Q: float = Faradey_functions[1]
    rho_U: float = Faradey_functions[2]
    rho_V: float = Faradey_functions[3]
    rho: float64 = norm(Faradey_functions[1:])

    if alpha_p != 0:

        I = (I_init * cosh(alpha_p * Aff_param) - dot(Init_Polarization_State[1:], Absorbtions_functions[1:]) / alpha_p * sinh(alpha_p * Aff_param)) * exp(-alpha_I * Aff_param) 
        I = I + dot(Emission_functions[1:], Absorbtions_functions[1:]) / (alpha_I**2 - alpha_p**2) * (-1 + (alpha_I * sinh(alpha_p * Aff_param) + alpha_p * cosh(alpha_p * Aff_param)) / alpha_p * exp(-alpha_I * Aff_param))
        I = I + alpha_I * j_I / (alpha_I**2 - alpha_p**2) * (1 - (alpha_I * cosh(alpha_p * Aff_param) + alpha_p * sinh(alpha_p * Aff_param)) / alpha_I * exp(-alpha_I * Aff_param))
        
        Q = (Q_init + alpha_Q * dot(Init_Polarization_State[1:], Absorbtions_functions[1:]) / alpha_p**2 * (cosh(alpha_p * Aff_param) - 1) - I_init * alpha_Q / alpha_p * sinh(alpha_p * Aff_param)) * exp(-alpha_I * Aff_param) 
        Q = Q + j_Q * (1 - exp(-alpha_I * Aff_param)) / alpha_I
        Q = Q + (dot(Emission_functions[1:], Absorbtions_functions[1:]) * alpha_Q / alpha_I / (alpha_I**2 - alpha_p**2) * 
                (1 - (1 - alpha_I**2 / alpha_Q**2 + alpha_I / alpha_Q**2 * (alpha_I * cosh(alpha_p * Aff_param) + alpha_p * sinh(alpha_p * Aff_param)))* exp(-alpha_I * Aff_param)))
        Q = Q + j_I * alpha_Q / alpha_p / (alpha_I**2 - alpha_p**2) * (-alpha_p + (alpha_p * cosh(alpha_p * Aff_param) + alpha_I * sinh(alpha_p * Aff_param)) * exp(-alpha_I * Aff_param))
        
        U = (U_init + alpha_U * dot(Init_Polarization_State[1:], Absorbtions_functions[1:]) / alpha_p**2 * (cosh(alpha_p * Aff_param) - 1) - I_init * alpha_U / alpha_p * sinh(alpha_p * Aff_param)) * exp(-alpha_I * Aff_param) 
        U = U + j_U * (1 - exp(-alpha_I * Aff_param)) / alpha_I
        U = U + (dot(Emission_functions[1:], Absorbtions_functions[1:]) * alpha_U / alpha_I / (alpha_I**2 - alpha_p**2) * 
                (1 - (1 - alpha_I**2 / alpha_Q**2 + alpha_I / alpha_Q**2 * (alpha_I * cosh(alpha_p * Aff_param) + alpha_p * sinh(alpha_p * Aff_param)))* exp(-alpha_I * Aff_param)))
        U = U + j_I * alpha_U / alpha_p / (alpha_I**2 - alpha_p**2) * (-alpha_p + (alpha_p * cosh(alpha_p * Aff_param) + alpha_I * sinh(alpha_p * Aff_param)) * exp(-alpha_I * Aff_param))
        
        V = (V_init + alpha_V * dot(Init_Polarization_State[1:], Absorbtions_functions[1:]) / alpha_p**2 * (cosh(alpha_p * Aff_param) - 1) - I_init * alpha_U / alpha_p * sinh(alpha_p * Aff_param)) * exp(-alpha_I * Aff_param) 
        V = V + j_V * (1 - exp(-alpha_I * Aff_param)) / alpha_I
        V = V + (dot(Emission_functions[1:], Absorbtions_functions[1:]) * alpha_V / alpha_I / (alpha_I**2 - alpha_p**2) * 
                (1 - (1 - alpha_I**2 / alpha_Q**2 + alpha_I / alpha_Q**2 * (alpha_I * cosh(alpha_p * Aff_param) + alpha_p * sinh(alpha_p * Aff_param)))* exp(-alpha_I * Aff_param)))
        V = V + j_I * alpha_V / alpha_p / (alpha_I**2 - alpha_p**2) * (-alpha_p + (alpha_p * cosh(alpha_p * Aff_param) + alpha_I * sinh(alpha_p * Aff_param)) * exp(-alpha_I * Aff_param))
        
    else:
        
        I = I_init * ones(len(Aff_param))
        
        Q = rho_Q / rho**2 * (j_Q * rho_Q + j_V * rho_V) * Aff_param - rho_V / rho**3 * (j_V * rho_Q - j_Q * rho_V) * sin(rho * Aff_param) - j_U * rho_V / rho**2 * (1 - cos(rho * Aff_param))
        
        U = (j_Q * rho_V - j_V * rho_Q) / rho**2 * (1 - cos(rho * Aff_param)) + j_U / rho * sin(rho * Aff_param)
        
        V = rho_V / rho**2 * (j_Q * rho_Q + j_V * rho_V) * Aff_param - rho_Q / rho**3 * (j_Q * rho_V - j_V * rho_Q) * sin(rho * Aff_param) + j_U * rho_Q / rho**2 * (1 - cos(rho * Aff_param))
            
    return I, Q, U, V

if __name__ == "__main__":
    
    params = {"ytick.color" : "black",
              "xtick.color" : "black",
              "axes.labelcolor" : "black",
              "axes.edgecolor" : "black",
              "text.usetex" : True,
              "font.family" : "serif",
              "font.serif" : ["Computer Modern Serif"]}
    
    plt.rcParams.update(params)
    
    parent_directory = os.path.abspath('...')
    
    Configure_and_make_Mjolnir_photon_log()
    
    Sim_parser = Simulation_Parser(parent_directory + "Reference_simulations\\Integration_tests\\Kerr_photon_log")
    
    Position_tuple, Momentum_tuple, Emission_tuple, Integration_step, Affine_param, Debug_tuple, Polarization_tuple = Sim_parser.get_photon_log()
    
    Affine_param = (array(Affine_param))
    
    # ===================================================================== Coordinate plotting =====================================================================
    
    t_coord, r_coord, theta_coord, phi_coord = Position_tuple
    _, (t_subplot, r_subplot, theta_subplot, phi_subplot) = plt.subplots(1, 4, gridspec_kw = {'width_ratios': [1, 1, 1, 1]}, constrained_layout = True)
        
    t_subplot.plot(Affine_param, t_coord)
    t_subplot.set_xlabel(r"Affine Parameter [M]")
    t_subplot.set_ylabel(r"$t$ coordinate [M]")
    
    r_subplot.plot(Affine_param, r_coord)
    r_subplot.set_xlabel(r"Affine Parameter [M]")
    r_subplot.set_ylabel(r"$r$ coordinate [M]")
    
    theta_subplot.plot(Affine_param, theta_coord)
    theta_subplot.set_xlabel(r"Affine Parameter [M]")
    theta_subplot.set_ylabel(r"$\theta$ coordinate [rad]")
    
    phi_subplot.plot(Affine_param, phi_coord)
    phi_subplot.set_xlabel(r"Affine Parameter [M]")
    phi_subplot.set_ylabel(r"$\phi$ coordinate [rad]")
    
    # ===================================================================== Momentum plotting =====================================================================
   
    p_t, p_r, p_theta, p_phi = Momentum_tuple
    _, (p_t_subplot, p_r_subplot, p_theta_subplot, p_phi_subplot) = plt.subplots(1, 4, gridspec_kw = {'width_ratios': [1, 1, 1, 1]}, constrained_layout = True)
        
    p_t_subplot.plot(Affine_param, p_t)
    p_t_subplot.set_xlabel(r"Affine Parameter [M]")
    p_t_subplot.set_ylabel(r"$p_t$ [-]")
    
    p_r_subplot.plot(Affine_param, p_r)
    p_r_subplot.set_xlabel(r"Affine Parameter [M]")
    p_r_subplot.set_ylabel(r"$p_r$ [-]")
    
    p_theta_subplot.plot(Affine_param, p_theta)
    p_theta_subplot.set_xlabel(r"Affine Parameter [M]")
    p_theta_subplot.set_ylabel(r"$p_\theta$ [rad/M]")
    
    p_phi_subplot.plot(Affine_param, p_phi)
    p_phi_subplot.set_xlabel(r"Affine Parameter [M]")
    p_phi_subplot.set_ylabel(r"$p_\phi$ [rad/M]")
    
    # ==================================================================  Polarization plotting ===================================================================
    
    Pol_x, Pol_y = Polarization_tuple
    _, (EVPA_subplot) = plt.subplots(1, 1, gridspec_kw = {'width_ratios': [1]}, constrained_layout = True)
    
    EVPA_subplot.plot(Affine_param, arctan(-array(Pol_x) / array(Pol_y)))
    
    # ===================================================================== Emission plotting =====================================================================
   
    I, Q, U, V = Emission_tuple
    I, Q, U, V = flip(I), flip(Q), flip(U), flip(V)
    _, (I_subplot, Q_subplot, U_subplot, V_subplot) = plt.subplots(1, 4, gridspec_kw = {'width_ratios': [1, 1, 1, 1]}, constrained_layout = True)
        
    I_subplot.plot(Affine_param, I, "r")
    I_subplot.set_xlabel(r"Affine Parameter [M]")
    I_subplot.set_ylabel(r"Stokes I [Jy/sRad]")
    
    Q_subplot.plot(Affine_param, Q)
    Q_subplot.set_xlabel(r"Affine Parameter [M]")
    Q_subplot.set_ylabel(r"Stokes Q [Jy/sRad]")
    
    U_subplot.plot(Affine_param, U)
    U_subplot.set_xlabel(r"Affine Parameter [M]")
    U_subplot.set_ylabel(r"Stokes U [Jy/sRad]")
    
    V_subplot.plot(Affine_param, V)
    V_subplot.set_xlabel(r"Affine Parameter [M]")
    V_subplot.set_ylabel(r"Stokes V [Jy/sRad]")
    
    # ===================================================================== Emission vs analytic emission plotting =====================================================================
   
    I_analytic, Q_analytic, U_analytic, V_analytic = Compute_analytic_radiative_transfer(Affine_parameter_log = Affine_param, Emission_functions = [0, 0.1, 0.1, 0.1], Absorbtions_functions = [0, 0, 0, 0], Faradey_functions = [0, 10, 0, -4])
    
    _, Subplots = plt.subplots(2, 4, gridspec_kw = {'width_ratios': [1, 1, 1, 1]}, constrained_layout = True)
    
    Subplots[0][0].plot(Affine_param, I, "r")
    Subplots[0][0].plot(Affine_param, I_analytic)
    Subplots[0][0].set_xlabel(r"Affine Parameter [M]")
    Subplots[0][0].set_ylabel(r"Stokes I [Jy/sRad]")
    
    Subplots[1][0].plot(Affine_param, I - I_analytic)
    Subplots[1][0].set_xlabel(r"Affine Parameter [M]")
    Subplots[1][0].set_ylabel(r"Stokes I Delta [Jy/sRad]")
    
    Subplots[0][1].plot(Affine_param, Q, "r")
    Subplots[0][1].plot(Affine_param, Q_analytic)
    Subplots[0][1].set_xlabel(r"Affine Parameter [M]")
    Subplots[0][1].set_ylabel(r"Stokes Q [Jy/sRad]")
    
    Subplots[1][1].plot(Affine_param, Q - Q_analytic)
    Subplots[1][1].set_xlabel(r"Affine Parameter [M]")
    Subplots[1][1].set_ylabel(r"Stokes Q Delta [Jy/sRad]")
    
    Subplots[0][2].plot(Affine_param, U, "r")
    Subplots[0][2].plot(Affine_param, U_analytic)
    Subplots[0][2].set_xlabel(r"Affine Parameter [M]")
    Subplots[0][2].set_ylabel(r"Stokes U [Jy/sRad]")
    
    Subplots[1][2].plot(Affine_param, U - U_analytic)
    Subplots[1][2].set_xlabel(r"Affine Parameter [M]")
    Subplots[1][2].set_ylabel(r"Stokes U Delta [Jy/sRad]")
    
    Subplots[0][3].plot(Affine_param, V, "r")
    Subplots[0][3].plot(Affine_param, V_analytic)
    Subplots[0][3].set_xlabel(r"Affine Parameter [M]")
    Subplots[0][3].set_ylabel(r"Stokes V [Jy/sRad]")
    
    Subplots[1][3].plot(Affine_param, V - V_analytic)
    Subplots[1][3].set_xlabel(r"Affine Parameter [M]")
    Subplots[1][3].set_ylabel(r"Stokes V Delta [Jy/sRad]")
    
 
    plt.show()