import sys
import os

from Support_functions.Parsers import Units_class
import xml.etree.cElementTree as ET
import xml.dom.minidom
import os
from numpy import pi, sqrt

""" Add the parent directory of this file to the search path, 
    so this script can be ran from the "Utilities" folder """
    
parent_directory = os.path.abspath('...')
sys.path.append(parent_directory)

class Integrator():

    __slots__ = ("init_stepsize",
                 "RK78_accuracy", 
                 "Step_controller_type",
                 "step_controller_safety_factor_1",
                 "step_controller_safety_factor_2",
                 "PID_controller_I_gain",
                 "PID_controller_P_gain",
                 "PID_controller_D_gain",
                 "Max_rel_step_increase",
                 "Min_rel_step_increase",
                 "Gustafsson_controller_k_1", 
                 "Gustafsson_controller_k_2", 
                 "max_integration_count",
                 "simpson_method_accuracy",
                 "max_affine_parameter",
                 "use_adaptive_step",
                 "max_stepsize",
                 "radiative_transfer_integrator_type")

class Disk_model():

    __slots__ = ("Ensamble_type",
                 "Disk_Model",
                 "Velocity_profile",
                 "Radial_velocity_fraction",
                 "Density_scale_factor",
                 "Temperature_scale_factor",
                 
                 "Opening_angle",
                 "Density_power_law_scale",
                 "Density_power_law_power",
                 "Temperature_power_law_scale",
                 "Temperature_power_law_power",
                 "Density_cutoff_radius",
                 "Density_cutoff_scale",
                 "Temperature_cutoff_radius",
                 "Temperature_cutoff_scale",
 
                 "Radial_scale",
                 "Vertical_scale",
                 
                 "Magnetization",
                 "Mag_field_geometry_r",
                 "Mag_field_geometry_theta",
                 "Mag_field_geometry_phi",
                 "Mag_field_magnitude_scale",
                 "Mag_field_power",
                 "Mag_field_radial_scale",
                 "Mag_field_geometry",
                 "Mag_field_magnitude_profile",

                 "Threshold_relative_density",
                 
                 "r_in_PT_disk",
                 "r_out_PT_disk")
    
class Hotspot_model():

    __slots__ = ("Ensamble_type",
                 "Density_profile",
                 "Temperature_profile",
                 "Velocity_profile",
                 "Radial_velocity_fraction",
                 
                 "Density_scale_factor",
                 "Density_spread",
                 "Density_power_law_power",
                 "Density_power_law_scale",
                 
                 "Temperature_scale_factor",
                 "Temperature_spread",
                 "Temperature_power_law_power",
                 "Temperature_power_law_scale",
                 
                 "Radius",
                 
                 "Mag_field_geometry_r",
                 "Mag_field_geometry_theta",
                 "Mag_field_geometry_phi",
                 "Mag_field_magnitude_scale",
                 "Mag_field_power",
                 "Mag_field_radial_scale",
                 "Mag_field_geometry",
                 "Mag_field_magnitude_profile",
                 
                 "Temporal_spread",
                 "Coord_time_at_max",
                 "Distance",
                 "Inclination",
                 "Azimuth",
                 "Magnetization",
                 "Threshold_relative_density")

class Metric_parameters():

    __slots__ = ("Mass",
                 "Spin",
                 "Horizon_radius",
                 "WH_redshift",
                 "WH_r_throat", 
                 "WH_stop_at_throat",
                 "RBH_param", 
                 "JNW_gamma", 
                 "EGB_gamma", 
                 "Halo_compactness", 
                 "Halo_mass",
                 "Metric_type", 
                 "Numerical_metric_spline_path",
                 "Numerical_metric_anzatz_type",
                 "Distance_to_singular_point",
                 "Scattering_radius")

class Observer():

    __slots__ = ("Init_time",
                 "Distance",
                 "Inclination",
                 "Azimuth",
                 "Cam_rotation_angle",
                 "Image_y_min",
                 "Image_y_max",
                 "Image_x_min",
                 "Image_x_max",
                 "Resolution_y",
                 "Resolution_x",
                 "Include_polarization",
                 "Obs_frequency")

class Emission_models():

    __slots__ = ("Emission_power_law",
                 "Source_f_power_law",
                 "Absorbtion_coeff",
                 "Emission_coeff",
                 "Kappa",
                 "Debug_j_I_value",
                 "Debug_j_Q_value",
                 "Debug_j_U_value",
                 "Debug_j_V_value",
                 "Debug_alpha_I_value",
                 "Debug_alpha_Q_value",
                 "Debug_alpha_U_value",
                 "Debug_alpha_V_value",
                 "Debug_rho_I_value",
                 "Debug_rho_Q_value",
                 "Debug_rho_U_value",
                 "Debug_rho_V_value",)  

class File_manager():

    __slots__ = ("Sim_mode_2_input_file_path", 
                 "Output_file_directory", 
                 "Common_file_names", 
                 "Vert_shader_path", 
                 "Frag_shader_path", 
                 "Truncate_files")

class Simulation_configurator:

    __slots__ = ("simulation_name", 
                 "integrator", 
                 "metric_parameters", 
                 "disk_model", 
                 "hotspot_model", 
                 "observer", 
                 "emission_models", 
                 "NT_model_params", 
                 "file_manager", 
                 "average_emission_pitch_angle", 
                 "thermalize_emission_medium",
                 "emission_pitch_angle_samples_to_average",
                 "object_mass",
                 "simulation_mode",
                 "sim_mode_2_param_value_number",
                 "sim_mode_3_X_init",
                 "sim_mode_3_Y_init",
                 "max_image_order")

    def __init__(self, 
                 Average_emission_pitch_angle: dict[str, int | str] = {"Value": 1, "Unit": "[-]"}, 
                 thermalize_emission_medium: dict[str, int | str] = {"Value": 0, "Unit": "[-]"}, 
                 emission_pitch_angle_samples_to_average: dict[str, int | str] = {"Value": 50, "Units": "[-]"},
                 object_mass: dict[str, float | str] = {"Value": 6.2e9, "Unit": "[M_sun]"},
                 simulation_name: dict[str, str] = {"Value": "Test_Simulation", "Unit": "[-]"},
                 simulation_mode: dict[str, int | str] = {"Value": 1, "Unit": "[-]"}, 
                 sim_mode_2_param_value_number: dict[str, int | str] = {"Value": 1, "Unit": "[-]"},
                 sim_mode_3_X_init: dict[str, float | str] = {"Value": 1, "Unit": "[M]"},
                 sim_mode_3_Y_init: dict[str, float | str] = {"Value": 1, "Unit": "[M]"},
                 max_image_order: dict[str, int | str] = {"Value": 3, "Unit": "[-]"},):

        self.average_emission_pitch_angle = Average_emission_pitch_angle
        self.thermalize_emission_medium = thermalize_emission_medium
        self.emission_pitch_angle_samples_to_average = emission_pitch_angle_samples_to_average
        self.simulation_name = simulation_name
        self.object_mass = object_mass
        self.simulation_mode = simulation_mode
        self.sim_mode_2_param_value_number = sim_mode_2_param_value_number
        self.sim_mode_3_X_init = sim_mode_3_X_init
        self.sim_mode_3_Y_init = sim_mode_3_Y_init
        self.max_image_order = max_image_order

        self._configure_integrator_settings()
        self._configure_observer()
        self._configure_disk_model()
        self._configure_hotspot_model()
        self._configure_file_manager()
        self._configure_metric_parameters()
        self._configure_emission_models()

    def _configure_integrator_settings(self, Init_stepsize: dict[str, float | str] = {"Value": 1e-5, "Unit": "[M]"},
                                             RK78_accuracy: dict[str, float | str] = {"Value": 1e-13, "Unit": "[-]"},
                                             Step_controller_type: dict[str, str] = {"Value": "Gustafsson", "Unit": "[-]"},
                                             Safety_factor_1: dict[str, float | str] = {"Value": 0.9, "Unit": "[-]"},
                                             Safety_factor_2: dict[str, float | str] = {"Value": 1e-35, "Unit": "[-]"},
                                             Max_rel_step_increase: dict[str, float | str] = {"Value": 2, "Unit": "[-]"},
                                             Min_rel_step_increase: dict[str, float | str] = {"Value": 0.01, "Unit": "[-]"},
                                             Step_controller_I_gain: dict[str, float | str] = {"Value": -0.58 / 7, "Unit": "[-]"},
                                             Step_controller_P_gain: dict[str, float | str] = {"Value": 0.21 / 7, "Unit": "[-]"},
                                             Step_controller_D_gain: dict[str, float | str] = {"Value": -0.1 / 7, "Unit": "[-]"},
                                             Gustafsson_controller_k_1: dict[str, float | str] = {"Value": -0.367 / 7, "Unit": "[-]"},
                                             Gustafsson_controller_k_2: dict[str, float | str] = {"Value": 0.268 / 7, "Unit": "[-]"},
                                             Max_integration_count: dict[str, float | str] = {"Value": 1e7, "Unit": "[-]"},
                                             simpson_method_accuracy: dict[str, float | str] = {"Value": 1e-6, "Unit": "[-]"},
                                             max_affine_parameter: dict[str, float | str] = {"Value": 1e6, "Unit": "[M]"},
                                             use_adaptive_step: dict[str, int | str] = {"Value": 1, "Unit": "[M]"},
                                             max_stepsize: dict[str, int | str] = {"Value": 5, "Unit": "[M]"},
                                             radiative_transfer_integrator_type: dict[str, str] = {"Value": "Implicit Trapezoid", "Unit": "[-]"}):

        self.integrator = Integrator()

        self.integrator.init_stepsize = Init_stepsize
        self.integrator.init_stepsize = Init_stepsize
        self.integrator.RK78_accuracy = RK78_accuracy
        self.integrator.Step_controller_type = Step_controller_type
        self.integrator.step_controller_safety_factor_1 = Safety_factor_1
        self.integrator.step_controller_safety_factor_2 = Safety_factor_2
        self.integrator.Max_rel_step_increase = Max_rel_step_increase
        self.integrator.Min_rel_step_increase = Min_rel_step_increase
        self.integrator.PID_controller_I_gain = Step_controller_I_gain
        self.integrator.PID_controller_P_gain = Step_controller_P_gain
        self.integrator.PID_controller_D_gain = Step_controller_D_gain
        self.integrator.Gustafsson_controller_k_1 = Gustafsson_controller_k_1
        self.integrator.Gustafsson_controller_k_2 = Gustafsson_controller_k_2
        self.integrator.max_integration_count  = Max_integration_count
        self.integrator.simpson_method_accuracy = simpson_method_accuracy
        self.integrator.max_affine_parameter = max_affine_parameter
        self.integrator.use_adaptive_step = use_adaptive_step
        self.integrator.max_stepsize = max_stepsize
        self.integrator.radiative_transfer_integrator_type = radiative_transfer_integrator_type

    def _configure_observer(self, Init_time:dict[str, float | str] = {"Value": 0, "Unit": "[M]"},
                                  Distance: dict[str, float | str] = {"Value": 1e4, "Unit": "[M]"},
                                  Inclination: dict[str, float | str] = {"Value": 160 / 180 * pi, "Unit": "[Rad]"},
                                  Azimuth: dict[str, float | str] = {"Value": 0.0, "Unit": "[Rad]"},
                                  Cam_rotation_angle: dict[str, float | str] = {"Value": 0.0, "Unit": "[Rad]"},
                                  Image_y_min: dict[str, float | str] = {"Value": -15, "Unit": "[M]"},
                                  Image_y_max: dict[str, float | str] = {"Value":  15, "Unit": "[M]"},
                                  Image_x_min: dict[str, float | str] = {"Value": -15, "Unit": "[M]"},
                                  Image_x_max: dict[str, float | str] = {"Value":  15, "Unit": "[M]"},
                                  Resolution_y: dict[str, int | str] = {"Value": 2048, "Unit": "[-]"},
                                  Resolution_x: dict[str, int | str] = {"Value": 2048, "Unit": "[-]"},
                                  Observation_frequency: dict[str, float | str] = {"Value": 230e9, "Unit": "[Hz]"},
                                  Include_polarization: dict[str, int | str] = {"Value": 0, "Unit": "[-]"}):
        
        self.observer = Observer()

        self.observer.Init_time = Init_time
        self.observer.Distance = Distance
        self.observer.Inclination = Inclination
        self.observer.Azimuth = Azimuth
        self.observer.Cam_rotation_angle = Cam_rotation_angle
        self.observer.Image_y_min  = Image_y_min
        self.observer.Image_y_max  = Image_y_max
        self.observer.Image_x_min  = Image_x_min
        self.observer.Image_x_max  = Image_x_max
        self.observer.Resolution_x = Resolution_x
        self.observer.Resolution_y = Resolution_y

        self.observer.Include_polarization = Include_polarization
        self.observer.Obs_frequency        = Observation_frequency

    def _configure_metric_parameters(self, Mass: dict[str, float | str] = {"Value": 1.0, "Unit": "[M]"}, 
                                           Spin: dict[str, float | str] = {"Value": 0.98, "Unit": "[M]"}, 
                                           Horizon_radius: dict[str, float | str] = {"Value": 2, "Unit": "[M]"}, 
                                           WH_redshift: dict[str, float | str] = {"Value": 2.0, "Unit": "[M]"},
                                           WH_r_throat: dict[str, float | str] = {"Value": 1.0, "Unit": "[M]"}, 
                                           WH_stop_at_throat: dict[str, float | str] = {"Value": 0, "Unit": "[-]"}, 
                                           RBH_param: dict[str, float | str] = {"Value": 0.5, "Unit": "[M]"}, 
                                           JNW_gamma: dict[str, float | str] = {"Value": 0.48, "Unit": "[-]"}, 
                                           EGB_gamma: dict[str, float | str] = {"Value": 1.15, "Unit": "[M^2]"}, 
                                           Halo_compactness: dict[str, float | str] = {"Value": 1e-4, "Unit": "[-]"},
                                           Halo_mass: dict[str, float | str] = {"Value": 1e4, "Unit": "[M]"},
                                           Metric_type: dict[str, float | str] = {"Value": "Kerr", "Unit": "[-]"},
                                           Numerical_metric_spline_path: str = "",
                                           Numerical_metric_anzatz_type: dict[str, str] = {"Value": "Anzatz_1", "Unit": "[-]"},
                                           Distance_to_singular_point: dict[str, float | str] = {"Value": 1e-4, "Unit": "[M]"},
                                           Scattering_radius: dict[str, float | str] = {"Value": 100, "Unit": "[M]"},):

        self.metric_parameters = Metric_parameters()

        self.metric_parameters.Mass = Mass
        self.metric_parameters.Spin = Spin
        self.metric_parameters.WH_redshift = WH_redshift
        self.metric_parameters.WH_r_throat = WH_r_throat
        self.metric_parameters.WH_stop_at_throat = WH_stop_at_throat
        self.metric_parameters.RBH_param = RBH_param
        self.metric_parameters.JNW_gamma = JNW_gamma
        self.metric_parameters.EGB_gamma = EGB_gamma
        self.metric_parameters.Halo_compactness = Halo_compactness
        self.metric_parameters.Halo_mass = Halo_mass
        self.metric_parameters.Metric_type = Metric_type
        self.metric_parameters.Numerical_metric_spline_path = Numerical_metric_spline_path
        self.metric_parameters.Numerical_metric_anzatz_type = Numerical_metric_anzatz_type
        self.metric_parameters.Scattering_radius = Scattering_radius
        self.metric_parameters.Distance_to_singular_point = Distance_to_singular_point

    def _configure_emission_models(self, Emission_power_law: dict[str, float | str] = {"Value": 0.0, "Unit": "[-]"},
                                         Source_f_power_law: dict[str, float | str] = {"Value": 2.5, "Unit": "[-]"},
                                         Absorbtion_coeff: dict[str, float | str] = {"Value": 1e5, "Unit": "[?]"},
                                         Emission_coeff: dict[str, float | str] = {"Value": 3e-18, "Unit": "[erg / (cm^3 s sr Hz)]"},
                                         Kappa: dict[str, float | str] = {"Value": 4.0, "Unit": "[-]"},
                                         Debug_j_I_value: dict[str, float | str] = {"Value": 0.0, "Unit": "[-]"},
                                         Debug_j_Q_value: dict[str, float | str] = {"Value": 0.0, "Unit": "[-]"},
                                         Debug_j_U_value: dict[str, float | str] = {"Value": 0.0, "Unit": "[-]"},
                                         Debug_j_V_value: dict[str, float | str] = {"Value": 0.0, "Unit": "[-]"},
                                         Debug_alpha_I_value: dict[str, float | str] = {"Value": 0.0, "Unit": "[-]"},
                                         Debug_alpha_Q_value: dict[str, float | str] = {"Value": 0.0, "Unit": "[-]"},
                                         Debug_alpha_U_value: dict[str, float | str] = {"Value": 0.0, "Unit": "[-]"},
                                         Debug_alpha_V_value: dict[str, float | str] = {"Value": 0.0, "Unit": "[-]"},
                                         Debug_rho_I_value: dict[str, float | str] = {"Value": 0.0, "Unit": "[-]"},
                                         Debug_rho_Q_value: dict[str, float | str] = {"Value": 0.0, "Unit": "[-]"},
                                         Debug_rho_U_value: dict[str, float | str] = {"Value": 0.0, "Unit": "[-]"},
                                         Debug_rho_V_value: dict[str, float | str] = {"Value": 0.0, "Unit": "[-]"}):
        
        self.emission_models = Emission_models()

        self.emission_models.Emission_power_law = Emission_power_law
        self.emission_models.Source_f_power_law = Source_f_power_law
        self.emission_models.Absorbtion_coeff   = Absorbtion_coeff
        self.emission_models.Emission_coeff     = Emission_coeff
        self.emission_models.Kappa              = Kappa
        
        self.emission_models.Debug_j_I_value = Debug_j_I_value
        self.emission_models.Debug_j_Q_value = Debug_j_Q_value
        self.emission_models.Debug_j_U_value = Debug_j_U_value
        self.emission_models.Debug_j_V_value = Debug_j_V_value
        self.emission_models.Debug_alpha_I_value = Debug_alpha_I_value
        self.emission_models.Debug_alpha_Q_value = Debug_alpha_Q_value
        self.emission_models.Debug_alpha_U_value = Debug_alpha_U_value
        self.emission_models.Debug_alpha_V_value = Debug_alpha_V_value
        self.emission_models.Debug_rho_I_value = Debug_rho_I_value
        self.emission_models.Debug_rho_Q_value = Debug_rho_Q_value
        self.emission_models.Debug_rho_U_value = Debug_rho_U_value
        self.emission_models.Debug_rho_V_value = Debug_rho_V_value\
            
    def _configure_disk_model(self, Ensamble_type: dict[str, str] = {"Value": "Thermal", "Unit": "[-]"},
                                    Disk_Model: dict[str, str] = {"Value": "Phenom_RIAF_1", "Unit": "[-]"},
                                    Velocity_profile: dict[str, str] = {"Value": "Theta Dependant", "Unit": "[-]"},
                                    Radial_velocity_fraction: dict[str, float | str] = {"Value": 0, "Unit": "[-]"},
                                    Density_scale_factor: dict[str, float | str] = {"Value": 1e5, "Unit": "[g/cm^3]"},
                                    Temperature_scale_factor: dict[str, float | str] = {"Value": 1e11, "Unit": "[K]"},
                                    
                                    Opening_angle: dict[str, float | str] = {"Value": 0.1, "Unit": "[tan(angle)]"},
                                    Density_power_law_scale: dict[str, float | str] = {"Value": 5.0, "Unit": "[M]"},
                                    Density_power_law_power : dict[str, float | str] = {"Value": 2.0, "Unit": "[-]"},
                                    Temperature_power_law_scale: dict[str, float | str] = {"Value": 5.0, "Unit": "[M]"},
                                    Temperature_power_law_power: dict[str, float | str] = {"Value": 1.0, "Unit": "[-]"},    
                                    Density_cutoff_radius: dict[str, float | str] = {"Value": 5.0, "Unit": "[M]"},
                                    Density_cutoff_scale: dict[str, float | str] = {"Value": 0.4, "Unit": "[M]"},
                                    Temperature_cutoff_radius: dict[str, float | str] = {"Value": 5.0, "Unit": "[M]"},
                                    Temperature_cutoff_scale: dict[str, float | str] = {"Value": 0.4, "Unit": "[M]"},
                                    
                                    Radial_scale: dict[str, float | str] = {"Value": 0.3, "Unit": "[M]"},
                                    Vertical_scale: dict[str, float | str] = {"Value": 10.0, "Unit": "[M]"},
                                    
                                    Magnetization: dict[str, float | str] = {"Value": 0.01, "Unit": "[-]"},
                                    Mag_field_geometry_r: dict[str, float | str] = {"Value": 0.5, "Unit": "[-]"},
                                    Mag_field_geometry_theta: dict[str, float | str] = {"Value": 0, "Unit": "[-]"},
                                    Mag_field_geometry_phi: dict[str, float | str] = {"Value": 0.87, "Unit": "[-]"},
                                    
                                    Mag_field_magnitude_scale: dict[str, float | str] = {"Value": 100, "Unit": "[G]"},
                                    Mag_field_power: dict[str, float | str] = {"Value": 1, "Unit": "[-]"},
                                    Mag_field_radial_scale: dict[str, float | str] = {"Value": 5, "Unit": "[M]"},
                                    Mag_field_geometry: dict[str, str] = {"Value": "Constant", "Unit": "[-]"},
                                    Mag_field_magnitude_profile: dict[str, str] = {"Value": "Magnetization_based", "Unit": "[-]"},
                                    
                                    Threshold_relative_density: dict[str, float | str] = {"Value": 1e-3, "Unit": "[-]"},
                                    
                                    r_in_PT_disk: dict[str, float | str] = {"Value": 6, "Unit": "[M]"},
                                    r_out_PT_disk: dict[str, float | str] = {"Value": 25, "Unit": "[M]"}):
        
        self.disk_model = Disk_model()

        self.disk_model.Ensamble_type  = Ensamble_type
        self.disk_model.Disk_Model     = Disk_Model
        self.disk_model.Velocity_profile         = Velocity_profile
        self.disk_model.Radial_velocity_fraction = Radial_velocity_fraction
        self.disk_model.Density_scale_factor     = Density_scale_factor
        self.disk_model.Temperature_scale_factor = Temperature_scale_factor

        self.disk_model.Mag_field_geometry_r     = Mag_field_geometry_r
        self.disk_model.Mag_field_geometry_theta = Mag_field_geometry_theta
        self.disk_model.Mag_field_geometry_phi   = Mag_field_geometry_phi
        self.disk_model.Magnetization            = Magnetization

        self.disk_model.Opening_angle = Opening_angle
        self.disk_model.Density_power_law_scale = Density_power_law_scale
        self.disk_model.Density_power_law_power = Density_power_law_power
        self.disk_model.Temperature_power_law_scale = Temperature_power_law_scale
        self.disk_model.Temperature_power_law_power = Temperature_power_law_power
        self.disk_model.Density_cutoff_radius     = Density_cutoff_radius
        self.disk_model.Density_cutoff_scale      = Density_cutoff_scale
        self.disk_model.Temperature_cutoff_radius = Temperature_cutoff_radius
        self.disk_model.Temperature_cutoff_scale  = Temperature_cutoff_scale

        self.disk_model.Radial_scale = Radial_scale
        self.disk_model.Vertical_scale = Vertical_scale

        self.disk_model.Mag_field_magnitude_scale = Mag_field_magnitude_scale
        self.disk_model.Mag_field_power        = Mag_field_power
        self.disk_model.Mag_field_radial_scale = Mag_field_radial_scale
        self.disk_model.Mag_field_geometry     = Mag_field_geometry
        self.disk_model.Mag_field_magnitude_profile  = Mag_field_magnitude_profile
        
        self.disk_model.Threshold_relative_density = Threshold_relative_density 
        
        self.disk_model.r_in_PT_disk = r_in_PT_disk
        self.disk_model.r_out_PT_disk = r_out_PT_disk

    def _configure_hotspot_model(self, Ensamble_type: dict[str, str] = {"Value": "Kappa", "Unit": "[-]"},
                                       Density_profile: dict[str, str] = {"Value": "Gaussian", "Unit": "[-]"},
                                       Temperature_profile: dict[str, str] = {"Value": "Gaussian", "Unit": "[-]"},
                                       Velocity_profile: dict[str, str] = {"Value": "Theta Dependant", "Unit": "[-]"},
                                       Radial_velocity_fraction: dict[str, float | str] = {"Value": 0, "Unit": "[-]"},
                                       
                                       Radius: dict[str, float | str] = {"Value": 1, "Unit": "[M]"},
                                       Temporal_spread: dict[str, float | str] = {"Value": 85, "Unit": "[GM/c^3]"},
                                       Density_spread: dict[str, float | str] = {"Value": 1.0, "Unit": "[M]"},
                                       Density_power_law_power: dict[str, float | str] = {"Value": 0.0, "Unit": "[M]"}, 
                                       Density_power_law_scale: dict[str, float | str] = {"Value": 1.0, "Unit": "[M]"}, 
                                       Temperature_spread: dict[str, float | str] = {"Value": 1.0, "Unit": "[M]"},
                                       Temperature_power_law_power: dict[str, float | str] = {"Value": 0.0, "Unit": "[M]"}, 
                                       Temperature_power_law_scale: dict[str, float | str] = {"Value": 1.0, "Unit": "[M]"}, 
                                       
                                       Density_scale_factor: dict[str, float | str] = {"Value": 1e6, "Unit": "[g/cm^3]"},
                                       Temperature_scale_factor: dict[str, float | str] = {"Value": 1e11, "Unit": "[K]"},
                                       
                                       Mag_field_magnitude_profile: dict[str, str] = {"Value": "Magnetization_based", "Unit": "[-]"},
                                       Mag_field_geometry_r: dict[str, float | str] = {"Value": 0.5, "Unit": "[-]"},
                                       Mag_field_geometry_theta: dict[str, float | str] = {"Value": 0, "Unit": "[-]"},
                                       Mag_field_geometry_phi: dict[str, float | str] = {"Value": 0.87, "Unit": "[-]"},
                                       Mag_field_magnitude_scale: dict[str, float | str] = {"Value": 100, "Unit": "[G]"},
                                       Mag_field_power: dict[str, float | str] = {"Value": 1, "Unit": "[-]"},
                                       Mag_field_radial_scale: dict[str, float | str] = {"Value": 5, "Unit": "[M]"},
                                       Mag_field_geometry: dict[str, str] = {"Value": "Constant", "Unit": "[-]"},
                                       Magnetization: dict[str, float | str] = {"Value": 1.0,  "Unit": "[-]"},
                                       
                                       Coord_time_at_max: dict[str, float | str] = {"Value": 0, "Unit": "[GM/c^3]"},
                                       Distance: dict[str, float | str] = {"Value": 8.0, "Unit": "[M]"},
                                       Inclination: dict[str, float | str] = {"Value": pi / 2,  "Unit": "[Rad]"},
                                       Azimuth: dict[str, float | str] = {"Value": 0,  "Unit": "[Rad]"},
                                       Threshold_relative_density: dict[str, float | str] = {"Value": 1e-3,  "Unit": "[-]"}):

        self.hotspot_model = Hotspot_model()

        self.hotspot_model.Ensamble_type       = Ensamble_type     
        self.hotspot_model.Density_profile     = Density_profile   
        self.hotspot_model.Temperature_profile = Temperature_profile
        
        self.hotspot_model.Velocity_profile         = Velocity_profile
        self.hotspot_model.Radial_velocity_fraction = Radial_velocity_fraction

        self.hotspot_model.Density_scale_factor     = Density_scale_factor     
        self.hotspot_model.Density_power_law_power  = Density_power_law_power
        self.hotspot_model.Density_power_law_scale  = Density_power_law_scale
        self.hotspot_model.Temperature_scale_factor = Temperature_scale_factor  
        self.hotspot_model.Temperature_power_law_power  = Temperature_power_law_power
        self.hotspot_model.Temperature_power_law_scale  = Temperature_power_law_scale

        self.hotspot_model.Mag_field_magnitude_profile = Mag_field_magnitude_profile
        self.hotspot_model.Mag_field_geometry_r = Mag_field_geometry_r
        self.hotspot_model.Mag_field_geometry_theta = Mag_field_geometry_theta
        self.hotspot_model.Mag_field_geometry_phi = Mag_field_geometry_phi
        self.hotspot_model.Mag_field_magnitude_scale = Mag_field_magnitude_scale
        self.hotspot_model.Mag_field_power = Mag_field_power
        self.hotspot_model.Mag_field_radial_scale = Mag_field_radial_scale
        self.hotspot_model.Mag_field_geometry = Mag_field_geometry

        self.hotspot_model.Density_spread     = Density_spread  
        self.hotspot_model.Temperature_spread = Temperature_spread
        self.hotspot_model.Radius             = Radius          
        self.hotspot_model.Temporal_spread    = Temporal_spread
        self.hotspot_model.Coord_time_at_max  = Coord_time_at_max
        self.hotspot_model.Distance           = Distance          
        self.hotspot_model.Inclination        = Inclination       
        self.hotspot_model.Azimuth            = Azimuth           
        self.hotspot_model.Magnetization      = Magnetization     
        
        self.hotspot_model.Threshold_relative_density = Threshold_relative_density         

    def _configure_file_manager(self, Vert_shader_path: str = "C:/Users/Valur/Documents/Repos/Mjolnir_GRRT/Libraries/shaders/default.vert",
                                      Frag_shader_path: str = "C:/Users/Valur/Documents/Repos/Mjolnir_GRRT/Libraries/shaders/default.frag",
                                      Output_file_directory: str = "C:/Users/Valur/Documents/Repos/Mjolnir_GRRT/Sim_Results",
                                      Common_file_names: str = "",
                                      Sim_mode_2_input_file_path: str = "",
                                      Truncate_files: int = 1):                       

        self.file_manager = File_manager()

        self.file_manager.Vert_shader_path = Vert_shader_path
        self.file_manager.Frag_shader_path = Frag_shader_path
        self.file_manager.Output_file_directory = Output_file_directory
        self.file_manager.Common_file_names = Common_file_names
        self.file_manager.Sim_mode_2_input_file_path = Sim_mode_2_input_file_path
        self.file_manager.Truncate_files = Truncate_files

    def generate_simulation_input(self, Path_to_input_dir: str, Input_file_name: str):

        Encoding = 'UTF-8'
        XML_root_node = ET.Element("Simulation_Input", {"Simulation_Name": self.simulation_name["Value"]})
        ET.SubElement(XML_root_node, "Simulation_mode", units = str(self.simulation_mode["Unit"])).text = "{}".format(self.simulation_mode["Value"])
        ET.SubElement(XML_root_node, "Thermalize_emission_medium", units = "[-]").text = "{}".format(self.thermalize_emission_medium["Value"])
        ET.SubElement(XML_root_node, "Average_emission_pitch_angle", units = "[-]").text = "{}".format(self.average_emission_pitch_angle["Value"])
        ET.SubElement(XML_root_node, "Emission_pitch_angle_samples_to_average", units = "[-]").text = "{}".format(self.emission_pitch_angle_samples_to_average["Value"])
        ET.SubElement(XML_root_node, "Central_object_mass", units = str(self.object_mass["Unit"])).text = "{}".format(self.object_mass["Value"])
        ET.SubElement(XML_root_node, "Sim_mode_2_param_value_number", units = str(self.sim_mode_2_param_value_number["Unit"])).text = "{}".format(self.sim_mode_2_param_value_number["Value"])
        ET.SubElement(XML_root_node, "Sim_mode_3_X_init", units = str(self.sim_mode_3_X_init["Unit"])).text = "{}".format(self.sim_mode_3_X_init["Value"])
        ET.SubElement(XML_root_node, "Sim_mode_3_Y_init", units = str(self.sim_mode_3_Y_init["Unit"])).text = "{}".format(self.sim_mode_3_Y_init["Value"])
        ET.SubElement(XML_root_node, "Max_image_order", units = str(self.max_image_order["Unit"])).text = "{}".format(self.max_image_order["Value"])
        
        # ============ Generate the metric XML section ============ #

        Metric_subelement = ET.SubElement(XML_root_node, "Metric")

        match self.metric_parameters.Metric_type["Value"]:
            
            case "Wormhole":
                ET.SubElement(Metric_subelement, "Metric_type", units = "[-]").text = "{}".format("Wormhole")
                ET.SubElement(Metric_subelement, "Spin_parameter", units = "[M]").text = "{}".format(self.metric_parameters.Spin["Value"])
                ET.SubElement(Metric_subelement, "WH_redshift", units = "[-]").text = "{}".format(self.metric_parameters.WH_redshift["Value"])
                ET.SubElement(Metric_subelement, "WH_r_throat", units = "[M]").text = "{}".format(self.metric_parameters.WH_r_throat["Value"])
                ET.SubElement(Metric_subelement, "WH_stop_at_throat", units = "[-]").text = "{}".format(self.metric_parameters.WH_stop_at_throat["Value"])

            case "Minkowski":
                ET.SubElement(Metric_subelement, "Metric_type", units = "[-]").text = "{}".format("Minkowski")

            case "Janis-Newman-Winicour":
                ET.SubElement(Metric_subelement, "Metric_type", units = "[-]").text = "{}".format("Janis-Newman-Winicour")
                ET.SubElement(Metric_subelement, "JNW_gamma", units = "[-]").text = "{}".format(self.metric_parameters.JNW_gamma["Value"])

            case "Einstein-Gauss-Bonnet":
                ET.SubElement(Metric_subelement, "Metric_type", units = "[-]").text = "{}".format("Einstein-Gauss-Bonnet")
                ET.SubElement(Metric_subelement, "EGB_gamma", units = "[-]").text = "{}".format(self.metric_parameters.EGB_gamma["Value"])

            case "Regular-Black-Hole":
                ET.SubElement(Metric_subelement, "Metric_type", units = "[-]").text = "{}".format("Regular-Black-Hole")
                ET.SubElement(Metric_subelement, "RBH_param", units = "[-]").text = "{}".format(self.metric_parameters.RBH_param["Value"])

            case "Black-Hole-w-Dark-Matter":
                ET.SubElement(Metric_subelement, "Metric_type", units = "[-]").text = "{}".format("Black-Hole-w-Dark-Matter")
                ET.SubElement(Metric_subelement, "Halo_compactness", units = "[-]").text = "{}".format(self.metric_parameters.Halo_compactness["Value"])    
                ET.SubElement(Metric_subelement, "Halo_mass", units = "[M]").text = "{}".format(self.metric_parameters.Halo_mass["Value"]) 
                
            case "Numerical":
                ET.SubElement(Metric_subelement, "Metric_type", units = "[-]").text = "{}".format("Numerical")
                ET.SubElement(Metric_subelement, "ADM_Mass", units = "[M]").text = "{}".format(self.metric_parameters.Mass["Value"])
                ET.SubElement(Metric_subelement, "Horizon_radius", units = "[G/c^2]").text = "{}".format(self.metric_parameters.Horizon_radius["Value"])
                ET.SubElement(Metric_subelement, "ADM_ang_momentum", units = "[M]").text = "{}".format(self.metric_parameters.Spin["Value"])
                ET.SubElement(Metric_subelement, "Numerical_metric_anzatz_type", units = "[M]").text = "{}".format(self.metric_parameters.Numerical_metric_anzatz_type["Value"])
                ET.SubElement(Metric_subelement, "Numerical_metric_spline_path").text = "{}".format(self.metric_parameters.Numerical_metric_spline_path)
                
            case _:
                ET.SubElement(Metric_subelement, "Metric_type", units = "[-]").text = "{}".format("Kerr")
                ET.SubElement(Metric_subelement, "Spin_parameter", units = "[M]").text = "{}".format(self.metric_parameters.Spin["Value"])
                
        ET.SubElement(Metric_subelement, "Scattering_radius", units = "[M]").text = "{}".format(self.metric_parameters.Scattering_radius["Value"])
        ET.SubElement(Metric_subelement, "Distance_to_singular_point", units = "[M]").text = "{}".format(self.metric_parameters.Distance_to_singular_point["Value"])
                

        # ============ Generate the observer XML section ============ #

        Observer_subelement = ET.SubElement(XML_root_node, "Observer")
        for Obs_attrib_name in self.observer.__slots__:
            Obs_attrib: dict[str, str | int | float] = getattr(self.observer, Obs_attrib_name)
            ET.SubElement(Observer_subelement, Obs_attrib_name, units = str(Obs_attrib["Unit"])).text = "{}".format(Obs_attrib["Value"])

        # ============ Generate the accretion disk XML section ============ #
        
        Page_Thorne_parameters = ["r_in",
                                  "r_out"]

        Common_RIAF_parameters = ["Opening_angle",
                                  "Density_power_law_scale",
                                  "Density_power_law_power",
                                  "Temperature_power_law_scale",
                                  "Temperature_power_law_power",
                                  "Density_cutoff_radius",
                                  "Density_cutoff_scale",
                                  "Temperature_cutoff_radius",
                                  "Temperature_cutoff_scale"]
        
        Colab_test_1_parameteres = ["Radial_scale",
                                    "Vertical_scale"]

        Common_slots = [slot for slot in self.disk_model.__slots__ if slot not in Common_RIAF_parameters + Colab_test_1_parameteres + ["Radial_velocity_fraction"] + Page_Thorne_parameters + ["Disk_Model"]]

        Disk_subelement = ET.SubElement(XML_root_node, "Accretion_Disk") 
        
        Sub_element = ET.SubElement(Disk_subelement, "Disk_Model", units = "[-]").text = "{}".format(self.disk_model.Disk_Model["Value"])
        
        if self.disk_model.Disk_Model["Value"] == "Page-Thorne":     
            Sub_element = ET.SubElement(Disk_subelement, "r_in", units = "[M]").text = "{}".format(self.disk_model.r_in_PT_disk["Value"])
            Sub_element = ET.SubElement(Disk_subelement, "r_out", units = "[M]").text = "{}".format(self.disk_model.r_out_PT_disk["Value"])
            Sub_element = ET.SubElement(Disk_subelement, "Mag_field_geometry", units = "[-]").text = "{}".format(self.disk_model.Mag_field_geometry["Value"])
            Sub_element = ET.SubElement(Disk_subelement, "Mag_field_geometry_r", units = "[M]").text = "{}".format(self.disk_model.Mag_field_geometry_r["Value"])
            Sub_element = ET.SubElement(Disk_subelement, "Mag_field_geometry_theta", units = "[M]").text = "{}".format(self.disk_model.Mag_field_geometry_theta["Value"])
            Sub_element = ET.SubElement(Disk_subelement, "Mag_field_geometry_phi", units = "[M]").text = "{}".format(self.disk_model.Mag_field_geometry_phi["Value"])
            
        else:
            # ------------- Common subsection
            Common_subelement = ET.SubElement(Disk_subelement, "Common_parameters") 
            for Disk_attrib_name in Common_slots:
                Disk_attrib: dict[str, str | int | float] = getattr(self.disk_model, Disk_attrib_name)
                Sub_element = ET.SubElement(Common_subelement, Disk_attrib_name, units = str(Disk_attrib["Unit"]))
                
                if Disk_attrib_name == "Velocity_profile":
                    ET.SubElement(Sub_element, "Type", units = "-").text = "{}".format(Disk_attrib["Value"])
                    ET.SubElement(Sub_element, "Radial_velocity_fraction", units = "-").text = "{}".format(self.disk_model.Radial_velocity_fraction["Value"])
                else:
                    Sub_element.text = "{}".format(Disk_attrib["Value"])

            match self.disk_model.Disk_Model["Value"]:

                case "Colab_test_1":

                    # ------------- Colab test 1 profile subsection
                    Colab_test_1_subelement = ET.SubElement(Disk_subelement, "Colab_test_1_profile") 
                    for Disk_attrib_name in Colab_test_1_parameteres:
                        Disk_attrib: dict[str, str | int | float] = getattr(self.disk_model, Disk_attrib_name)
                        ET.SubElement(Colab_test_1_subelement, Disk_attrib_name, units = str(Disk_attrib["Unit"])).text = "{}".format(Disk_attrib["Value"])

                case _:

                    # ------------- Phenomenological RIAF subsection
                    Power_law_subelement = ET.SubElement(Disk_subelement, "Common_RIAF_profile") 
                    for Disk_attrib_name in Common_RIAF_parameters:
                        Disk_attrib: dict[str, str | int | float] = getattr(self.disk_model, Disk_attrib_name)
                        ET.SubElement(Power_law_subelement, Disk_attrib_name, units = str(Disk_attrib["Unit"])).text = "{}".format(Disk_attrib["Value"])

        # ============ Generate the hotspot XML section ============ #

        Gaussian_density_slots = ["Density_spread"]
        Gaussian_temperature_slots = ["Temperature_spread"]
        
        Hybrid_density_slots = ["Density_spread", "Density_power_law_power", "Density_power_law_scale"]
        Hybrid_temperature_slots = ["Temperature_spread", "Temperature_power_law_power", "Temperature_power_law_scale"]
        
        Sphere_slots = ["Radius"]

        Common_slots = [slot for slot in self.hotspot_model.__slots__ if slot not in Gaussian_density_slots + 
                                                                                     Gaussian_temperature_slots + 
                                                                                     Sphere_slots + 
                                                                                     Hybrid_density_slots +
                                                                                     Hybrid_temperature_slots + 
                                                                                   ["Radial_velocity_fraction"]]
        
        Hotspot_subelement = ET.SubElement(XML_root_node, "Hotspot") 
        for Hotspot_attrib_name in Common_slots:
            Hotspot_attrib: dict[str, str | int | float] = getattr(self.hotspot_model, Hotspot_attrib_name)
            Sub_element = ET.SubElement(Hotspot_subelement, Hotspot_attrib_name, units = str(Hotspot_attrib["Unit"]))
 
            if Hotspot_attrib_name == "Velocity_profile":
                ET.SubElement(Sub_element, "Type", units = "-").text = "{}".format(Hotspot_attrib["Value"])
                ET.SubElement(Sub_element, "Radial_velocity_fraction", units = "-").text = "{}".format(self.hotspot_model.Radial_velocity_fraction["Value"])
            else:
                Sub_element.text = "{}".format(Hotspot_attrib["Value"])

        match self.hotspot_model.Density_profile["Value"]:

            case "Gaussian":
                # ------------- Gaussian profile subsection
                Gaussian_subelement = ET.SubElement(Hotspot_subelement, "Gaussian_profile") 
                for Hotspot_attrib_name in Gaussian_density_slots:
                    Hotspot_attrib: dict[str, str | int | float] = getattr(self.hotspot_model, Hotspot_attrib_name)
                    ET.SubElement(Gaussian_subelement, Hotspot_attrib_name, units = str(Hotspot_attrib["Unit"])).text = "{}".format(Hotspot_attrib["Value"])

            case "Hybrid_power_law_gaussian":
                # ------------- Hybrid profile subsection
                Hybrid_subelement = ET.SubElement(Hotspot_subelement, "Hybrid_power_law_gaussian_profile") 
                for Hotspot_attrib_name in Hybrid_density_slots:
                    Hotspot_attrib: dict[str, str | int | float] = getattr(self.hotspot_model, Hotspot_attrib_name)
                    ET.SubElement(Hybrid_subelement, Hotspot_attrib_name, units = str(Hotspot_attrib["Unit"])).text = "{}".format(Hotspot_attrib["Value"])
            
            case _:
                # ------------- Spherical profile subsection
                Spherical_subelement = ET.SubElement(Hotspot_subelement, "Spherical_profile") 
                for Hotspot_attrib_name in Sphere_slots:
                    Hotspot_attrib: dict[str, str | int | float] = getattr(self.hotspot_model, Hotspot_attrib_name)
                    ET.SubElement(Spherical_subelement, Hotspot_attrib_name, units = str(Hotspot_attrib["Unit"])).text = "{}".format(Hotspot_attrib["Value"])

        match self.hotspot_model.Temperature_profile["Value"]:

            case "Gaussian":
                # ------------- Gaussian profile subsection
                if "Gaussian_subelement" not in locals():
                    Gaussian_subelement = ET.SubElement(Hotspot_subelement, "Gaussian_profile") 
                
                for Hotspot_attrib_name in Gaussian_temperature_slots:
                    Hotspot_attrib: dict[str, str | int | float] = getattr(self.hotspot_model, Hotspot_attrib_name)
                    ET.SubElement(Gaussian_subelement, Hotspot_attrib_name, units = str(Hotspot_attrib["Unit"])).text = "{}".format(Hotspot_attrib["Value"]) # type: ignore

            case "Hybrid_power_law_gaussian":
                
                if "Hybrid_subelement" not in locals():
                    Hybrid_subelement = ET.SubElement(Hotspot_subelement, "Hybrid_power_law_gaussian_profile") 

                # ------------- Hybrid profile subsection
                for Hotspot_attrib_name in Hybrid_temperature_slots:
                    Hotspot_attrib: dict[str, str | int | float] = getattr(self.hotspot_model, Hotspot_attrib_name)
                    ET.SubElement(Hybrid_subelement, Hotspot_attrib_name, units = str(Hotspot_attrib["Unit"])).text = "{}".format(Hotspot_attrib["Value"]) # type: ignore

            case _:
            
                # ------------- Spherical profile subsectiontry: 
                if "Spherical_subelement" not in locals():
                    
                    Spherical_subelement = ET.SubElement(Hotspot_subelement, "Spherical_profile") 
                    
                    for Hotspot_attrib_name in Sphere_slots:
                        Hotspot_attrib: dict[str, str | int | float] = getattr(self.hotspot_model, Hotspot_attrib_name)
                        ET.SubElement(Spherical_subelement, Hotspot_attrib_name, units = str(Hotspot_attrib["Unit"])).text = "{}".format(Hotspot_attrib["Value"])
    

        # ============ Generate the emission models XML section ============ #

        Emission_subelement = ET.SubElement(XML_root_node, "Emission_models")

        if ((self.hotspot_model.Ensamble_type["Value"] == "Phenomenological" and self.hotspot_model.Density_scale_factor["Value"] != 0) 
            or (self.disk_model.Ensamble_type["Value"] == "Phenomenological" and self.disk_model.Density_scale_factor["Value"] != 0)):

            for Emission_attrib_name in self.emission_models.__slots__:

                attrib_mask = ['Emission_power_law', 'Source_f_power_law', 'Absorbtion_coeff', 'Emission_coeff']

                if Emission_attrib_name in attrib_mask:
                    Emission_attrib: dict[str, str | int | float] = getattr(self.emission_models, Emission_attrib_name)
                    ET.SubElement(Emission_subelement, Emission_attrib_name, units = str(Emission_attrib["Unit"])).text = "{}".format(Emission_attrib["Value"])

        if ((self.hotspot_model.Ensamble_type["Value"] == "Kappa" and self.hotspot_model.Density_scale_factor["Value"] != 0) 
            or (self.disk_model.Ensamble_type["Value"] == "Kappa" and self.disk_model.Density_scale_factor["Value"] != 0)):
            
            ET.SubElement(Emission_subelement, "Kappa", units = "[-]").text = "{}".format(getattr(self.emission_models, "Kappa")["Value"])

        if (self.disk_model.Ensamble_type["Value"] == "Debug_constant_functions"):
            
            attrib_mask = ['Emission_power_law', 'Source_f_power_law', 'Absorbtion_coeff', 'Emission_coeff', 'Kappa']
            
            for Emission_attrib_name in self.emission_models.__slots__:

                if Emission_attrib_name not in attrib_mask:
                    Emission_attrib: dict[str, str | int | float] = getattr(self.emission_models, Emission_attrib_name)
                    ET.SubElement(Emission_subelement, Emission_attrib_name, units = str(Emission_attrib["Unit"])).text = "{}".format(Emission_attrib["Value"])


        # ============ Generate the integrator XML section ============ #

        Integrator_subelement = ET.SubElement(XML_root_node, "Integrator")
        for Integrator_attrib_name in self.integrator.__slots__:
            Integrator_attrib: dict[str, str | int | float] = getattr(self.integrator, Integrator_attrib_name)
            ET.SubElement(Integrator_subelement, Integrator_attrib_name, units = str(Integrator_attrib["Unit"])).text = "{}".format(Integrator_attrib["Value"])

        # ============ Generate the file paths XML section ============ #

        Files_subelement = ET.SubElement(XML_root_node, "File_Manager")
        for Files_attrib_name in self.file_manager.__slots__:
            Files_attrib: dict[str, str | int | float] = getattr(self.file_manager, Files_attrib_name)
            ET.SubElement(Files_subelement, Files_attrib_name).text = "{}".format(Files_attrib)

        # ========================================================== #

        XML_struct = xml.dom.minidom.parseString(ET.tostring(XML_root_node))
        formatted_XML_string = XML_struct.toprettyxml()
        Header, Body = formatted_XML_string.split('?>')

        if not os.path.exists(Path_to_input_dir):
            os.makedirs(Path_to_input_dir)

        with open(Path_to_input_dir + "\\" + Input_file_name, 'w') as xfile:
            xfile.write(Header + 'encoding=\"{}\"?>\n'.format(Encoding) + Body)
            xfile.close()


if __name__ == "__main__":

    Units_class_instance = Units_class()

    Sim_config = Simulation_configurator()
    
    Sim_config.metric_parameters.Numerical_metric_spline_path = "C:/Users/Valur/Documents/Repos/Mjolnir_GRRT/Utilities/Galin_numerical_config_II.XML"

    Sim_config.simulation_mode = {"Value": 1, "Unit": "[-]"}

    Sim_config.object_mass = {"Value": 6.2e9, "Unit": "[M_sun]"}

    # ================================================== Metric ================================================== #

    Sim_config.metric_parameters.Metric_type    = {"Value": "Kerr", "Unit": "[-]"}
    Sim_config.metric_parameters.Mass           = {"Value": 0.881990876889021, "Unit": "[M]"}
    Sim_config.metric_parameters.Horizon_radius = {"Value": 0.01, "Unit": "[G/c^2]"}
    Sim_config.metric_parameters.Spin           = {"Value": 0, "Unit": "[M]"}
    Sim_config.metric_parameters.Numerical_metric_anzatz_type = {"Value": "Anzatz_1", "Unit": "[M]"} 
    
    Sim_config.metric_parameters.Scattering_radius = {"Value": 400, "Unit": "[M]"} 
    
    # ================================================== Observer ================================================== #

    Sim_config.observer.Resolution_x = {"Value": 1024, "Unit": "[-]"}
    Sim_config.observer.Resolution_y = {"Value": 1024, "Unit": "[-]"}
    
    Sim_config.observer.Distance    = {"Value": 200, "Unit": "[M]"}
    Sim_config.observer.Inclination = {"Value": 90 * pi / 180, "Unit": "[Rad]"}
    Sim_config.observer.Obs_frequency = {"Value": 230e9, "Unit": "[Hz]"}
    Sim_config.observer.Cam_rotation_angle = {"Value": 0, "Unit": "[Hz]"}

    # ================================================== Disk ================================================== #
    
    Sim_config.disk_model.Ensamble_type = {"Value": "Thermal",   "Unit": "[-]"}
    Sim_config.disk_model.Disk_Model    = {"Value": "Phenom_RIAF_1", "Unit": "[-]"}
    Sim_config.disk_model.Mag_field_geometry = {"Value": "Constant", "Unit": "[-]"}
    
    Sim_config.disk_model.Density_scale_factor = {"Value": 500000, "Unit": "[g/cm^3]"}
    Sim_config.disk_model.Temperature_scale_factor = {"Value": 4.1e+10, "Unit": "[K]"}
            
    Sim_config.disk_model.Density_cutoff_radius = {"Value": 5, "Unit": "[M]"}
    Sim_config.disk_model.Temperature_cutoff_radius = {"Value": 5, "Unit": "[M]"}

    Sim_config.disk_model.Density_power_law_scale     = {"Value": 5, "Unit": "[M]"}
    Sim_config.disk_model.Temperature_power_law_scale = {"Value": 5, "Unit": "[M]"}

    Sim_config.disk_model.Opening_angle = {"Value": 0.4, "Unit": "[tan(angle)]"}
    
    Sim_config.disk_model.Density_power_law_power     = {"Value": 2.0, "Unit": "[-]"}
    Sim_config.disk_model.Temperature_power_law_power = {"Value": 1.0, "Unit": "[-]"}
    
    Sim_config.disk_model.Velocity_profile = {"Value": "Theta Dependant", "Unit": "[-]"}
    
    Sim_config.observer.Image_y_min = {"Value": -10, "Unit": "[M]"}
    Sim_config.observer.Image_y_max = {"Value":  10, "Unit": "[M]"}
    Sim_config.observer.Image_x_min = {"Value": -10, "Unit": "[M]"}
    Sim_config.observer.Image_x_max = {"Value":  10, "Unit": "[M]"}
        
    Sim_config.integrator.RK78_accuracy      = {"Value": 1e-10, "Unit": "[-]"}
    Sim_config.observer.Include_polarization = {"Value": 0, "Unit": "[-]"}
    Sim_config.integrator.Step_controller_type = {"Value": "PID", "Unit": "[-]"}
    Sim_config.integrator.max_integration_count = {"Value": 1000000, "Unit": "[-]"}
    Sim_config.integrator.max_affine_parameter = {"Value": 1000000, "Unit": "[-]"}
    Sim_config.metric_parameters.Distance_to_singular_point = {"Value": 1e-2, "Unit": "[M]"}
    
    Sim_config.integrator.Max_rel_step_increase = {"Value": 5, "Unit": "[-]"}
    # ================================================== Hotspot ================================================== #

    Sim_config.hotspot_model.Density_scale_factor = {"Value": 0, "Unit": "[g / cm^3]"}
    
    """ The simulation name and input file path """
    Sim_config.simulation_name = {"Value": "Reference_Simulation_2", "Unit": "[-]"}

    """ The simulation output file path """
    Sim_config.file_manager.Output_file_directory = parent_directory + "Reference_simulations"
    Sim_config.simulation_name = {"Value": "Old_wormhole_sanity_check", "Unit": "[-]"}

    Sim_config.generate_simulation_input(Path_to_input_dir = "Reference_simulations\\Old_wormhole_sanity_check",
                                         Input_file_name = "Old_wormhole_sanity_check.XML")

    import subprocess

    filename = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\Utilities\\Reference_simulations\\Old_wormhole_sanity_check\\Old_wormhole_sanity_check.xml"
    args = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\x64\\Release\\Mjolnir_GRRT.exe -in " + filename + " -print_to_console 1"
    
    subprocess.call(args, shell = True)