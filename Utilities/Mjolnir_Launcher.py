from numpy import pi
import numpy as np

import xml.etree.cElementTree as ET
import xml.dom.minidom

from Support_functions.Parsers import Units_class

class Integrator():

    __slots__ = ("init_stepsize",
                 "RK45_accuracy", 
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
                 "simpson_method_accuracy")

class Disk_model():

    __slots__ = ("Ensamble_type",
                 "Density_profile",
                 "Temperature_profile",
                 "Velocity_profile",
                 "Density_scale_factor",
                 "Temperature_scale_factor",
                 "Mag_field_geometry_X",
                 "Mag_field_geometry_Y",
                 "Mag_field_geometry_Z",
                 "Opening_angle",

                 "Density_radial_power_law",
                 "Density_cutoff_scale",
                 "Density_r_cutoff",
                 "Density_r_0",
                 "Temperature_r_0",
                 "Temperature_r_cutoff",
                 "Temperature_cutoff_scale",
                 "Temperature_radial_power_law",

                 "Magnetization",
                 
                 "Density_exp_height_scale",
                 "Density_exp_radial_scale",
                 "Temperature_exp_height_scale",
                 "Temperature_exp_radial_scale")
    
class Hotspot_model():

    __slots__ = ("Ensamble_type",
                 "Density_profile",
                 "Temperature_profile",
                 "Velocity_profile",
                 "Density_scale_factor",
                 "Temperature_scale_factor",
                 "Radius",
                 "Mag_field_geometry_X",
                 "Mag_field_geometry_Y",
                 "Mag_field_geometry_Z",
                 "Density_spread",
                 "Temperature_spread",
                 "Temporal_spread",
                 "Coord_time_at_max",
                 "Distance",
                 "Inclination",
                 "Azimuth",
                 "Magnetization")

class Metric_parameters():

    __slots__ = ("Mass",
                 "Spin",
                 "WH_redshift",
                 "WH_r_throat", 
                 "WH_stop_at_throat",
                 "RBH_param", 
                 "JNW_gamma", 
                 "EGB_gamma", 
                 "Halo_compactness", 
                 "Halo_mass",
                 "Metric_type")

class Observer():

    __slots__ = ("Distance",
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
                 "Kappa")  

class NT_model_params():

    __slots__ = ("r_in", 
                 "r_out", 
                 "Evaluate_NT_disk")

class File_manager():

    __slots__ = ("Sim_mode_2_input_file_path", "Output_file_directory", "Common_file_names", "Vert_shader_path", "Frag_shader_path", "Truncate_files")

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
                 "emission_pitch_angle_samples_to_average",
                 "object_mass",
                 "simulation_mode",
                 "sim_mode_2_param_value_number",
                 "sim_mode_3_X_init",
                 "sim_mode_3_Y_init")

    def __init__(self, 
                 Average_emission_pitch_angle: dict = {"Value": 1, "Unit": "[-]"}, 
                 emission_pitch_angle_samples_to_average: dict = {"Value": 50, "Units": "[-]"},
                 object_mass: dict = {"Value": 6.2e9, "Unit": "[M_sun]"},
                 simulation_name: dict = {"Value": "Test_Simulation", "Unit": "[-]"},
                 simulation_mode: dict = {"Value": 1, "Unit": "[-]"}, 
                 sim_mode_2_param_value_number: dict = {"Value": 1, "Unit": "[-]"},
                 sim_mode_3_X_init: dict = {"Value": 1, "Unit": "[M]"},
                 sim_mode_3_Y_init: dict = {"Value": 1, "Unit": "[M]"}):

        self.average_emission_pitch_angle = Average_emission_pitch_angle
        self.emission_pitch_angle_samples_to_average = emission_pitch_angle_samples_to_average
        self.simulation_name = simulation_name
        self.object_mass = object_mass
        self.simulation_mode = simulation_mode
        self.sim_mode_2_param_value_number = sim_mode_2_param_value_number
        self.sim_mode_3_X_init = sim_mode_3_X_init
        self.sim_mode_3_Y_init = sim_mode_3_Y_init

        self._configure_integrator_settings()
        self._configure_observer()
        self._configure_NT_model()
        self._configure_disk_model()
        self._configure_hotspot_model()
        self._configure_file_manager()
        self._configure_metric_parameters()
        self._configure_emission_models()

    def _configure_integrator_settings(self, Init_stepsize: dict = {"Value": 1e-5, "Unit": "[M]"},
                                             RK45_accuracy: dict = {"Value": 1e-13, "Unit": "[-]"},
                                             Step_controller_type: dict = {"Value": "Gustafsson", "Unit": "[-]"},
                                             Safety_factor_1: dict = {"Value": 0.8, "Unit": "[-]"},
                                             Safety_factor_2: dict = {"Value": 1e-25, "Unit": "[-]"},
                                             Max_rel_step_increase: dict = {"Value": 20, "Unit": "[-]"},
                                             Min_rel_step_increase: dict = {"Value": 0.1, "Unit": "[-]"},
                                             Step_controller_I_gain: dict = {"Value": 0.117, "Unit": "[-]"},
                                             Step_controller_P_gain: dict = {"Value": -0.042, "Unit": "[-]"},
                                             Step_controller_D_gain: dict = {"Value": 0.02, "Unit": "[-]"},
                                             Gustafsson_controller_k_1: dict = {"Value": 0.0734, "Unit": "[-]"},
                                             Gustafsson_controller_k_2: dict = {"Value": 0.1136, "Unit": "[-]"},
                                             Max_integration_count: dict = {"Value": 1e7, "Unit": "[-]"},
                                             simpson_method_accuracy: dict = {"Value": 1e-6, "Unit": "[-]"}):

        self.integrator = Integrator()

        self.integrator.init_stepsize = Init_stepsize
        self.integrator.init_stepsize = Init_stepsize
        self.integrator.RK45_accuracy = RK45_accuracy
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

    def _configure_observer(self, Distance: dict = {"Value": 1e4, "Unit": "[M]"},
                                  Inclination: dict = {"Value": 160 / 180 * pi, "Unit": "[Rad]"},
                                  Azimuth: dict = {"Value": 0.0, "Unit": "[Rad]"},
                                  Cam_rotation_angle: dict = {"Value": 0.0, "Unit": "[Rad]"},
                                  Image_y_min: dict = {"Value": -15, "Unit": "[M]"},
                                  Image_y_max: dict = {"Value":  15, "Unit": "[M]"},
                                  Image_x_min: dict = {"Value": -15, "Unit": "[M]"},
                                  Image_x_max: dict = {"Value":  15, "Unit": "[M]"},
                                  Resolution_y: dict = {"Value": 2048, "Unit": "[-]"},
                                  Resolution_x: dict = {"Value": 2048, "Unit": "[-]"},
                                  Observation_frequency: dict = {"Value": 230e9, "Unit": "[Hz]"},
                                  Include_polarization: dict = {"Value": 0, "Unit": "[-]"}):
        
        self.observer = Observer()

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

    def _configure_metric_parameters(self, Mass: dict = {"Value": 1.0, "Unit": "[M]"}, 
                                           Spin: dict = {"Value": 0.98, "Unit": "[M]"}, 
                                           WH_redshift: dict = {"Value": 2.0, "Unit": "[M]"},
                                           WH_r_throat: dict = {"Value": 1.0, "Unit": "[M]"}, 
                                           WH_stop_at_throat: dict = {"Value": 0, "Unit": "[-]"}, 
                                           RBH_param: dict = {"Value": 0.5, "Unit": "[M]"}, 
                                           JNW_gamma: dict = {"Value": 0.48, "Unit": "[-]"}, 
                                           EGB_gamma: dict = {"Value": 1.15, "Unit": "[M^2]"}, 
                                           Halo_compactness: dict = {"Value": 1e-4, "Unit": "[-]"},
                                           Halo_mass: dict = {"Value": 1e4, "Unit": "[M]"},
                                           Metric_type: dict = {"Value": "Kerr", "Unit": "[-]"}):

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

    def _configure_NT_model(self, r_in: dict = {"Value": 6, "Unit": "[M]"},
                                  r_out: dict = {"Value": 50, "Unit": "[M]"},
                                  evaluate_NT_disk: dict = {"Value": 1, "Unit": "[-]"}):
        
        self.NT_model_params = NT_model_params()

        self.NT_model_params.r_in = r_in
        self.NT_model_params.r_out = r_out
        self.NT_model_params.Evaluate_NT_disk = evaluate_NT_disk

    def _configure_emission_models(self, Emission_power_law: dict = {"Value": 0.0, "Unit": "[-]"},
                                         Source_f_power_law: dict = {"Value": 2.5, "Unit": "[-]"},
                                         Absorbtion_coeff: dict = {"Value": 1e5, "Unit": "[?]"},
                                         Emission_coeff: dict = {"Value": 3e-18, "Unit": "[erg / (cm^3 s sr Hz)]"},
                                         Kappa: dict = {"Value": 4.0, "Unit": "[-]"}):
        
        self.emission_models = Emission_models()

        self.emission_models.Emission_power_law = Emission_power_law
        self.emission_models.Source_f_power_law = Source_f_power_law
        self.emission_models.Absorbtion_coeff   = Absorbtion_coeff
        self.emission_models.Emission_coeff     = Emission_coeff
        self.emission_models.Kappa              = Kappa

    def _configure_disk_model(self, Ensamble_type: dict = {"Value": "Thermal", "Unit": "[-]"},
                                    Density_profile: dict = {"Value": "Power Law", "Unit": "[-]"},
                                    Temperature_profile: dict = {"Value": "Power Law", "Unit": "[-]"},
                                    Velocity_profile: dict = {"Value": "Theta Dependant", "Unit": "[-]"},
                                    Density_scale_factor: dict = {"Value": 1e5, "Unit": "[g/cm^3]"},
                                    Temperature_scale_factor: dict = {"Value": 1e11, "Unit": "[K]"},
                                    Mag_field_geometry_X: dict = {"Value": 0.5, "Unit": "[-]"},
                                    Mag_field_geometry_Y: dict = {"Value": 0, "Unit": "[-]"},
                                    Mag_field_geometry_Z: dict = {"Value": 0.87, "Unit": "[-]"},
                                    Opening_angle: dict = {"Value": 0.1, "Unit": "[tan(angle)]"},
                                    Density_radial_power_law : dict = {"Value": 2.0, "Unit": "[-]"},
                                    Density_cutoff_scale: dict = {"Value": 0.4, "Unit": "[M]"},
                                    Density_r_cutoff: dict = {"Value": 5.0, "Unit": "[M]"},
                                    Density_r_0: dict = {"Value": 5.0, "Unit": "[M]"},
                                    Magnetization: dict = {"Value": 0.01, "Unit": "[-]"},
                                    Temp_r_0: dict = {"Value": 5.0, "Unit": "[M]"},
                                    Temp_r_cutoff: dict = {"Value": 5.0, "Unit": "[M]"},
                                    Temp_cutoff_scale: dict = {"Value": 0.4, "Unit": "[M]"},
                                    Temp_radial_power_law: dict = {"Value": 1.0, "Unit": "[-]"},
                                    Density_exp_height_scale: dict = {"Value": 0.3, "Unit": "[M]"},
                                    Density_exp_radial_scale: dict = {"Value": 10.0, "Unit": "[M]"},
                                    Temperature_exp_height_scale: dict = {"Value": 0.3, "Unit": "[M]"},
                                    Temperature_exp_radial_scale: dict = {"Value": 10.0, "Unit": "[M]"}):
        
        self.disk_model = Disk_model()

        self.disk_model.Ensamble_type       = Ensamble_type
        self.disk_model.Density_profile     = Density_profile
        self.disk_model.Temperature_profile = Temperature_profile
        self.disk_model.Velocity_profile    = Velocity_profile

        self.disk_model.Density_scale_factor     = Density_scale_factor
        self.disk_model.Temperature_scale_factor = Temperature_scale_factor

        self.disk_model.Mag_field_geometry_X = Mag_field_geometry_X
        self.disk_model.Mag_field_geometry_Y = Mag_field_geometry_Y
        self.disk_model.Mag_field_geometry_Z = Mag_field_geometry_Z
        self.disk_model.Magnetization        = Magnetization

        self.disk_model.Opening_angle            = Opening_angle
        self.disk_model.Density_r_0              = Density_r_0
        self.disk_model.Density_r_cutoff         = Density_r_cutoff
        self.disk_model.Density_cutoff_scale     = Density_cutoff_scale
        self.disk_model.Density_radial_power_law = Density_radial_power_law

        self.disk_model.Density_exp_height_scale = Density_exp_height_scale
        self.disk_model.Density_exp_radial_scale = Density_exp_radial_scale
        self.disk_model.Temperature_exp_height_scale = Temperature_exp_height_scale
        self.disk_model.Temperature_exp_radial_scale = Temperature_exp_radial_scale

        self.disk_model.Temperature_r_0              = Temp_r_0
        self.disk_model.Temperature_r_cutoff         = Temp_r_cutoff
        self.disk_model.Temperature_cutoff_scale     = Temp_cutoff_scale
        self.disk_model.Temperature_radial_power_law = Temp_radial_power_law

    def _configure_hotspot_model(self, Ensamble_type: dict = {"Value": "Kappa", "Unit": "[-]"},
                                       Density_profile: dict = {"Value": "Gaussian", "Unit": "[-]"},
                                       Temperature_profile: dict = {"Value": "Gaussian", "Unit": "[-]"},
                                       Velocity_profile: dict = {"Value": "Theta Dependant", "Unit": "[-]"},
                                       Radius: dict = {"Value": 1, "Unit": "[M]"},
                                       Density_scale_factor: dict = {"Value": 1e6, "Unit": "[g/cm^3]"},
                                       Temperature_scale_factor: dict = {"Value": 1e11, "Unit": "[K]"},
                                       Mag_field_geometry_X: dict = {"Value": 0.5, "Unit": "[-]"},
                                       Mag_field_geometry_Y: dict = {"Value": 0, "Unit": "[-]"},
                                       Mag_field_geometry_Z: dict = {"Value": 0.87, "Unit": "[-]"},
                                       Density_spread: dict = {"Value": 1.0, "Unit": "[M]"},
                                       Temperature_spread: dict = {"Value": 1.0, "Unit": "[M]"},
                                       Temporal_spread: dict = {"Value": 85, "Unit": "[GM/c^3]"},
                                       Coord_time_at_max: dict = {"Value": 1000, "Unit": "[GM/c^3]"},
                                       Distance: dict = {"Value": 8.0, "Unit": "[M]"},
                                       Inclination: dict = {"Value": pi / 2,  "Unit": "[Rad]"},
                                       Azimuth: dict = {"Value": -pi / 2,  "Unit": "[Rad]"},
                                       Magnetization: dict = {"Value": 1.0,  "Unit": "[-]"}):

        self.hotspot_model = Hotspot_model()

        self.hotspot_model.Ensamble_type       = Ensamble_type     
        self.hotspot_model.Density_profile     = Density_profile   
        self.hotspot_model.Temperature_profile = Temperature_profile
        self.hotspot_model.Velocity_profile    = Velocity_profile

        self.hotspot_model.Density_scale_factor     = Density_scale_factor     
        self.hotspot_model.Temperature_scale_factor = Temperature_scale_factor 

        self.hotspot_model.Mag_field_geometry_X = Mag_field_geometry_X
        self.hotspot_model.Mag_field_geometry_Y = Mag_field_geometry_Y
        self.hotspot_model.Mag_field_geometry_Z = Mag_field_geometry_Z

        self.hotspot_model.Density_spread     = Density_spread  
        self.hotspot_model.Temperature_spread = Temperature_spread
        self.hotspot_model.Radius             = Radius          
        self.hotspot_model.Temporal_spread    = Temporal_spread
        self.hotspot_model.Coord_time_at_max  = Coord_time_at_max
        self.hotspot_model.Distance           = Distance          
        self.hotspot_model.Inclination        = Inclination       
        self.hotspot_model.Azimuth            = Azimuth           
        self.hotspot_model.Magnetization      = Magnetization     

    def _configure_file_manager(self, Vert_shader_path: str = "C:/Users/Valur/Documents/Repos/Mjolnir_GRRT/Libraries/shaders/default.vert",
                                      Frag_shader_path: str = "C:/Users/Valur/Documents/Repos/Mjolnir_GRRT/Libraries/shaders/default.frag",
                                      Output_file_directory: str = "C:/Users/Valur/Documents/Repos/Mjolnir_GRRT/Sim_Results",
                                      Common_file_names: str = "",
                                      Sim_mode_2_input_file_path: str = "",
                                      Truncate_files: bool = 1):
                            

        self.file_manager = File_manager()

        self.file_manager.Vert_shader_path = Vert_shader_path
        self.file_manager.Frag_shader_path = Frag_shader_path
        self.file_manager.Output_file_directory = Output_file_directory
        self.file_manager.Common_file_names = Common_file_names
        self.file_manager.Sim_mode_2_input_file_path = Sim_mode_2_input_file_path
        self.file_manager.Truncate_files = Truncate_files

    def generate_simulation_input(self):

        Encoding = 'UTF-8'
        XML_root_node = ET.Element("Simulation_Input", {"Simulation_Name": self.simulation_name["Value"]})
        ET.SubElement(XML_root_node, "Simulation_mode", units = self.simulation_mode["Unit"]).text = "{}".format(self.simulation_mode["Value"])
        ET.SubElement(XML_root_node, "Average_emission_pitch_angle", units = "[-]").text = "{}".format(self.average_emission_pitch_angle["Value"])
        ET.SubElement(XML_root_node, "Emission_pitch_angle_samples_to_average", units = "[-]").text = "{}".format(self.emission_pitch_angle_samples_to_average["Value"])
        ET.SubElement(XML_root_node, "Central_object_mass", units = self.object_mass["Unit"]).text = "{}".format(self.object_mass["Value"])
        ET.SubElement(XML_root_node, "Sim_mode_2_param_value_number", units = self.sim_mode_2_param_value_number["Unit"]).text = "{}".format(self.sim_mode_2_param_value_number["Value"])
        ET.SubElement(XML_root_node, "Sim_mode_3_X_init", units = self.sim_mode_3_X_init["Unit"]).text = "{}".format(self.sim_mode_3_X_init["Value"])
        ET.SubElement(XML_root_node, "Sim_mode_3_Y_init", units = self.sim_mode_3_Y_init["Unit"]).text = "{}".format(self.sim_mode_3_Y_init["Value"])

        # ============ Generate the observer XML section ============ #

        Metric_subelement = ET.SubElement(XML_root_node, "Metric")

        match self.metric_parameters.Metric_type["Value"]:
            
            case "Wormhole":
                ET.SubElement(Metric_subelement, "Metric_type", units = "[-]").text = "{}".format("Wormhole")
                ET.SubElement(Metric_subelement, "Spin_parameter", units = "[M]").text = "{}".format(self.metric_parameters.Spin["Value"])
                ET.SubElement(Metric_subelement, "WH_redshift", units = "[-]").text = "{}".format(self.metric_parameters.WH_redshift["Value"])
                ET.SubElement(Metric_subelement, "WH_r_throat", units = "[M]").text = "{}".format(self.metric_parameters.WH_r_throat["Value"])
                ET.SubElement(Metric_subelement, "WH_stop_at_throat", units = "[-]").text = "{}".format(self.metric_parameters.WH_stop_at_throat["Value"])

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

            case _:
                ET.SubElement(Metric_subelement, "Metric_type", units = "[-]").text = "{}".format("Kerr")
                ET.SubElement(Metric_subelement, "Spin_parameter", units = "[M]").text = "{}".format(self.metric_parameters.Spin["Value"])

        # ============ Generate the observer XML section ============ #

        Observer_subelement = ET.SubElement(XML_root_node, "Observer")
        for Obs_attrib_name in self.observer.__slots__:
            Obs_attrib = getattr(self.observer, Obs_attrib_name)
            ET.SubElement(Observer_subelement, Obs_attrib_name, units = Obs_attrib["Unit"]).text = "{}".format(Obs_attrib["Value"])

        # ============ Generate the accretion disk XML section ============ #

        Density_Power_law_slots = ["Density_radial_power_law",
                                   "Density_cutoff_scale",
                                   "Density_r_cutoff",
                                   "Density_r_0",
                                   "Temperature_radial_power_law",
                                   "Temperature_cutoff_scale",
                                   "Temperature_r_cutoff",
                                   "Temperature_r_0",
                                   "Opening_angle"]
        
        Temperature_Power_law_slots = ["Temperature_radial_power_law",
                                       "Temperature_cutoff_scale",
                                       "Temperature_r_cutoff",
                                       "Temperature_r_0"]
        
        Density_exponential_law_slots = ["Density_exp_height_scale",
                                        "Density_exp_radial_scale"]
        
        Temperature_exponential_law_slots = ["Temperature_exp_height_scale",
                                             "Temperature_exp_radial_scale"]

        Common_slots = [slot for slot in self.disk_model.__slots__ if slot not in Density_Power_law_slots + Temperature_Power_law_slots + Density_exponential_law_slots + Temperature_exponential_law_slots]

        Disk_subelement = ET.SubElement(XML_root_node, "Accretion_Disk") 

        # ------------- Common subsection
        Common_subelement = ET.SubElement(Disk_subelement, "Common_parameters") 
        for Disk_attrib_name in Common_slots:
            Disk_attrib = getattr(self.disk_model, Disk_attrib_name)
            ET.SubElement(Common_subelement, Disk_attrib_name, units = Disk_attrib["Unit"]).text = "{}".format(Disk_attrib["Value"])

        match self.disk_model.Density_profile["Value"]:

            case "Exponential Law":

                # ------------- Exponential law subsection
                Exponential_law_subelement = ET.SubElement(Disk_subelement, "Exponential_law_profile") 
                for Disk_attrib_name in Density_exponential_law_slots:
                    Disk_attrib = getattr(self.disk_model, Disk_attrib_name)
                    ET.SubElement(Exponential_law_subelement, Disk_attrib_name, units = Disk_attrib["Unit"]).text = "{}".format(Disk_attrib["Value"])

            case _:

                # ------------- Power law subsection
                Power_law_subelement = ET.SubElement(Disk_subelement, "Power_law_profile") 
                for Disk_attrib_name in Density_Power_law_slots:
                    Disk_attrib = getattr(self.disk_model, Disk_attrib_name)
                    ET.SubElement(Power_law_subelement, Disk_attrib_name, units = Disk_attrib["Unit"]).text = "{}".format(Disk_attrib["Value"])

        match self.disk_model.Temperature_profile["Value"]:

            case "Exponential Law":

                # ------------- Exponential law subsection

                try: 
                    Exponential_law_subelement
                except:
                    Exponential_law_subelement = ET.SubElement(Disk_subelement, "Exponential_law_profile") 

                for Disk_attrib_name in Temperature_exponential_law_slots:
                    Disk_attrib = getattr(self.disk_model, Disk_attrib_name)
                    ET.SubElement(Exponential_law_subelement, Disk_attrib_name, units = Disk_attrib["Unit"]).text = "{}".format(Disk_attrib["Value"])

            case _:

                # ------------- Power law subsection

                try: 
                    Power_law_subelement
                except:
                    Power_law_subelement = ET.SubElement(Disk_subelement, "Power_law_profile") 

                for Disk_attrib_name in Temperature_Power_law_slots:
                    Disk_attrib = getattr(self.disk_model, Disk_attrib_name)
                    ET.SubElement(Power_law_subelement, Disk_attrib_name, units = Disk_attrib["Unit"]).text = "{}".format(Disk_attrib["Value"])


        # ============ Generate the hotspot XML section ============ #

        Gaussian_density_slots = ["Density_spread"]
        Gaussian_temperature_slots = ["Temperature_spread"]
        
        Sphere_slots = ["Radius"]

        Common_slots = [slot for slot in self.hotspot_model.__slots__ if slot not in Gaussian_density_slots + Gaussian_temperature_slots + Sphere_slots]
        
        Hotspot_subelement = ET.SubElement(XML_root_node, "Hotspot") 
        for Hotspot_attrib_name in Common_slots:
            Hotspot_attrib = getattr(self.hotspot_model, Hotspot_attrib_name)
            ET.SubElement(Hotspot_subelement, Hotspot_attrib_name, units = Hotspot_attrib["Unit"]).text = "{}".format(Hotspot_attrib["Value"])

        match self.hotspot_model.Density_profile["Value"]:

            case "Gaussian":
                # ------------- Gaussian profile subsection
                Gaussian_subelement = ET.SubElement(Hotspot_subelement, "Gaussian_profile") 
                for Hotspot_attrib_name in Gaussian_density_slots:
                    Hotspot_attrib = getattr(self.hotspot_model, Hotspot_attrib_name)
                    ET.SubElement(Gaussian_subelement, Hotspot_attrib_name, units = Hotspot_attrib["Unit"]).text = "{}".format(Hotspot_attrib["Value"])

            case _:
            
                # ------------- Spherical profile subsection
                Spherical_subelement = ET.SubElement(Hotspot_subelement, "Spherical_profile") 
                for Hotspot_attrib_name in Sphere_slots:
                    Hotspot_attrib = getattr(self.hotspot_model, Hotspot_attrib_name)
                    ET.SubElement(Spherical_subelement, Hotspot_attrib_name, units = Hotspot_attrib["Unit"]).text = "{}".format(Hotspot_attrib["Value"])

        match self.hotspot_model.Temperature_profile["Value"]:

            case "Gaussian":
                # ------------- Gaussian profile subsection
                try: 
                    Gaussian_subelement
                except:
                    Gaussian_subelement = ET.SubElement(Hotspot_subelement, "Gaussian_profile") 
                
                for Hotspot_attrib_name in Gaussian_temperature_slots:
                    Hotspot_attrib = getattr(self.hotspot_model, Hotspot_attrib_name)
                    ET.SubElement(Gaussian_subelement, Hotspot_attrib_name, units = Hotspot_attrib["Unit"]).text = "{}".format(Hotspot_attrib["Value"])

            case _:
            
                # ------------- Spherical profile subsectiontry: 
                try:
                    Spherical_subelement
                except:
                    Spherical_subelement = ET.SubElement(Hotspot_subelement, "Spherical_profile") 
                    
                    for Hotspot_attrib_name in Sphere_slots:
                        Hotspot_attrib = getattr(self.hotspot_model, Hotspot_attrib_name)
                        ET.SubElement(Spherical_subelement, Hotspot_attrib_name, units = Hotspot_attrib["Unit"]).text = "{}".format(Hotspot_attrib["Value"])
    
        # ============ Generate the emission models XML section ============ #

        Emission_subelement = ET.SubElement(XML_root_node, "Emission_models")

        if ((self.hotspot_model.Ensamble_type["Value"] == "Phenomenological" and self.hotspot_model.Density_scale_factor["Value"] != 0) 
            or (self.disk_model.Ensamble_type["Value"] == "Phenomenological" and self.disk_model.Density_scale_factor["Value"] != 0)):

            for Emission_attrib_name in self.emission_models.__slots__:

                if Emission_attrib_name != "Kappa":
                    Emission_attrib = getattr(self.emission_models, Emission_attrib_name)
                    ET.SubElement(Emission_subelement, Emission_attrib_name, units = Emission_attrib["Unit"]).text = "{}".format(Emission_attrib["Value"])

        if ((self.hotspot_model.Ensamble_type["Value"] == "Kappa" and self.hotspot_model.Density_scale_factor["Value"] != 0) 
            or (self.disk_model.Ensamble_type["Value"] == "Kappa" and self.disk_model.Density_scale_factor["Value"] != 0)):
            
            ET.SubElement(Emission_subelement, "Kappa", units = "[-]").text = "{}".format(getattr(self.emission_models, "Kappa")["Value"])

        # ============ Generate the Novikov-Thorne models XML section ============ #

        NT_subelement = ET.SubElement(XML_root_node, "Novikov_Thorne_disk")

        NT_attrib = getattr(self.NT_model_params, "Evaluate_NT_disk")
        ET.SubElement(NT_subelement, "Evaluate_NT_disk", units = NT_attrib["Unit"]).text = "{}".format(NT_attrib["Value"])

        if self.NT_model_params.Evaluate_NT_disk["Value"]:

            NT_attrib = getattr(self.NT_model_params, "r_in")
            ET.SubElement(NT_subelement, "r_in", units = NT_attrib["Unit"]).text = "{}".format(NT_attrib["Value"])
            NT_attrib = getattr(self.NT_model_params, "r_out")
            ET.SubElement(NT_subelement, "r_out", units = NT_attrib["Unit"]).text = "{}".format(NT_attrib["Value"])

        # ============ Generate the integrator XML section ============ #

        Integrator_subelement = ET.SubElement(XML_root_node, "Integrator")
        for Integrator_attrib_name in self.integrator.__slots__:
            Integrator_attrib = getattr(self.integrator, Integrator_attrib_name)
            ET.SubElement(Integrator_subelement, Integrator_attrib_name, units = Integrator_attrib["Unit"]).text = "{}".format(Integrator_attrib["Value"])

        # ============ Generate the file paths XML section ============ #

        Files_subelement = ET.SubElement(XML_root_node, "File_Manager")
        for Files_attrib_name in self.file_manager.__slots__:
            Files_attrib = getattr(self.file_manager, Files_attrib_name)
            ET.SubElement(Files_subelement, Files_attrib_name).text = "{}".format(Files_attrib)

        # ========================================================== #

        XML_struct = xml.dom.minidom.parseString(ET.tostring(XML_root_node))
        formatted_XML_string = XML_struct.toprettyxml()
        Header, Body = formatted_XML_string.split('?>')

        with open("FILE.xml", 'w') as xfile:
            xfile.write(Header + 'encoding=\"{}\"?>\n'.format(Encoding) + Body)
            xfile.close()


Units_class_instance = Units_class()

Sim_config = Simulation_configurator()

Sim_config.object_mass = {"Value": 4.2e6, "Unit": "[M_sun]"}

# ================================================== Metric ================================================== #

Sim_config.metric_parameters.Metric_type = {"Value": "Wormhole", "Unit": "[-]"}
Sim_config.metric_parameters.Spin = {"Value": 0.9, "Unit": "[M]"}

# ================================================== Observer ================================================== #

Sim_config.observer.Resolution_x = {"Value": 256, "Unit": "[-]"}
Sim_config.observer.Resolution_y = {"Value": 256, "Unit": "[-]"}
Sim_config.observer.Distance = {"Value": 1e4, "Unit": "[M]"}
Sim_config.observer.Inclination = {"Value": 80 * pi / 180, "Unit": "[Rad]"}
Sim_config.observer.Obs_frequency = {"Value": 230e9, "Unit": "[Hz]"}

# ================================================== Disk ================================================== #
Sim_config.disk_model.Ensamble_type = {"Value": "Phenomenological", "Unit": "[-]"}
Sim_config.disk_model.Density_profile = {"Value": "Exponential Law", "Unit": "[-]"}
Sim_config.disk_model.Temperature_profile = {"Value": "Exponential Law", "Unit": "[-]"}
Sim_config.disk_model.Temperature_scale_factor = {"Value": 5.85e10, "Unit": "[K]"}
Sim_config.disk_model.Density_scale_factor = {"Value": 500000, "Unit": "[g / cm^3]"}

Sim_config.disk_model.Density_r_cutoff = {"Value": 4.5, "Unit": "[M]"}
Sim_config.disk_model.Temperature_r_cutoff = {"Value": 4.5, "Unit": "[M]"}

Sim_config.disk_model.Density_r_0     = {"Value": 4.5, "Unit": "[M]"}
Sim_config.disk_model.Temperature_r_0 = {"Value": 4.5, "Unit": "[M]"}

Sim_config.disk_model.Opening_angle = {"Value": 0.1, "Unit": "[tan(angle)]"}

# ================================================== Hotspot ================================================== #

Sim_config.hotspot_model.Density_scale_factor = {"Value": 0, "Unit": "[g / cm^3]"}
Sim_config.hotspot_model.Temperature_scale_factor = {"Value": 9.03e10, "Unit": "[K]"}
Sim_config.hotspot_model.Ensamble_type = {"Value": "Kappa", "Unit": "[-]"}
Sim_config.hotspot_model.Magnetization = {"Value": 0.01, "Unit": "[-]"}
Sim_config.emission_models.Kappa = {"Value": 5, "Unit": "[-]"}
Sim_config.hotspot_model.Distance = {"Value": 9, "Unit": "[M]"}
Sim_config.hotspot_model.Velocity_profile = {"Value":"Keplarian", "Unit": "[-]"}

Sim_config.hotspot_model.Temperature_profile = {"Value": "Sphere", "Unit": "[-]"}
Sim_config.hotspot_model.Density_profile = {"Value": "Sphere", "Unit": "[-]"}
Sim_config.hotspot_model.Radius = {"Value": 1, "Unit": "[M]"}
# ================================================== Novikov - Thorne Disk ================================================== #

from Support_functions import Spacetimes

Wormhole_class = Spacetimes.Wormhole(r_throat = 1, parameter = 2)

Sim_config.NT_model_params.Evaluate_NT_disk = {"Value": 1, "Unit": "[-]"}
Sim_config.NT_model_params.r_in = {"Value": 6, "Unit": "[M]"}
Sim_config.NT_model_params.r_out = {"Value": 30, "Unit": "[M]"}
# ================================================== Integrator ================================================== #

# Sim_config.integrator.step_controller_I_gain = {"Value": 0.18, "Unit": "[-]"}
Sim_config.integrator.RK45_accuracy = {"Value": 1e-12, "Unit": "[-]"}

Sim_config.file_manager.Sim_mode_2_input_file_path = "C:/Users/Valur/Documents/University stuff/General Relativity/Polarization/Schwarzschild_Impact_parameters/Direct_image/geodesic_data_20_deg_Sch_r6_500_photons.txt"

Sim_config.generate_simulation_input()


import subprocess
filename = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\Utilities\\FILE.xml"
args = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\x64\\Release\\Mjolnir_GRRT.exe -in " + filename

# for i in range(19):

subprocess.call(args, shell=True)

print("kek")