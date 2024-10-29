from numpy import pi
import numpy as np

import xml.etree.cElementTree as ET
import xml.dom.minidom

from Support_functions.Parsers import Units_class

class Integrator():

    __slots__ = ("init_stepsize",
                 "RK45_accuracy", 
                 "step_controller_safety_factor_1",
                 "step_controller_safety_factor_2",
                 "step_controller_I_gain",
                 "step_controller_P_gain",
                 "step_controller_D_gain", 
                 "max_integration_count")

class Disk_model():

    __slots__ = ("Ensamble_type",
                 "Density_profile",
                 "Temperature_profile",
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
                 "Density_scale_factor",
                 "Temperature_scale_factor",
                 "Mag_field_geometry_X",
                 "Mag_field_geometry_Y",
                 "Mag_field_geometry_Z",
                 "Density_spread",
                 "Temperature_spread",
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

class File_paths():

    __slots__ = ("Sim_mode_2_input_file_path", "Output_file_path", "Common_file_names", "Vert_shader_path", "Frag_shader_path")

class Simulation_configurator:

    __slots__ = ("simulation_name", "integrator", "metric_parameters", "disk_model", "hotspot_model", "observer", "emission_models", "NT_model_params", "file_paths", "Average_emission_pitch_angle")

    def __init__(self, Average_emission_pitch_angle: dict = {"Value": 1, "Unit": "[-]"}):

        self. Average_emission_pitch_angle = Average_emission_pitch_angle

        self._configure_integrator_settings()
        self._configure_observer()
        self._configure_NT_model()
        self._configure_disk_model()
        self._configure_hotspot_model()
        self._configure_file_paths()
        self._configure_metric_parameters()
        self._configure_emission_models()

    def _configure_integrator_settings(self, Init_stepsize: dict = {"Value": 1e-5, "Unit": "[M]"},
                                             RK45_accuracy: dict = {"Value": 1e-12, "Unit": "[-]"},
                                             Safety_factor_1: dict = {"Value": 0.8, "Unit": "[-]"},
                                             Safety_factor_2: dict = {"Value": 1e-25, "Unit": "[-]"},
                                             Step_controller_I_gain: dict = {"Value": 0.117, "Unit": "[-]"},
                                             Step_controller_P_gain: dict = {"Value": -0.042, "Unit": "[-]"},
                                             Step_controller_D_gain: dict = {"Value": 0.02, "Unit": "[-]"},
                                             Max_integration_count: dict = {"Value": 1e7, "Unit": "[-]"},):

        self.integrator = Integrator()

        self.integrator.init_stepsize = Init_stepsize
        self.integrator.init_stepsize = Init_stepsize
        self.integrator.RK45_accuracy = RK45_accuracy
        self.integrator.step_controller_safety_factor_1 = Safety_factor_1
        self.integrator.step_controller_safety_factor_2 = Safety_factor_2
        self.integrator.step_controller_I_gain = Step_controller_I_gain
        self.integrator.step_controller_P_gain = Step_controller_P_gain
        self.integrator.step_controller_D_gain = Step_controller_D_gain
        self.integrator.max_integration_count  = Max_integration_count

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
                                       Density_scale_factor: dict = {"Value": 1e6, "Unit": "[g/cm^3]"},
                                       Temperature_scale_factor: dict = {"Value": 1e11, "Unit": "[K]"},
                                       Mag_field_geometry_X: dict = {"Value": 0.5, "Unit": "[-]"},
                                       Mag_field_geometry_Y: dict = {"Value": 0, "Unit": "[-]"},
                                       Mag_field_geometry_Z: dict = {"Value": 0.87, "Unit": "[-]"},
                                       Density_spread: dict = {"Value": 1.0, "Unit": "[M]"},
                                       Temperature_spread: dict = {"Value": 1.0, "Unit": "[M]"},
                                       Distance: dict = {"Value": 8.0, "Unit": "[M]"},
                                       Inclination: dict = {"Value": pi / 2,  "Unit": "[Rad]"},
                                       Azimuth: dict = {"Value": -pi / 2,  "Unit": "[Rad]"},
                                       Magnetization: dict = {"Value": 1.0,  "Unit": "[-]"}):

        self.hotspot_model = Hotspot_model()

        self.hotspot_model.Ensamble_type      = Ensamble_type     
        self.hotspot_model.Density_profile    = Density_profile   
        self.hotspot_model.Temperature_profile = Temperature_profile
        self.hotspot_model.Density_scale_factor = Density_scale_factor     
        self.hotspot_model.Temperature_scale_factor  = Temperature_scale_factor 
        self.hotspot_model.Mag_field_geometry_X = Mag_field_geometry_X
        self.hotspot_model.Mag_field_geometry_Y = Mag_field_geometry_Y
        self.hotspot_model.Mag_field_geometry_Z = Mag_field_geometry_Z
        self.hotspot_model.Density_spread     = Density_spread  
        self.hotspot_model.Temperature_spread = Temperature_spread          
        self.hotspot_model.Distance           = Distance          
        self.hotspot_model.Inclination        = Inclination       
        self.hotspot_model.Azimuth            = Azimuth           
        self.hotspot_model.Magnetization      = Magnetization     

    def _configure_file_paths(self, Vert_shader_path: str = "C:/Users/Valur/Documents/Repos/Gravitational_Lenser/Libraries/shaders/default.vert",
                                    Frag_shader_path: str = "C:/Users/Valur/Documents/Repos/Gravitational_Lenser/Libraries/shaders/default.frag",
                                    Output_file_path: str = "C:/Users/Valur/Documents/Repos/Gravitational_Lenser/Sim_Results",
                                    Common_file_names: str = "",
                                    Sim_mode_2_input_file_path: str = ""):
                            

        self.file_paths = File_paths()

        self.file_paths.Vert_shader_path = Vert_shader_path
        self.file_paths.Frag_shader_path = Frag_shader_path
        self.file_paths.Output_file_path = Output_file_path
        self.file_paths.Common_file_names = Common_file_names
        self.file_paths.Sim_mode_2_input_file_path = Sim_mode_2_input_file_path

    def generate_simulation_input(self):

        Encoding = 'UTF-8'
        XML_root_node = ET.Element("Simulation_Input", {"Simulation_Name": "INSERT_NAME_HERE"})
        ET.SubElement(XML_root_node, "Average_emission_pitch_angle", units = "[-]").text = "{}".format(self.Average_emission_pitch_angle["Value"])

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

        Gaussin_slots = ["Density_spread", "Temperature_spread"]

        Common_slots = [slot for slot in self.hotspot_model.__slots__ if slot not in Gaussin_slots]
        
        Hotspot_subelement = ET.SubElement(XML_root_node, "Hotspot") 
        for Hotspot_attrib_name in Common_slots:
            Hotspot_attrib = getattr(self.hotspot_model, Hotspot_attrib_name)
            ET.SubElement(Hotspot_subelement, Hotspot_attrib_name, units = Hotspot_attrib["Unit"]).text = "{}".format(Hotspot_attrib["Value"])

        # ------------- Gaussian profile subsection
        Gaussian_subelement = ET.SubElement(Hotspot_subelement, "Gaussian_profile") 
        for Hotspot_attrib_name in Gaussin_slots:
            Hotspot_attrib = getattr(self.hotspot_model, Hotspot_attrib_name)
            ET.SubElement(Gaussian_subelement, Hotspot_attrib_name, units = Hotspot_attrib["Unit"]).text = "{}".format(Hotspot_attrib["Value"])

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

        Files_subelement = ET.SubElement(XML_root_node, "File_Paths")
        for Files_attrib_name in self.file_paths.__slots__:
            Files_attrib = getattr(self.file_paths, Files_attrib_name)
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

Sim_config.observer.Distance = {"Value": 1e3, "Unit": "[M]"}
Sim_config.observer.Inclination = {"Value": 60 * pi / 180, "Unit": "[Rad]"}
Sim_config.disk_model.Ensamble_type = {"Value": "Phenomenological", "Unit": "[-]"}
Sim_config.disk_model.Density_profile = {"Value": "Exponential Law", "Unit": "[-]"}
Sim_config.metric_parameters.Spin = {"Value": 0.9, "Unit": "[M]"}
Sim_config.hotspot_model.Density_scale_factor = {"Value": 0, "Unit": "[-]"}
# Sim_config.disk_model.Temperature_profile = {"Value": "Exponential Law", "Unit": "[-]"}

Sim_config.generate_simulation_input()
