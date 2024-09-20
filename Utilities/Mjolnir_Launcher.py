class Simulation_configurator:

    def __init__(self):

        self._configure_integrator_settings()
        self._configure_observer()
        self._configure_metric()
        self._configure_NT_model()
        self._configure_GOT_model()
        self._configure_shader_paths()

    def _configure_integrator_settings(self, init_stepsize: float = 1e-5,
                                             integral_accuracy: float = 1e-6,
                                             RK45_accuracy: float = 1e-12,
                                             safety_factor_1: float = 0.5,
                                             safety_factor_2: float = 1e-25,
                                             step_controller_I_gain: float = 0.117,
                                             step_controller_P_gain: float = -0.042,
                                             step_controller_D_gain: float = 0.02,
                                             max_integration_count: int = 1e7):

        self.init_stepsize = init_stepsize
        self.integral_accuracy = integral_accuracy
        self.RK45_accuracy = RK45_accuracy
        self.safety_factor_1 = safety_factor_1
        self.safety_factor_2 = safety_factor_2
        self.step_controller_I_gain = step_controller_I_gain
        self.step_controller_P_gain = step_controller_P_gain
        self.step_controller_D_gain = step_controller_D_gain
        self.max_integration_count  = max_integration_count

    def _configure_observer(self, r_obs: float = 1e4,
                                  theta_obs: float = 2.792,
                                  phi_obs: float = 0.0,
                                  obs_cam_rotation_angle: float = -2.00675,
                                  image_y_min: float = -15,
                                  image_y_max: float = 15,
                                  image_x_min: float = -15,
                                  image_x_max: float = 15,
                                  resolution_y: int = 2048,
                                  resolution_x: int = 2048):
        
        self.r_obs = r_obs
        self.theta_obs = theta_obs
        self.phi_obs = phi_obs
        self.obs_cam_rotation_angle = obs_cam_rotation_angle
        self.image_y_min = image_y_min
        self.image_y_max = image_y_max
        self.image_x_min = image_x_min
        self.image_x_max = image_x_max
        self.resolution_x = resolution_x
        self.resolution_y = resolution_y

    def _configure_metric(self, spin: float = 0.98,
                                wh_redshift: float = 2.0,
                                wh_r_throat: float = 1.0,
                                rbh_param: float = 0.5,
                                jnw_gamma: float = 0.48,
                                gb_gamma: float = 1.15,
                                compactness: float = 1e-4,
                                m_halo: float = 1e4,
                                metric: str = "Kerr"):
        
        self.spin = spin
        self.wh_redshift = wh_redshift
        self.wh_r_throat = wh_r_throat
        self.rbh_param = rbh_param
        self.jnw_gamma = jnw_gamma
        self.gb_gamma = gb_gamma
        self.compactness = compactness
        self.m_halo = m_halo
        self.metric = metric

    def _configure_NT_model(self, r_in: float = 6,
                                  r_out: float = 35,
                                  evaluate_NT_disk: bool = True):
        
        self.NT_r_in = r_in
        self.NT_r_out = r_out
        self.evaluate_NT_disk = evaluate_NT_disk

    def _configure_GOT_model(self, emission_model: str = "Synchotron Exact",
                                   density_model: str = "Power Law",
                                   opening_angle: float = 0.1,
                                   cutoff_scale: float = 0.4,
                                   r_cutoff: float = 5.0,
                                   r_0: float = 5.0,
                                   disk_magnetization: float = 0.01,
                                   mag_field_geometry: list = [0.87, 0, 0.5],
                                   electron_density_scale: float = 5e5,
                                   electron_T_scale: float = 5.1e10,
                                   hotspot_rel_scale: float = 0,
                                   hotspot_scale: float = 1,
                                   hotspot_coord: list = [8.0, 0, 0],                             
                                   emission_power_law: float = 0.0,
                                   source_f_power_law: float = 2.5,
                                   absorbtion_coeff: float = 1.5,
                                   emission_scale_coeff: float = 3e-18,
                                   observation_frequency: float = 230e9,
                                   average_emission_pitch_angle: bool = True,
                                   include_polarization: bool = False):
        
        self.emission_model = emission_model
        self.density_model = density_model
        self.opening_angle = opening_angle
        self.cutoff_scale = cutoff_scale
        self.r_cutoff = r_cutoff
        self.r_0 = r_0
        self.disk_magnetization = disk_magnetization
        self.mag_field_geometry = mag_field_geometry
        self.electron_density_scale = electron_density_scale
        self.electron_T_scale = electron_T_scale
        self.hotspot_rel_scale = hotspot_rel_scale
        self.hotspot_scale = hotspot_scale
        self.hotspot_coord = hotspot_coord
        self.emission_power_law = emission_power_law
        self.source_f_power_law = source_f_power_law
        self.absorbtion_coeff = absorbtion_coeff
        self.emission_scale_coeff = emission_scale_coeff
        self.average_emission_pitch_angle = average_emission_pitch_angle
        self.include_polarization = include_polarization
        self.observation_frequency = observation_frequency

    def _configure_shader_paths(self, vert_shader_path: str = "C:\\Users\\Valur\\Documents\\Repos\\Gravitational_Lenser\\Libraries\\shaders\\default.vert",
                                      frag_shader_path: str = "C:\\Users\\Valur\\Documents\\Repos\\Gravitational_Lenser\\Libraries\\shaders\\default.frag"):
        
        self.vert_shader_path = vert_shader_path
        self.frag_shader_path = frag_shader_path

    def create_input_file(self):





        pass