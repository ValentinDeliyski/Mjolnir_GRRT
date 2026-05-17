#pragma once
#include "Model_Structs.h"
#include "Enumerations.h"
#include <string>

/* Declaring these here, so I dont have to include their headers,
   which hold a whole bunch of function declarations that are not needed in this header. */

class Spacetime_Base_Class;
class Emission_models_class;
class Novikov_Thorne_Model_class;
class Observer_class;
class File_manager_class;
class Geodesic_Integrator_class;
struct Disk_model_type;

struct Hotspot_profile_parameters_type {

    /* ======== Gaussian parameters ======== */

    double Gaussian_variable{};
    double Gaussian_spread{};
    double Gaussian_mean{};

    /* ======== Spherical parameters ======== */

    double Distance_from_sphere_center{};
    double Sphere_radius{};

    /* ======== Power law parameters ======== */

    double Power_law_variable{};
    double Power_law_power{};
    double Power_law_scale{};

};

struct Disk_profile_parameters_type {

    /*! @brief The power law variable ~ (power_law_scale / radial_coordinate)^power */
    double radial_coordinate{};

    /* */

    double theta_coordinate{};

    /*! @brief The power law scale ~ (power_law_scale / radial_coordinate)^power */
    double power_law_scale{};

    /*! @brief The power law exponent ~ (power_law_scale / radial_coordinate)^power */
    double power{};

    /*! @brief The variable of the vertical gaussian ~ exp(-(gauss_var - gauss_mean)^2 / gauss_std^2 / 2) */
    double gaussian_variable{};

    /*! @brief The mean of the vertical gaussian ~ exp(-(gauss_var - gauss_mean)^2 / gauss_std^2 / 2) */
    double gaussian_mean{};

    /*! @brief The standard deviation of the vertical gaussian ~ exp(-(gauss_var - gauss_mean)^2 / gauss_std^2 / 2) */
    double gaussian_std{};

    /*! @brief The cutoff radius for the radial gaussian ~ exp(-(r - cutoff_radius)^2 / cutoff_scale^2 / 2) */
    double cutoff_radius{};

    /*! @brief The cutoff scale for the radial gaussian ~ exp(-(r - cutoff_radius)^2 / cutoff_scale^2 / 2) */
    double cutoff_scale{};

};

struct Disk_model_parameters_type {

    bool Enable_flag{};

    double r_ISCO;

    double Ang_momentum_below_ISCO;

    double Ang_momentum_exponent;

    /*! @brief Specifies the used model of the disk. */
    Disk_model_enums e_Disk_model{};

    /*! @brief Specifies the statistical ensamble of the disk. */
    Ensamble_enums Ensamble_type{};

    /*! @brief Specifies the velocity profile of the disk. */
    Velocity_enums Velocity_profile_type{};

    /*! @brief Specifies the relative density at which we start evaluating the emission of the disk. */
    double Threshold_relative_density{};

    double Max_disk_density{};

    /* ========= The disk profile parameters ========= */

    Common_RIAF_params_type Common_RIAF_params{};

    Colab_test_1_params_type Colab_test_1_params{};

    Novikov_Thorne_params_type Novikov_Thorne_params{};

    Numerical_disk_params_type Numerical_disk_params{};

};

struct Magnetic_fields_type {

    /*! @brief The magnetic field 4-vector in the plasma frame in units of [G]. */
    double B_field_plasma_frame[4]{};

    /*! @brief The magnetic field 4-vector in the coordinate frame in units of [G]. */
    double B_field_eulerian_frame[4]{};

    /*! @brief Unit vector that specifies the direction of the magnetic field in the plasma frame. */
    double Mag_field_geometry_vector[3]{};

    /*! @brief The magnitude of the magnetic field in the plasma frame in units of [G]. */
    double B_field_plasma_frame_norm{};

};

struct Hotspot_model_parameters_type {

    bool Enable_flag{};

    /*! @brief Specifies the statistical ensamble of the hotspot. */
    Ensamble_enums Ensamble_type{};

    /*! @brief Specifies the density profile of the hotspot. The current supported profiles are:
       - Gaussian
       - Sphere with a constant Radius  */
    Profile_enums Density_profile_type{};

    /*! @brief Specifies the temperature profile of the hotspot. The current supported profiles are:
       - Gaussian
       - Sphere with a constant Radius */
    Profile_enums Temperature_profile_type{};

    /*! @brief Specifies the velocity profile of the hotspot. */
    Velocity_enums Velocity_profile_type{};

    /*! @brief The hotspot initial position, specified as [Time, Distance, Polar Angle, Azimuth Angle]. */
    double Init_Position[4]{};

    /*! @brief Struct that holds all the hotspot model parameters. */
    Hotspot_model_params_type Profile_params{};

    /*! @brief The peak density value in [g / cm^3]. */
    double Electron_density_scale{};

    /*! @brief The peak temperature value in [K]. */
    double Electron_temperature_scale{};

    /*! @brief Specifies the direction of the magnetic field. */
    Magnetic_field_geometry_enums e_Mag_field_geometry{};

    /*! @brief Specifies how the magnitude of the magnetic field is calculated. */
    Magnetic_field_magnitude_enums e_Mag_field_magnitude_profile{};

    /* The hotspot magnetization value [-]. */
    double Magnetization{};

    /*! @brief The constant magnetic field geometry in the Eularian frame.
       The components are specified as [B_r, B_theta, B_phi].
       This vector gets normalized with the metric before use. */
    double Mag_field_geometry[3]{};

    /*! @brief The density at which we start evaluating the emission. */
    double Threshold_relative_density{};

    /*! @brief This is the overall scale factor for the power-law magnetic field profile, which is given by Mag_field_magnitude_scale * pow(Mag_field_radial_scale / State_Vector[e_r], Mag_field_power) */
    double Mag_field_B_0{};

    /*! @brief This determines the power of the power-law magnetic field profile, which is given by Mag_field_magnitude_scale * pow(Mag_field_radial_scale / State_Vector[e_r], Mag_field_power) */
    double Mag_field_power{};

    /*! @brief This affects the power-law magnetic field profile, which is given by Mag_field_magnitude_scale * pow(Mag_field_radial_scale / State_Vector[e_r], Mag_field_power) */
    double Mag_field_r_0{};

};

struct Emission_medium_state_type {

    /*! @brief The number density of the emission medium at the current point in units of [1 / cm^3]. */
    double Density{};

    /*! @brief The temperature of the emission medium at the current point in units of [K]. */
    double Temperature{};

    /*! @brief Pointer to the plasma velocity array at the current point in geometric units and local coordinates. */
    double Plasma_Velocity[4];

    /*! @brief Struct that holds the magnetic field in local coordinates. */
    Magnetic_fields_type Magnetic_fields{};

    /*! @brief Enum that specifies which emission medium is being considered. */
    Ensamble_enums Ensamble_type{};

};

struct Hotspot_position_type {

    /*! @brief The distance to the hotspot center in geometric units. */
    double Distance{};

    /*! @brief The polar angle to the hotspot center in units of [Rad]. */
    double Inclination{};

    /*! @brief The azimuthal angle to the hotspot center in units of [Rad].*/
    double Azimuth{};

    /*! @brief The X coordinate of the hotspot center in geometric units. */
    double x{};

    /*! @brief The Y coordinate of the hotspot center in geometric units. */
    double y{};

    /*! @brief The Z coordinate of the hotspot center in geometric units. */
    double z{};

};

struct Emission_model_parameters_type {

    /* ================= Thermal Synchrotron Model ================== */
    /* It is fully determined by the electron density and temperature */

    /* ============= Phenomenological Synchrotron Model ============= */

     /*! @brief The emission is scales linearly with this parameter. It has units of [erg / s / sr / Hz]. */
    double Phenomenological_emission_coeff{};

    /*! @brief The absorbtion scales linearly with this parameter. It has units of [cm^-1]. */
    double Phenomenological_absorbtion_coeff{};

    /*! @brief The emission scales as pow(redshift, Phenomenological_emission_power_law). */
    double Phenomenological_emission_power_law{};

    /*! @brief The source function scales as pow(redshift, Phenomenological_source_f_power_law). */
    double Phenomenological_source_f_power_law{};

    /* ================== Kappa Synchrotron Model ================== */

    /*! @brief The free dimentionless parameter for the kappa distribution. */
    double Kappa{};

    /* ==================== Debug Emission model ==================== */

    /*! @brief Constant value used for testing in the e_Debug_constant_functions emission model. */
    double Debug_j_I_value{};

    /*! @brief Constant value used for testing in the e_Debug_constant_functions emission model. */
    double Debug_j_Q_value{};

    /*! @brief Constant value used for testing in the e_Debug_constant_functions emission model. */
    double Debug_j_U_value{};

    /*! @brief Constant value used for testing in the e_Debug_constant_functions emission model. */
    double Debug_j_V_value{};

    /*! @brief Constant value used for testing in the e_Debug_constant_functions emission model. */
    double Debug_alpha_I_value{};

    /*! @brief Constant value used for testing in the e_Debug_constant_functions emission model. */
    double Debug_alpha_Q_value{};

    /*! @brief Constant value used for testing in the e_Debug_constant_functions emission model. */
    double Debug_alpha_U_value{};

    /*! @brief Constant value used for testing in the e_Debug_constant_functions emission model. */
    double Debug_alpha_V_value{};

    /*! @brief Constant value used for testing in the e_Debug_constant_functions emission model. */
    double Debug_rho_I_value{};

    /*! @brief Constant value used for testing in the e_Debug_constant_functions emission model. */
    double Debug_rho_Q_value{};

    /*! @brief Constant value used for testing in the e_Debug_constant_functions emission model. */
    double Debug_rho_U_value{};

    /*! @brief Constant value used for testing in the e_Debug_constant_functions emission model. */
    double Debug_rho_V_value{};

};

struct Numerical_metric_potentials_type {

    double exp_2F_0{};
    double exp_2F_1{};
    double exp_2F_2{};
    double W{};

};

struct Numerical_metric_params_type {

    Spline_selection_enums e_Spline_type{};

    /*! @brief The ADM mass. */
    double M_ADM{};

    /*! @brief The horizon radius. NOTE: This is NOT normalized to the mass. */ 
    double Horizon_radius{};

    /*! @brief The ADM rotation parameter. NOTE: This is NOT normalized to the mass. */ 
    double a_ADM{};

    /* ============ Pointers to the numerical metric potentials control vector arrays ============ = */

    double* F_0_control_vector{};
    double* F_1_control_vector{};
    double* F_2_control_vector{};
    double* W_control_vector{};
    int Control_vector_size{};

    double* Raw_F_0_data{};
    double* Raw_F_1_data{};
    double* Raw_F_2_data{};
    double* Raw_W_data{};

    /* ============ Pointers to the coordinate grid, its step and control vector arrays ============ = */

    double* Compactified_radial_grid{};
    double* Compactified_radial_grid_steps{};
    double* Compactified_radial_grid_control_vector{};
    long long Radial_grid_size{};

    double* Theta_grid{};
    double* Theta_grid_steps{};
    double* Theta_grid_control_vector{};
    long long Theta_grid_size{};

    /*! @brief Enum that decides which anzatz to use for the numerical metric line element. */
    Numerical_Anzatz_enums e_Anzatz{};

    /* @brief String with the full path to the metric spline XML. */
    std::string Metric_file_path{};

};

struct Spline_arguments_type {

    /*! @brief Real value in the range [0, 1] that interpolates between radial grid points. */
    double Radial_natural_param{};

    /*! @brief Real value in the range [0, 1] that interpolates between theta grid points. */
    double Theta_natural_param{};

    /*! @brief The index of the desired grid point in the Compactified_radial_grid array. */
    long long Radial_idx{};

    /*! @brief The index of the desired grid point in the Theta_grid array. */
    long long Theta_idx{};

    /*! @brief The interpolated radial coordinate. */
    double r_coord{};

    /*! @brief The interpolated compactified radial coordinate. */
    double r_coord_compactified{};

    double Theta_coord;

};

struct Delta_coeffs_type {

    double a_coeff{};

    double b_coeff{};

    double c_coeff{};

    double d_coeff{};

    double e_coeff{};

    double f_coeff{};

};

struct Metric_parameters_type {

    /*! @brief Enum that specifies the active spacetime. */
    Spacetime_enums e_Spacetime{};

    /* ============ Wormhole Specific Parameters ============ = */

     /*! @brief The Wormhole lapse function is given by exp(-M / r - Redshift_param * (M / r)^2). */
    double Redshift_Parameter{};

    /*! @brief The Wormhole g_rr function is given by 1. / (1 - R_throat / r). */
    double R_throat{};

    /*! @brief Boolean that decides weather to allow photons to cross the Worhmole throat. */
    bool Stop_At_Throat{};

    /* ============ Janis-Newman-Winicour Specific Parameters ============ = */

    double JNW_Gamma_Parameter{};

    /* ============ Gauss-Bonnet Specific Parameters ============ = */

    double GB_Gamma_Parameter{};

    /* ============ Regular Black Hole Specific Parameters ============ = */

    double RBH_Parameter{};

    /* ============ Black Hole w Dark Matter Halo Specific Parameters ============ = */

    double Compactness{};
    double Halo_Mass{};

    Numerical_metric_params_type Numerical_metric_params{};

    /* ============ Generic Parameters ============ = */

    double Mass{};

    double Spin{}; // Only affects Kerr and the Wormhole

    /*! @brief Radial distance after which we stop the integration. */
    double Scattering_radius{};

    /*! @brief The minimum distance to special surfaces (event horizons, WH throat, ...), after which we stop the integration. */
    double Min_distance_to_singular_point{};

};

struct Precomputed_e_pitch_angles_type {

    double* sin_electron_pitch_angles{};
    double* cos_electron_pitch_angles{};

    // === Used in the thermal synchotron emission functions === //

    /*! @brief 1. / sqrt(sin) */
    double* one_over_sqrt_sin{};

    /*! @brief 1. / cbrt(sin) */
    double* one_over_cbrt_sin{};

    // === Used in the thermal synchotron Faradey functions === //

     /*! @brief 1. / pow(sin, 0.5175) */
    double* one_over_sin_to_0_p_5175{};

    /*! @brief 1. / pow(sin, 0.6) */
    double* one_over_sin_to_0_p_6{};

    /*! @brief 1. / pow(sin, 0.7515) */
    double* one_over_sin_to_0_p_7515{};

    // === Used in the kappa synchotron emission functions === //

    /*! @brief 1. / pow(sin, 7. / 20) */
    double* one_over_sin_to_7_over_20{};

};

struct Thermal_transfer_f_arguments_type {

    /*! @brief The dimensionless parameter that the thermal fit functions depend on. */
    double X{};

    /*! @brief sqrt(X). */
    double sqrt_X{};

    /*! @brief cbrt(X). */
    double cbrt_X{};

    /*! @brief pow(X, 0.5175). */
    double X_to_0_p_5175{};

    /*! @brief pow(X, 0.6). */
    double X_to_0_p_6{};

    /*! @brief pow(X, 0.7515). */
    double X_to_0_p_7515{};

    /*! @brief The dimensionless electron temperature at the current point. */
    double T_electron_dim{};

    /*! @brief pow(X, 24. / 25). */
    double T_electron_dim_to_24_25{};

    /*! @brief Sin of the emission angle. */
    double sin_pitch_angle{};

    /*! @brief Cos of the emission angle. */
    double cos_pitch_angle{};

    /*! @brief The current photon frequency in units of [Hz]. */
    double frequency{};

};

struct Kappa_transfer_f_arguments_type {

    /*! @brief The dimensionless parameter that the kappa fit functions depend on. */
    double X{};

    /*! @brief sqrt(X). */
    double sqrt_X{};

    /*! @brief cbrt(X). */
    double cbrt_X{};

    /*! @brief pow(X, 7. / 20). */
    double X_to_7_over_20{};

    /*! @brief The dimensionless electron temperature at the current point. */
    double T_electron_dim{};

    /*! @brief Sin of the emission angle. */
    double sin_emission_angle{};

    /*! @brief Cos of the emission angle. */
    double cos_emission_angle{};

    /*! @brief The kappa parameter. */
    double kappa{};

    /*! @brief The current photon frequency in units of [Hz]. */
    double frequency{};

};

struct Phenomenological_transfer_f_arguments_type {

    /*! @brief The redshift at the current point. */
    double redshift{};

    /*! @brief The cyclotron frequency at the current point in units of [Hz]. */
    double f_cyclo{};

    /*! @brief The current photon frequency in units of [Hz]. */
    double frequency{};

};

struct Transfer_functions_type {

    /*! @brief Array that holds the current emission functions for all polarizations. */
    double Emission_functions[e_Stokes_param_num]{};

    /*! @brief Array that holds the current Faradey functions for all polarizations. */
    double Faradey_functions[e_Stokes_param_num]{};

    /*! @brief Array that holds the current absorbtion functions for all polarizations. */
    double Absorbtion_functions[e_Stokes_param_num]{};

};

struct Metric_type {

    /*! @brief The metric tensor. */
    double Metric[4][4]{};

    /*! @brief The metric lapse function. */
    double Lapse_function{};

    /*! @brief The metric shift function. */
    double Shift_function{};

};

struct Step_Controller_parameters_type {

    /*! @brief The integral gain of the PID controller for the explicit Runge-Kutta integrator. */
    double RK_PID_gain_I{};

    /*! @brief The proportional gain of the PID controller for the explicit Runge-Kutta integrator. */
    double RK_PID_gain_P{};

    /*! @brief The derivative gain of the PID controller for the explicit Runge-Kutta integrator. */
    double RK_PID_gain_D{};

    /*! @brief The k1 gain of the Gustafsson controller for the explicit Runge-Kutta integrator. */
    double RK_Gustafsson_k1{};

    /*! @brief The k2 gain of the Gustafsson controller for the explicit Runge-Kutta integrator. */
    double RK_Gustafsson_k2{};

    /*! @brief The integral gain of the PID controller for the ESDIRK54 integrator. */
    double ESDIRK54_PID_gain_I{};

    /*! @brief The proportional gain of the PID controller for the ESDIRK54 integrator. */
    double ESDIRK54_PID_gain_P{};

    /*! @brief The derivative gain of the PID controller for the ESDIRK54 integrator. */
    double ESDIRK54_PID_gain_D{};

    /*! @brief The k1 gain of the Gustafsson controller for the ESDIRK54 integrator. */
    double ESDIRK54_Gustafsson_k1{};

    /*! @brief The k2 gain of the Gustafsson controller for the ESDIRK54 integrator. */
    double ESDIRK54_Gustafsson_k2{};

    /*! @brief The maximum allowed relative step increase. */
    double Max_rel_step_increase{};

    /*! @brief The minimum allowed relative step increase. */
    double Min_rel_step_increase{};

    /*! @brief The initial stepsize. */
    double Init_stepzie{};

    /*! @brief The Maximum allowed stepsize, when using the adaptive integrator. Also the constant value for the fixed step. */
    double Max_upper_stepsize{};

    double Min_upper_stepsize{};

    /*! @brief A multiplicative factor forr the integration step in the range (0, 1] that makes the integrator more stable. */
    double Safety_1{};

    /*! @brief A additive factor <<1 that ensures no "division by 0" problems occur when calculating the new stepsizes. */
    double Safety_2{};

    /*! @brief Enum that decides weather to use an adaptive step (with the selected step controller), or a fixed one, determined by the initial step. */
    bool Use_adaptive_step{};

    /*! @brief Enum that decides which step controller to use. */
    Step_controller_type_enums Controller_type{};

    /*! @brief The adaptive explicit Kunge-Kutta integrator absolute error threshold parameter. */
    double RK_abs_accuracy{};

    /*! @brief The adaptive explicit Kunge-Kutta integrator relative error threshold parameter. */
    double RK_rel_accuracy{};

    /*! @brief The adaptive ESDIRK5(4) absolute error threshold parameter. */
    double ESDIRK54_abs_accuracy{};

    /*! @brief The adaptive ESDIRK5(4) relative error threshold parameter. */
    double ESDIRK54_rel_accuracy{};

    double Dist_to_Observer{};

    double Dist_at_min_upper_stepsize{};

    double Max_step_b_coeff{};

};

struct Integrator_parameters_type {

    Step_Controller_parameters_type Geodesic_Step_Controller_Params;

    /*! @brief The Maximum allowed affine parameter value, before terminating the integration.
        NOTE: This is taken by absolute value. */
    double Max_affine_param{};

    /*! @brief The maximum allowed (accepted) integration steps. */
    int Max_integration_count{};

    /*! @brief Enum that selects which radiative transfer integrator to use. */
    Integrator_enums e_Radiative_transfer_integrator{};

    /*! @brief Enum that selects which radiative transfer integrator to use. */
    Integrator_enums e_Parallel_transport_integrator{};

    /*! @brief Enum that selects which geodesic integrator to use by default. When that one fails to integrate a given geodesic, we switch to ESDIRK54. */
    Integrator_enums e_Default_geodesic_integrator{};

};

struct von_Zeipel_cylinder_condition_wrapper_struct {

    Disk_model_type* Disk_model{};

    Metric_type Metric{};

};

struct RHS_wrapper_struct {

    /*! @brief Pointer to the inegrator class instance. This exists so I can call member functions from the RHS wrapper in the run_ESDIRK54 function.  */
    Geodesic_Integrator_class* Integrator{};

    /*! @brief Pointer to the current iteration number. This is the gsl way of passing parameters to functions.  */
    void* p_Iteration_number{};

};

struct Observer_parameters_type {

    /*! @brief The moment of observation, relative to the central object dynamics, in geometric units. */
    double init_time{};

    /*! @brief The distance to the observer in geometric units. */
    double distance{};

    /*! @brief The observer's polar angle in units of [Rad]. */
    double inclination{};

    /*! @brief The azimuthal angle of the observer in units of [Rad]. */
    double azimuth{};

    /*! @brief The minimum linear X coordinate of the observer's image. */
    double x_min{};

    /*! @brief The maximum linear X coordinate of the observer's image. */
    double x_max{};

    /*! @brief The minimum linear Y coordinate of the observer's image. */
    double y_min{};

    /*! @brief The maximum linear Y coordinate of the observer's image. */
    double y_max{};

    /*! @brief The minimum angular X coordinate of the observer's image. */
    double x_angle_min{};

    /*! @brief The maximum angular X coordinate of the observer's image. */
    double x_angle_max{};

    /*! @brief The minimum angular Y coordinate of the observer's image. */
    double y_angle_min{};

    /*! @brief The maximum angular Y coordinate of the observer's image. */
    double y_angle_max{};

    /*! @brief Flag that determines which coordinates (linear or angular) to use for setting up the viewing window. */
    bool Use_angular_coords{};

    /*! @brief The observation frequency in units of [Hz]. */
    double obs_frequency{};

    /*! @brief The observer's rotation angle along the line of sight in units of [Rad].*/
    double cam_rotation_angle{};

    /*! @brief The resolution along the X axis of the observer's image. */
    int resolution_x{};

    /*! @brief The resolution along the Y axis of the observer's image. */
    int resolution_y{};

    /*! @brief Flag that decides weather to include polarization in the radiative transfer. */
    bool include_polarization{};

};

struct File_manager_parameters {

    /*! @brief The full path for the input file for simulation mode 2. */
    std::string Sim_mode_1_imput_path{};

    /*! @brief The full path for the output file directory. The actual simulation results get stored in the
       ...Output_file_directory\\Simulation_name folder. */
    std::string Output_file_directory{};

    /*! @brief The common part of the names for all output files. These get _nX appended to them, where X is the image order. */
    std::string Common_file_names{};

    /*! @brief The full path to the vertex shader file. */
    std::string Vert_shader_path{};

    /*! @brief The full path to the fragment shader file. */
    std::string Frag_shader_path{};

    /*! @brief THe simulation name. */
    std::string Simulation_name{};

    /*! @brief Flag that decieds weather to trncate the output files or not. */
    bool Truncate_files{};

};

struct Initial_conditions_type {

    /*! @brief Struct that holds the file manager parameters. */
    File_manager_parameters File_manager_params{};

    /*! @brief Struct that holds the disk model parameters. */
    Disk_model_parameters_type Disk_params{};

    /*! @brief Struct that holds the inital metric (evaluated at the observer). */
    Metric_type Init_metric{};

    /*! @brief Struct that holds the hotspot model parameters. */
    Hotspot_model_parameters_type Hotspot_params{};

    /*! @brief Struct that holds the integrator parameters. */
    Integrator_parameters_type Integrator_params{};

    /*! @brief Struct that holds the observer parameters. */
    Observer_parameters_type Observer_params{};

    /*! @brief Struct that holds the metric parameters. */
    Metric_parameters_type Metric_parameters{};

    /*! @brief Struct that holds the emission model parameters. */
    Emission_model_parameters_type Emission_params{};

    Order_counging_enums e_Order_counting_scheme{};

    /*! @brief Initial condition for simulation mode 3. The image X coordiante is used to compute the initial azimuthal momentum. */
    double Sim_mode_2_X_init{};

    /*! @brief Initial condition for simulation mode 3. The image Y coodirnate is used to compute the initial polar momentum. */
    double Sim_mode_2_Y_init{};

    /*! @brief The initial photon momentum. */
    double Init_Momentum[4]{};

    /*! @brief The central object mass in units of [M_sun]. */
    double central_object_mass{};

    /*! @brief Enum that specifies the currently active simulation mode. */
    Simulation_mode_enums Simulation_mode{};

    /*! @brief The number of metric parameter values that will be used during simulation mode 2.
       I usually use this mode to do "parameter sweeps", so this comes in handy when parsing the results. */
    int Sim_mode_1_param_value_number{};

    /*! @brief The number of samples of the electron pitch angle to use when averaging the emission.
       NOTE: This only has an effect if Average_electron_pitch_angle is set to "ture". */
    int Emission_pitch_angle_samples_to_average{};

    /*! @brief Boolean flag that decides weather to print the simulation metadata to the console. */
    bool Print_to_console{};

    /*! @brief Boolean flag that decides weather to average the emission over the electron pitch angle. */
    bool Average_electron_pitch_angle{};

    /*! @brief The maximum image order for which the radiative transfer will be evaluated.
       NOTE: For image orders above this, the ray is still propagated, along with the parallel transport
       of the polarization vector */
    int Max_order{};

    /*! @brief The minimum image order for which the radiative transfer will be evaluated.
       NOTE: For image orders blow this, the ray is still propagated, along with the parallel transport
       of the polarization vector (though the stokes vector will be zero, so the "if" statement in which
       the parallel transport is placed won't pass, so really only the ray will be propagated). */
    int Min_order{};

};

struct Simulation_Context_type {

    /*! @brief Pointer to the struct that holds the initial conditions. */
    Initial_conditions_type* p_Init_Conditions{};

    /*! @brief Pointer to the class that holds all the spacetime related functions. */
    Spacetime_Base_Class* p_Spacetime{};

    /*! @brief Pointer to the class that holds all the observer related functions. */
    Observer_class* p_Observer{};

    /*! @brief Pointer to the class that holds all the emission medium related functions. */
    Emission_models_class* p_Emission_Model{};

    /*! @brief Pointer to the class that holds all the Novikov-Thorne related functions. */
    Novikov_Thorne_Model_class* p_NT_model{};

};

struct Ray_log_type {

    /* Pointer to the arrays that hold the Stokes vector along the photon trajectory. */
    double* Ray_emission_log[e_Stokes_param_num];

    double* Ray_polarization_log[2];

    /* Pointer to the array that holds the entire photon trajectory in local coordinates. */
    double* Ray_path_log_local;

    /* Pointer to the array that holds the entire photon trajectory in global coordinates. */
    double* Ray_path_log_global;

    size_t Log_offet_at_disk_edge;

    /* Int that specifies where in the log to write.
       This exists for the sole purpose of minimizing the number of arguments in the functions that write to the photon log. */
    size_t Log_offset;

    /* The length of the photon log. */
    size_t Log_length;

};

struct Adaptive_RK_Integrator_debug_type {

    /*! @brief Pointer to the log of the state error. Used only in sim mode Make_geodesic_log. */
    double* State_error_history{};

    /*! @brief Pointer to the log of the number of rejected steps. Used only in sim mode Make_geodesic_log. */
    double* N_steps_rejected{};

};

struct Polarization_debug_type {

    double* Polarization_fourvec[4];

    double* PW_constant[2];

};

struct Results_type {

    /*! @brief The struct that tholds the metric parameters. */
    Metric_parameters_type Metric_parameters{};

    /*! @brief The struct that holds the photon log. */
    Ray_log_type Ray_log_struct{};

    /*! @brief The struct that holds the adaptive integrator debug parameters log. */
    Adaptive_RK_Integrator_debug_type RK_integrator_debug_log{};

    Polarization_debug_type Polarization_debug_log{};

    /*! @brief Array that holds the integrated intensity for each polarization component. */
    double Intensity[e_Stokes_param_num]{};

    /*! @brief Array that holds the photon state vector when the integration is terminated. */
    double Final_State_Vector[e_Full_state_size]{};

    /*! @brief Array that holds the state vector at the emission point of the Novikov-Thorne disk.
       This exists for use in simulation mode 2. */
    double Thin_Disk_State_Vector[e_Dynamic_state_size]{};

    /*! @brief Array that holds the phi and theta coordinates of the ray when it reaches the scattering radius. */
    double Celestial_sphere_crossing_coords[4]{};

    /*! @brief Boolean that keeps track weather the Novikov-Thorne disk was found by the ray (going from the observer backwards).
        This exists so the disk image orders overlap correctly in the final image. */
    bool NT_Disk_found{};

    /*! @brief Array that holds the Novikov-Thorne disk flux. */
    double Flux_NT{};

    /*! @brief Array that holds the Novikov-Thorne disk polariation vector (transported to the observer) in the ZAMO frame. */
    double Projected_polarization_vector[2]{};

    /*! @brief Array that holds the Novikov-Thorne disk redshift. */
    double Redshift_NT{};

    /*! @brief Placeholder for an array that will hold the integrated optical depth. */
    double Optical_Depth{};

    /*! @brief Array that holds the coordinates of the image on the observer plane.
        NOTE: These get affected by the cam_rotation_angle parameter of the observer. */
    double Image_Coords[2]{};

};

struct Debug_mode_struct {

    int Array_length;

    double* NT_Flux_integral_array;
    double* NT_Flux_r_coord_array;


};