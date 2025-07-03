#pragma once
#include "Model_Structs.h"
#include "Enumerations.h"
#include <string>

/* Declaring these here, so I dont have to include their headers,
   which hold a whole bunch of function declarations that are not needed in this header. */

class Spacetime_Base_Class;
class Emission_models_class;
class Page_Thorne_Model_class;
class Observer_class;
class File_manager_class;

struct Hotspot_profile_parameters_type {

    /* ======== Gaussian parameters ======== */

    double Gaussian_variable;
    double Gaussian_spread;
    double Gaussian_mean;

    /* ======== Spherical parameters ======== */

    double Distance_from_sphere_center;
    double Sphere_radius;

    /* ======== Power law parameters ======== */

    double Power_law_variable;
    double Power_law_power;
    double Power_law_scale;

};

struct Disk_profile_parameters_type {

    /*! The power law variable ~ (power_law_scale / radial_coordinate)^power */
    double radial_coordinate;

    /*! The power law scale ~ (power_law_scale / radial_coordinate)^power */
    double power_law_scale;

    /*! The power law exponent ~ (power_law_scale / radial_coordinate)^power */
    double power;

    /*! The variable of the vertical gaussian ~ exp(-(gauss_var - gauss_mean)^2 / gauss_std^2 / 2) */
    double gaussian_variable;

    /*! The mean of the vertical gaussian ~ exp(-(gauss_var - gauss_mean)^2 / gauss_std^2 / 2) */
    double gaussian_mean;

    /*! The standard deviation of the vertical gaussian ~ exp(-(gauss_var - gauss_mean)^2 / gauss_std^2 / 2) */
    double gaussian_std;

    /*! The cutoff radius for the radial gaussian ~ exp(-(r - cutoff_radius)^2 / cutoff_scale^2 / 2) */
    double cutoff_radius;

    /*! The cutoff scale for the radial gaussian ~ exp(-(r - cutoff_radius)^2 / cutoff_scale^2 / 2) */
    double cutoff_scale;

};

struct Disk_model_parameters_type {

    /*! Specifies the used model of the disk. */
    Disk_model_enums e_Disk_model;

    /*! Specifies the statistical ensamble of the disk. */
    Ensamble_enums Ensamble_type; 

    /*! Specifies the velocity profile of the disk. */
    Velocity_enums Velocity_profile_type;  

    /*! Specifies the magnitude of the radial velocity component. I use this to interpolate the circular velocity profile, 
        specified by the "Velocity_profile_type" enum, with a purely radial profile. The range is [0, 1]. */
    double Radial_velocity_fraction;

    /*! The peak density value in [g / cm^3]. */
    double Electron_density_scale;

    /*! The peak temperature value in [K]. */
    double Electron_temperature_scale;

    /*! Specifies the direction of the magnetic field. */
    Magnetic_field_geometry_enums e_Mag_field_geometry;

    /*! Specifies how the magnitude of the magnetic field is calculated. */
    Magnetic_field_magnitude_enums e_Mag_field_magnitude_profile;

    /*! The disk magnetization value [-]. */
    double Magnetization;

    /*! The constant magnetic field geometry in the Eularian frame.
        The components are specified as [B_r, B_theta, B_phi].
        This vector gets normalized with the metric before use. */
    double Mag_field_geometry[3];

    double Mag_field_magnitude_scale;

    double Mag_field_power;

    double Mag_field_radial_scale;

    /*! Specifies the relative density at which we start evaluating the emission of the disk. */
    double Threshold_relative_density;

    /* ========= The density profile parameters ========= */

    Common_RIAF_params_type Common_RIAF_params;

    Colab_test_1_params_type Colab_test_1_params;

    Page_Thorne_params_type Page_Thorne_params;

};

struct Magnetic_fields_type {

    /*! The magnetic field 4-vector in the plasma frame in units of [G]. */
    double B_field_plasma_frame[4];

    /*! The magnetic field 4-vector in the coordinate frame in units of [G]. */
    double B_field_eularian_frame[4];

    /*! Unit vector that specifies the direction of the magnetic field in the plasma frame. */
    double Mag_field_geometry_vector[3];

    /*! The magnitude of the magnetic field in the plasma frame in units of [G]. */
    double B_field_plasma_frame_norm;

    Magnetic_field_geometry_enums e_Mag_field_geometry;

    Magnetic_field_magnitude_enums e_Mag_field_magnitude_profile;

    double Mag_field_magnitude_scale;

    double Mag_field_power;

    double Mag_field_radial_scale;

};

struct Hotspot_model_parameters_type {

    /*! Specifies the statistical ensamble of the hotspot. */
    Ensamble_enums Ensamble_type; 

    /*! Specifies the density profile of the hotspot. The current supported profiles are:
        - Gaussian
        - Sphere with a constant Radius  */
    Profile_enums Density_profile_type; 

    /*! Specifies the temperature profile of the hotspot. The current supported profiles are:
        - Gaussian
        - Sphere with a constant Radius */
    Profile_enums Temperature_profile_type; 

    /*! Specifies the velocity profile of the hotspot. */
    Velocity_enums Velocity_profile_type;  

    /*! Specifies the magnitude of the radial velocity component. I use this to interpolate the circular velocity profile, 
       specified by the "Velocity_profile_type" enum, with a purely radial profile. The range is [0, 1]. */
    double Radial_velocity_fraction;

    /*! The hotspot initial position, specified as [Time, Distance, Polar Angle, Azimuth Angle] */
    double Position[4]; 

    Hotspot_model_params_type Profile_params;

    /*! The peak density value in [g / cm^3]. */
    double Electron_density_scale;     

    /*! The peak temperature value in [K]. */
    double Electron_temperature_scale;

    /*! Specifies the direction of the magnetic field. */
    Magnetic_field_geometry_enums e_Mag_field_geometry;

    /*! Specifies how the magnitude of the magnetic field is calculated. */
    Magnetic_field_magnitude_enums e_Mag_field_magnitude_profile;
    
    /* The hotspot magnetization value [-]. */
    double Magnetization; 

    /*! The constant magnetic field geometry in the Eularian frame.
        The components are specified as [B_r, B_theta, B_phi]. 
        This vector gets normalized with the metric before use. */
    double Mag_field_geometry[3];

    /*! The density at which we start evaluating the emission. */
    double Threshold_relative_density;
    
    double Mag_field_magnitude_scale;

    double Mag_field_power;

    double Mag_field_radial_scale;

};

struct Emission_medium_state_type {

    /*! The density of the emission medium at the current point in units of [g / cm^3]. */
    double Density;

    /*! The temperature of the emission medium at the current point in units of [K]. */
    double Temperature;

    /*! The magnetization of the emission medium at the current point. */
    double Magnetization;

    /*! Pointer to the plasma velocity array at the current point in geometric units. */
    double Plasma_Velocity[4];

    /*! Struct that holds the magnetic field. */
    Magnetic_fields_type Magnetic_fields;

    /*! Enum that specifies which emission medium is being considered. */
    Ensamble_enums Ensamble_type;

};

struct Hotspot_position_type {

    /*! The distance to the hotspot center in geometric units. */
    double Distance;

    /*! The polar angle to the hotspot center in units of [Rad]. */
    double Inclination;

    /*! The azimuthal angle to the hotspot center in units of [Rad].*/
    double Azimuth;

    /*! The X coordinate of the hotspot center in geometric units. */
    double x;

    /*! The Y coordinate of the hotspot center in geometric units. */
    double y;

    /*! The Z coordinate of the hotspot center in geometric units. */
    double z;

};

struct Emission_model_parameters_type {

    /* ================= Thermal Synchrotron Model ================== */
    /* It is fully determined by the electron density and temperature */

    /* ============= Phenomenological Synchrotron Model ============= */

    /*! The emission is scales linearly with this parameter. It has units of [erg / s / sr / Hz]. */
    double Phenomenological_emission_coeff;
    
    /*! The absorbtion scales linearly with this parameter. It has units of [cm^-1]. */
    double Phenomenological_absorbtion_coeff;

    /*! The emission scales as pow(redshift, Phenomenological_emission_power_law). */
    double Phenomenological_emission_power_law;

    /*! The source function scales as pow(redshift, Phenomenological_source_f_power_law). */
    double Phenomenological_source_f_power_law; 

    /* ================== Kappa Synchrotron Model ================== */

    /*! The free dimentionless parameter for the kappa distribution. */
    double Kappa;
};

struct Numerical_metric_potentials_type {

    double F_0;
    double F_1;
    double F_2;
    double W;

};

struct Numerical_metric_params_type {

    double M_ADM;

    // NOTE: This is NOT normalized to the mass.
    double Horizon_radius;

    // NOTE: This NOT normalized to the mass.
    double Horizon_radius_BL;

    // NOTE: This NOT normalized to the mass.
    double a_ADM;

    /* ============ Pointers to the numerical metric potentials control vector arrays ============ = */

    double* F_0_control_vector;
    double* F_1_control_vector;
    double* F_2_control_vector;
    double* W_control_vector;
    int Control_vector_size;

    /* ============ Pointers to the coordinate grid and its control vector arrays ============ = */

    double* Compactified_radial_grid;
    double* Compactified_radial_grid_control_vector;
    int Radial_grid_size;

    double* Theta_grid;
    double* Theta_grid_control_vector;
    int Theta_grid_size;

    Numerical_Anzatz_enums e_Anzatz;

};

struct Metric_parameters_type {

    /*! Enum that specifies the active spacetime. */
    Spacetime_enums e_Spacetime;

    /* ============ Wormhole Specific Parameters ============ = */

    /*! The Wormhole lapse function is given by exp(-M / r - Redshift_param * (M / r)^2). */
    double Redshift_Parameter;

    /*! The Wormhole g_rr function is given by 1. / (1 - R_throat / r). */
    double R_throat;

    /*! Boolean that decides weather to allow photons to cross the Worhmole throat. */
    bool Stop_At_Throat;

    /* ============ Janis-Newman-Winicour Specific Parameters ============ = */

    double JNW_Gamma_Parameter;
    
    /* ============ Gauss-Bonnet Specific Parameters ============ = */

    double GB_Gamma_Parameter;

    /* ============ Regular Black Hole Specific Parameters ============ = */

    double RBH_Parameter;

    /* ============ Black Hole w Dark Matter Halo Specific Parameters ============ = */

    double Compactness;
    double Halo_Mass;

    /* ============ Generic Parameters ============ = */

    double Spin; // Only affects Kerr and the Wormhole

    Numerical_metric_params_type Numerical_metric_params;
};

struct Precomputed_e_pitch_angles_type {

    double* sin_electron_pitch_angles;
    double* cos_electron_pitch_angles;

    // === Used in the thermal synchotron emission functions === //

    /*! 1. / sqrt(sin) */
    double* one_over_sqrt_sin;

    /*! 1. / cbrt(sin) */
    double* one_over_cbrt_sin;

    // === Used in the thermal synchotron Faradey functions === //

    /*! 1. / pow(sin, 0.5175) */
    double* one_over_sin_to_0_p_5175;

    /*! 1. / pow(sin, 0.6) */
    double* one_over_sin_to_0_p_6;

    /*! 1. / pow(sin, 0.7515) */
    double* one_over_sin_to_0_p_7515;

    // === Used in the kappa synchotron emission functions === //

    /*! 1. / pow(sin, 7. / 20) */
    double* one_over_sin_to_7_over_20;
};

struct Thermal_transfer_f_arguments_type {

    /*! The dimensionless parameter that the thermal fit functions depend on. */
    double X;

    /*! sqrt(X). */
    double sqrt_X;

    /*! cbrt(X). */
    double cbrt_X;

    /*! pow(X, 0.5175). */
    double X_to_0_p_5175;

    /*! pow(X, 0.6). */
    double X_to_0_p_6;

    /*! pow(X, 0.7515). */
    double X_to_0_p_7515;

    /*! The dimensionless electron temperature at the current point. */
    double T_electron_dim;

    /*! pow(X, 24. / 25). */
    double T_electron_dim_to_24_25;

    /*! Sin of the emission angle. */
    double sin_pitch_angle;

    /*! Cos of the emission angle. */
    double cos_pitch_angle;

    /*! The current photon frequency in units of [Hz]. */
    double frequency;
};

struct Kappa_transfer_f_arguments_type {

    /*! The dimensionless parameter that the kappa fit functions depend on. */
    double X;

    /*! sqrt(X). */
    double sqrt_X;

    /*! cbrt(X). */
    double cbrt_X;

    /*! pow(X, 7. / 20). */
    double X_to_7_over_20;

    /*! The dimensionless electron temperature at the current point. */
    double T_electron_dim;

    /*! Sin of the emission angle. */
    double sin_emission_angle;

    /*! The kappa parameter. */
    double kappa;

    /*! The current photon frequency in units of [Hz]. */
    double frequency;
};

struct Phenomenological_transfer_f_arguments_type {

    /*! The redshift at the current point. */
    double redshift;

    /*! The cyclotron frequency at the current point in units of [Hz]. */
    double f_cyclo;

    /*! The current photon frequency in units of [Hz]. */
    double frequency;

};

struct Transfer_functions_type {

    /*! Array that holds the current emission functions for all polarizations. */
    double Emission_functions[e_Stokes_param_num];

    /*! Array that holds the current Faradey functions for all polarizations. */
    double Faradey_functions[e_Stokes_param_num];

    /*! Array that holds the current absorbtion functions for all polarizations. */
    double Absorbtion_functions[e_Stokes_param_num];

};

struct Metric_type {

    /*! The metric tensor. */
    double Metric[4][4];

    /*! The metric lapse function. */
    double Lapse_function;

    /*! The metric shift function. */
    double Shift_function;

};

struct Integrator_parameters_type {

    /*! The integral gain for the PID controller. */
    double PID_gain_I;

    /*! The proportional gain for the PID controller. */
    double PID_gain_P;

    /*! The derivative gain for the PID controller. */
    double PID_gain_D;

    /*! The k1 gain for the Gustafsson controller. */
    double Gustafsson_k1;

    /*! The k2 gain for the Gustafsson controller. */
    double Gustafsson_k2;

    /*! The maximum allowed relative step increase. */
    double Max_rel_step_increase;

    /*! The minimum allowed relative step increase. */
    double Min_rel_step_increase;

    /*! The initial stepsize. */
    double Init_stepzie;

    /*! The Dormond-Prince error threshold parameter. Right now, this is used as both an absolute and relative thresholds. */
    double RK_45_accuracy;

    /*! A multiplicative factor forr the integration step in the range (0, 1] that makes the integrator more stable. */
    double Safety_1;

    /*! A additive factor <<1 that ensures no "division by 0" problems occur when calculating the new stepsizes. */
    double Safety_2;

    /*! The error theshold parameter for the Simpson integral solving method. Currently this is only used in the Page-Thorne flux integral. */
    double Simpson_accuracy;

    /*! The maximum allowed attempted integration steps. */
    int Max_integration_count;

    /*! Enum that decides which step controller to use. */
    Step_controller_type_enums Controller_type;

};

struct Observer_parameters_type {

    /*! The moment of observation, relative to the central object dynamics, in geometric units. */
    double init_time;

    /*! The distance to the observer in geometric units. */
    double distance;

    /*! The observer's polar angle in units of [Rad]. */
    double inclination;

    /*! The azimuthal angle of the observer in units of [Rad]. */
    double azimuth;

    /*! The minimum X coordinate of the observer's image. */
    double x_min;

    /*! The maximum X coordinate of the observer's image. */
    double x_max;

    /*! The minimum Y coordinate of the observer's image. */
    double y_min;

    /*! The maximum Y coordinate of the observer's image. */
    double y_max;

    /*! The observation frequency in units of [Hz]. */
    double obs_frequency;

    /*! The observer's rotation angle along the line of sight in units of [Rad].*/
    double cam_rotation_angle;

    /*! The resolution along the X axis of the observer's image. */
    int resolution_x;

    /*! The resolution along the Y axis of the observer's image. */
    int resolution_y;

    /*! Flag that decides weather to include polarization in the radiative transfer. */
    bool include_polarization;

};

struct File_manager_parameters {

    /*! The full path for the input file for simulation mode 2. */
    std::string Sim_mode_2_imput_path;

    /*! The full path for the output file directory. The actual simulation results get stored in the
        ...Output_file_directory\\Simulation_name folder. */
    std::string Output_file_directory;

    /*! The common part of the names for all output files. These get _nX appended to them, where X is the image order. */
    std::string Common_file_names;

    /*! The full path to the vertex shader file. */
    std::string Vert_shader_path;

    /*! The full path to the fragment shader file. */
    std::string Frag_shader_path;

    /*! THe simulation name. */
    std::string Simulation_name;

    /*! Flag that decieds weather to trncate the output files or not. */
    bool Truncate_files;

};

struct Initial_conditions_type {

    /*! Struct that holds the file manager parameters. */
    File_manager_parameters File_manager_params;

    /*! Struct that holds the disk model parameters. */
    Disk_model_parameters_type Disk_params;

    /*! Struct that holds the inital metric (evaluated at the observer). */
    Metric_type Init_metric;

    /*! Struct that holds the hotspot model parameters. */
    Hotspot_model_parameters_type Hotspot_params;

    /*! Struct that holds the integrator parameters. */
    Integrator_parameters_type Integrator_params;

    /*! Struct that holds the observer parameters. */
    Observer_parameters_type Observer_params;

    /*! Struct that holds the metric parameters. */
    Metric_parameters_type Metric_parameters;

    /*! Struct that holds the emission model parameters. */
    Emission_model_parameters_type Emission_params;

    /*! Initial condition for simulation mode 3. The image X coordiante is used to compute the initial azimuthal momentum. */
    double Sim_mode_3_X_init;

    /*! Initial condition for simulation mode 3. The image Y coodirnate is used to compute the initial polar momentum. */
    double Sim_mode_3_Y_init;

    /*! The initial photon momentum. */
    double Init_Momentum[4];

    /*! The central object mass in units of [M_sun]. */
    double central_object_mass;

    /*! Integer that specifies the currently active simulation mode. */
    int Simulation_mode;

    /*! The number of metric parameter values that will be used during simulation mode 2.
        I usually use this mode to do "parameter sweeps", so this comes in handy when parsing the results. */
    int Sim_mode_2_param_value_number;

    /*! The number of samples of the electron pitch angle to use when averaging the emission.
        NOTE: This only has an effect if Average_electron_pitch_angle is set to "ture". */
    int Emission_pitch_angle_samples_to_average;

    /*! Boolean flag that decides weather to print the simulation metadata to the console. */
    bool Print_to_console;

    /*! Boolean flag that decides weather to average the emission over the electron pitch angle. */
    bool Average_electron_pitch_angle;

    /*! The maximum image order for which the radiative transfer will be evaluated. 
        NOTE: For image orders above this, the ray is still propagated, along with the parallel transport 
        of the polarization vector */
    int Max_order;

    /*! Boolean flag that controls weather to simply add the hotspot and disk density and temperature 
        (effectively treating them as one single medium). */
    bool Thermalize_emission_medium;

};

struct Simulation_Context_type {

    /*! Pointer to the struct that holds the initial conditions. */
    Initial_conditions_type* p_Init_Conditions;

    /*! Pointer to the class that holds all the spacetime related functions. */
    Spacetime_Base_Class* p_Spacetime;

    /*! Pointer to the class that holds all the observer related functions. */
    Observer_class* p_Observer;

    /*! Pointer to the class that holds all the emission medium related functions. */
    Emission_models_class* p_Emission_Model;

    /*! Pointer to the class that holds all the Page-Thorne related functions. */
    Page_Thorne_Model_class* p_PT_model;

    /*! Pointer to the class that holds all the file manager related functions. */
    File_manager_class* File_manager;

};

struct s_Ray_log_type {

    /* Pointer to the arrays that hold the intensity and optical depth along the photon trajectory. */
    double* Ray_emission_log[4];

    /* Pointer to the array that holds the entire photon trajectory. */
    double *Ray_path_log;

    /* Int that specifies where in the log to write. 
       This exists for the sole purpose of minimizing the number of arguments in the functions that write to the photon log. */
    int Log_offset;

    /* The length of the phoyon log. */
    int Log_length;

};

struct Results_type {

    /* The struct that tholds the metric parameters. */
    Metric_parameters_type Metric_parameters;

    /* The struct that holds the photon log. */
    s_Ray_log_type Ray_log_struct;

    /* Array that holds the integrated intensity for each image order and polarization. */
    double Intensity[e_order_number][e_Stokes_param_num]{};

    /* Array that holds the source coordinates from the Page-Thorne disk for each image order.
       This exists for use in simulation mode 2. */
    double Source_Coords[4][e_order_number]{};

    /* Array that holds the photon momentum at the source from the Page-Thorne disk for each image order.
       This exists for use in simulation mode 2. */
    double Photon_Momentum[4][e_order_number]{};

    /* Array that holds the Page-Thorne disk flux for each image order. */
    double Flux_PT[e_order_number]{};

    /* Array that holds the Page-Thorne disk redshift for each image order. */
    double Redshift_PT[e_order_number]{};

    /* Placeholder for an array that will hold the integrated optical depth for each image order and polarization. */
    double Optical_Depth{};

    /* Array that holds the coordinates of the image on the observer plane. 
       NOTE: These get affected by the cam_rotation_angle parameter of the observer. */
    double Image_Coords[2]{};

    double Length_scale;

    double Intensity_scale;

};