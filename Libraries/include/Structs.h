#pragma once
#include "Enumerations.h"
#include <string>

struct Disk_model_parameters_type {

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

    /*! The peak density value in [g / cm^3]. */
    double Electron_density_scale;

    /*! The peak temperature value in [K]. */
    double Electron_temperature_scale;

    /*! The hotspot magnetization value [-]. */
    double Magnetization;

    /*! The constant magnetic field geometry in the plasma rest frame.
        The components are specified as [B_r, B_theta, B_phi].
        This vector gets normalized when read from the input XML. */
    double Mag_field_geometry[3];

    /* ========= Power law density profile parameters ========= */

    /*! The vertical density profile scales asa exp(-(cotan(theta) / 2. / opening_angle)^2). */
    double Power_law_disk_opening_angle;

    /*! The radial density profile scales as pow(r / R_0, radial_power_lawa). */
    double Power_law_density_R_0;

    /*! The radial density profile scales as pow(r / R_0, radial_power_lawa). */
    double Power_law_density_radial_power_law;

    /*! Under the cutoff radius, the radial density profile gains an additional factor of 
        exp(-(r - r_cutoff)^2 / cutoff_scale^2). */
    double Power_law_density_R_cutoff;

    /*! Under the cutoff radius, the radial density profile gains an additional factor of 
        exp(-(r - r_cutoff)^2 / cutoff_scale^2). */
    double Power_law_density_cutoff_scale;

    /* ========= Exponential law density profile parameters ========= */

    /*! The vertical density profile scales as exp(-(cos(theta) / height_scale)^2). */
    double Exp_law_density_height_scale;

    /*! The radial density profile scales as exp(-(r / radial_scale)^2). */
    double Exp_law_density_radial_scale;

    /* ========= Power law temperature profile parameters ========= */

    /*! The radial temperature profile scales as pow(r / R_0, radial_power_lawa). */
    double Power_law_temperature_R_0;

    /*! The radial temperature profile scales as pow(r / R_0, radial_power_lawa). */
    double Power_law_temperature_radial_power_law;

    /*! Under the cutoff radius, the radial temperature profile gains an additional factor of 
        exp(-(r - r_cutoff)^2 / cutoff_scale^2). */
    double Power_law_temperature_R_cutoff;

    /* Under the cutoff radius, the radial temperature profile gains an additional factor of 
       exp(-(r - r_cutoff)^2 / cutoff_scale^2). */
    double Power_law_temperature_cutoff_scale;

    /* ========= Exponential law temperature profile parameters ========= */

    /*! The vertical temperature profile scales as exp(-(cos(theta) / height_scale)^2). */
    double Exp_law_temperature_height_scale;

    /*! The radial temperature profile scales as exp(-(r / radial_scale)^2). */
    double Exp_law_temperature_radial_scale;

};

struct Magnetic_fields_type {

    double Magnetic_field_geometry[3];
    double B_field_plasma_frame[4];
    double B_field_coord_frame[4];
    double B_field_plasma_frame_norm;

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

    /*! The hotspot position, specified as [Distance, Polar Angle, Azimuth Angle] */
    double Position[3]; 

    /* The hotspot is modelled as a localized Gaussian overdensity.
     * The density, temperature and overall time evolution profiles are specified with
     * their respective standard deviations.
     */

     /*! Standard deviation of the Gaussian density profile. */
    double Density_spread;     

    /*! Standard deviation of the Gaussian temperature profile. */
    double Temperature_spread; 

    /*! Standard deviation of the Gaussian temporal profile. Setting this to zero ignores the time 
       evolution of the hotspot profile. */
    double Temporal_spread;   

    /*! Radius of the hotspot. Only affects the Spherical profile. */
    double Radius; 

    /*! Coordinate time of maximum hotspot density */
    double Coord_time_at_max; 

    /*! The peak density value in [g / cm^3]. */
    double Electron_density_scale;     

    /*! The peak temperature value in [K]. */
    double Electron_temperature_scale; 
    
    double Magnetization; /* The hotspot magnetization value [-]. */

    /*! The constant magnetic field geometry in the plasma rest frame.
        The components are specified as [B_r, B_theta, B_phi]. 
        This vector gets normalized when read from the input XML. */
    double Mag_field_geometry[3];

};

struct Emission_medium_state_type {

    double Density;
    double Temperature;
    double Magnetization;
    double* Plasma_Velocity;
    Magnetic_fields_type Magnetic_fields;

    Ensamble_enums Ensamble_type;

};

struct Hotspot_position_type {

    double Distance;
    double Inclination;
    double Azimuth;

    double x;
    double y;
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

    double X;
    double sqrt_X;
    double cbrt_X;
    double X_to_0_p_5175;
    double X_to_0_p_6;
    double X_to_0_p_7515;

    double T_electron_dim;
    double T_electron_dim_to_24_25;

    double sin_pitch_angle;
    double cos_pitch_angle;

    double frequency;
};

struct Kappa_transfer_f_arguments_type {

    double X;
    double sqrt_X;
    double cbrt_X;
    double X_to_7_over_20;

    double T_electron_dim;

    double sin_emission_angle;

    double kappa;

    double frequency;
};

struct Phenomenological_transfer_f_arguments_type {

    double redshift;
    double f_cyclo;
    double frequency;

};

struct Transfer_functions_type {

    double Emission_functions[STOKES_PARAM_NUM];
    double Faradey_functions[STOKES_PARAM_NUM];
    double Absorbtion_functions[STOKES_PARAM_NUM];

};

struct Metric_type {

    double Metric[4][4];
    double Lapse_function;
    double Shift_function;

};

struct Integrator_parameters_type {

    Step_controller_type_enums Controller_type;

    double PID_gain_I;
    double PID_gain_P;
    double PID_gain_D;
    double Gustafsson_k1;
    double Gustafsson_k2;
    double Max_rel_step_increase;
    double Min_rel_step_increase;
    double Init_stepzie;
    double RK_45_accuracy;
    double Safety_1;
    double Safety_2;
    double Simpson_accuracy;
    int Max_integration_count;

};

struct Observer_parameters_type {

    double distance;
    double inclination;
    double azimuth;

    double x_min;
    double x_max;
    double y_min;
    double y_max;

    int resolution_x;
    int resolution_y;

    double obs_frequency;
    double cam_rotation_angle;

    bool include_polarization;

};

struct NT_parameters_type {

    double r_in;
    double r_out;
    bool evaluate_NT_disk;

};

struct File_manager_parameters {

    std::string Sim_mode_2_imput_path;
    std::string Output_file_directory;
    std::string Common_file_names;
    std::string Vert_shader_path;
    std::string Frag_shader_path;
    std::string Simulation_name;
    bool Truncate_files;

};

struct Initial_conditions_type {

    int Simulation_mode;
    bool Print_to_console;
    int Sim_mode_2_param_value_number;
    int Emission_pitch_angle_samples_to_average;
    bool Average_electron_pitch_angle;
    double Sim_mode_3_X_init;
    double Sim_mode_3_Y_init;

    double init_metric[4][4];
    double init_metric_Redshift_func;
    double init_metric_Shitft_func;
    double init_Three_Momentum[4];
    double central_object_mass;

    Disk_model_parameters_type Disk_params;
    Hotspot_model_parameters_type Hotspot_params;
    Emission_model_parameters_type Emission_params;
    Metric_parameters_type Metric_params;
    Integrator_parameters_type Integrator_params;
    Observer_parameters_type Observer_params;
    NT_parameters_type NT_params;
    File_manager_parameters File_manager_params;

};

class Spacetime_Base_Class;
class Generic_Optically_Thin_Model;
class Novikov_Thorne_Model;
class Observer_class;
class File_manager_class;

struct Simulation_Context_type {

    Initial_conditions_type* p_Init_Conditions;
    Spacetime_Base_Class* p_Spacetime;
    Observer_class*       p_Observer;
    Generic_Optically_Thin_Model* p_GOT_Model;
    Novikov_Thorne_Model* p_NT_model;
    File_manager_class* File_manager;

};

struct s_Ray_log_type {

    double *Ray_path_log;
    double *Ray_emission_log[4];
    int Log_offset;
    int Log_length;

};

struct Results_type {

    double Flux_NT[ORDER_NUM]{};
    double Redshift_NT[ORDER_NUM]{};

    double Intensity[ORDER_NUM][STOKES_PARAM_NUM]{};
    double Optical_Depth{};

    double Source_Coords[4][ORDER_NUM]{};
    double Photon_Momentum[4][ORDER_NUM]{};

    double Image_Coords[2]{};

    s_Ray_log_type Ray_log_struct;

    Metric_parameters_type Parameters{};

};