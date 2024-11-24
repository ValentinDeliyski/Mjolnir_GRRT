#pragma once
#include "Enumerations.h"
#include <string>

struct Disk_model_parameters_type {

    Ensamble_enums Ensamble_type;
    Profile_enums Density_profile_type;
    Profile_enums Temperature_profile_type;
    Velocity_enums Velocity_profile_type;

    double Electron_density_scale;
    double Electron_temperature_scale;
    double Magnetization;

    double Mag_field_geometry[3];

    /* ----------- Power law density profile parameters ----------- */

    double Power_law_disk_opening_angle;
    double Power_law_density_R_0;
    double Power_law_density_R_cutoff;
    double Power_law_density_cutoff_scale;
    double Power_law_density_radial_power_law;

    /* -------- Exponential law density profile parameters ------- */

    double Exp_law_density_height_scale;
    double Exp_law_density_radial_scale;

    /* --------- Power law temperature profile parameters -------- */

    double Power_law_temperature_R_0;
    double Power_law_temperature_R_cutoff;
    double Power_law_temperature_cutoff_scale;
    double Power_law_temperature_radial_power_law;

    /* -------- Exponential law temperature profile parameters ------- */

    double Exp_law_temperature_height_scale;
    double Exp_law_temperature_radial_scale;

};

struct Magnetic_fields_type {

    double B_field_plasma_frame[4];
    double B_field_coord_frame[4];
    double B_field_plasma_frame_norm;

};

struct Hotspot_model_parameters_type {

    /* Specifies the statistical ensamble of the hotspot. */
    Ensamble_enums Ensamble_type; 

    /* Specifies the density profile of the hotspot. The current supported profiles are:
        - Gaussian
        - Sphere with a constant Radius  */
    Profile_enums Density_profile_type; 

    /* Specifies the temperature profile of the hotspot. The current supported profiles are:
        - Gaussian
        - Sphere with a constant Radius */
    Profile_enums Temperature_profile_type; 

    /* Specifies the velocity profile of the hotspot. */
    Velocity_enums Velocity_profile_type;  

    /* The hotspot potision, specified as [Distance, Polar Angle, Azimuth Angle] */
    double Position[3]; 

    /* The hotspot is modelled as a localized Gaussian overdensity.
     * The density, temperature and overall time evolution profiles are specified with
     * their respective standard deviations.
     */

     /* Standard deviation of the Gaussian density profile. */
    double Density_spread;     

    /* Standard deviation of the Gaussian temperature profile. */
    double Temperature_spread; 

    /* Standard deviation of the Gaussian temporal profile. Setting this to zero ignores the time 
       evolution of the hotspot profile. */
    double Temporal_spread;   

    /* Radius of the hotspot. Only affects the Spherical profile. */
    double Radius; 

    /* Coordinate time of maximum hotspot emission */
    double Coord_time_at_max; 

    /* The peak density value. */
    double Electron_density_scale;     

    /* The peak temperature value. */
    double Electron_temperature_scale; 

    /* The hotspot magnetization value. */
    double Magnetization;   

    /* The constant magnetic field geometry in the plasma rest frame.
       The components are specified as [B_r, B_theta, B_phi]. */
    double Mag_field_geometry[3];

};

struct Emission_model_parameters_type {

    // --------------- Thermal Synchotron Model --------------- //
    // It is fully determined by the electron density and temperature

    // ----------- Phenomenological Synchotron Model ---------- //

    double Phenomenological_emission_coeff;
    double Phenomenological_absorbtion_coeff;
    double Phenomenological_emission_power_law;  // emission   ~ pow( redshift, EMISSION_POWER_LAW )
    double Phenomenological_source_f_power_law;  // absorbtion ~ pow( redshift, SOURCE_F_POWER_LAW + EMISSION_POWER_LAW )

    // ---------------- Kappa Synchotron Model ---------------- //

    double Kappa;
};

struct Metric_parameters_type {

    Spacetime_enums e_Spacetime;

    /* ============ Wormhole Specific Parameters ============ = */

    double Redshift_Parameter;
    double R_throat;
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

struct Precomputed_e_pitch_angles {

    double* sin_electron_pitch_angles;
    double* cos_electron_pitch_angles;

    // Used in the thermal synchotron emission functions

    double* one_over_sqrt_sin;
    double* one_over_cbrt_sin;

    // Used in the thermal synchotron Faradey functions

    double* one_over_sin_to_1_point_035;
    double* one_over_sin_to_1_point_2_over_2;

    // Used in the kappa synchotron emission functions

    double* one_over_sin_to_7_over_20;
};

struct Thermal_emission_f_arguments {

    double X;
    double sqrt_X;
    double cbrt_X;
    double frequency;

};

struct Thermal_faradey_f_arguments {

    double X;
    double X_to_1_point_035;
    double X_to_1_point_2;
    double frequency;

};

struct Kappa_transfer_f_arguments {

    double X;
    double sqrt_X;
    double cbrt_X;
    double X_to_7_over_20;
    double kappa;
    double sin_emission_angle;
    double T_electron_dim;

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