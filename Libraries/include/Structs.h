#pragma once

#ifndef STRUCTS

    #define STRUCTS
    #include "Constants.h"
    #include <vector>

    struct Disk_model_parameters_type {

        Ensamble_enums Ensamble_type;
        Profile_enums Density_profile_type;
        Profile_enums Temperature_profile_type;

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

    struct Hotspot_model_parameters_type {

        Ensamble_enums Ensamble_type;
        Profile_enums Density_profile_type;
        Profile_enums Temperature_profile_type;

        double Position[3];

        // The hotspot is modelled as a Gaussian with an std = Spread
        double Density_spread;
        double Temperature_spread;

        double Electron_density_scale;
        double Electron_temperature_scale;
        double Magnetization;

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

        double sin_electron_pitch_angles[NUM_SAMPLES_TO_AVG]{};
        double cos_electron_pitch_angles[NUM_SAMPLES_TO_AVG]{};

        // Used in the thermal synchotron emission functions

        double one_over_sqrt_sin[NUM_SAMPLES_TO_AVG];
        double one_over_cbrt_sin[NUM_SAMPLES_TO_AVG];

        // Used in the thermal synchotron Faradey functions

        double one_over_sin_to_1_point_035[NUM_SAMPLES_TO_AVG];
        double one_over_sin_to_1_point_2_over_2[NUM_SAMPLES_TO_AVG];

        // Used in the kappa synchotron emission functions

        double one_over_sin_to_7_over_20[NUM_SAMPLES_TO_AVG];
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

        double Gain_I;
        double Gain_P;
        double Gain_D;
        double Init_stepzie;
        double RK_45_accuracy;
        double Safety_1;
        double Safety_2;
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

        double resolution_x;
        double resolution_y;

        double obs_frequency;
        double cam_rotation_angle;

        bool include_polarization;

    };

    struct NT_parameters_type {

        double r_in;
        double r_out;
        bool evaluate_NT_disk;

    };

    struct File_manager_paths_type {

        std::string Sim_mode_2_imput_path;
        std::string Output_file_path;
        std::string Common_file_names;
        std::string Vert_shader_path;
        std::string Frag_shader_path;

    };

    struct Initial_conditions_type {

        double init_metric[4][4];
        double init_metric_Redshift_func;
        double init_metric_Shitft_func;
        double init_Three_Momentum[3];

        Disk_model_parameters_type Disk_params;
        Hotspot_model_parameters_type Hotspot_params;
        Emission_model_parameters_type Emission_params;
        Metric_parameters_type Metric_params;
        Integrator_parameters_type Integrator_params;
        Observer_parameters_type Observer_params;
        NT_parameters_type NT_params;
        File_manager_paths_type File_manager_paths;

        bool Average_electron_pitch_angle;

    };

    class Spacetime_Base_Class;
    class Generic_Optically_Thin_Model;
    class Novikov_Thorne_Model;
    class Observer_class;
    class File_manager_class;

    struct Simulation_Context_type {

        Initial_conditions_type* p_Init_Conditions;

        Spacetime_enums       e_Spacetime;
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

        double Source_Coords[3][ORDER_NUM]{};
        double Three_Momentum[3][ORDER_NUM]{};

        double Image_Coords[2]{};

        s_Ray_log_type Ray_log_struct;

        Metric_parameters_type Parameters{};

    };

#endif