#pragma once

#ifndef STRUCTS

    #define STRUCTS
    #include "Constants.h"
    #include <vector>

    struct Disk_model_parameters {

        /* --- Hotspot Parameters --- */

        double Hotspot_position[3];
        double Hotspot_spread;
        double Hotspot_scale;

        /* --- Density Parameters --- */

        double Density_scale;

        /* Power Law Profile Parameters */

        double Disk_opening_angle;
        double Disk_cutoff_scale;
        double Disk_r_cutoff;
        double Power_law_radial_scale;

        /* Exponential Law Profile Parameters */

        double Exp_law_height_scale;
        double Exp_law_radial_scale;

        /* Temperature Parameters */

        double Temperature_scale;

        /* Magnetic Field Parameters */

        double Magnetization;

    };

    struct Metric_Parameters_type {

        /* ============ Wormhole Specific Parameters ============ = */

        double Redshift_Parameter;
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

    struct Emission_law_parameters {

        /* Phenomenological emission model parameters */

        double Emission_scale;
        double Absorbtion_coeff;    
        double Emission_power_law;  // emission   ~ pow( redshift, EMISSION_POWER_LAW )
        double Source_f_power_law;  // absorbtion ~ pow( redshift, SOURCE_F_POWER_LAW + EMISSION_POWER_LAW )
        
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
    };

    struct Thermal_emission_f_arguments {

        double X;
        double sqrt_X;
        double cbrt_X;

    };

    struct Thermal_faradey_f_arguments {

        double X;
        double X_to_1_point_035;
        double X_to_1_point_2;

    };

    struct Metric_type {

        double Metric[4][4];
        double Lapse_function;
        double Shift_function;

    };

    struct Initial_conditions_type {

        double init_metric[4][4];
        double init_metric_Redshift_func;
        double init_metric_Shitft_func;

        double init_Pos[3];
        double init_Three_Momentum[3];

        Disk_model_parameters Disk_params;
        Emission_law_parameters Emission_params;
        Metric_Parameters_type Metric_Parameters;

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

        double Ray_path_log[MAX_INTEGRATION_COUNT * e_path_log_number];
        double Ray_emission_log[MAX_INTEGRATION_COUNT * 2][STOKES_PARAM_NUM];
        int Log_offset{};
        int Log_length{};

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

        Metric_Parameters_type Parameters{};

    };

#endif