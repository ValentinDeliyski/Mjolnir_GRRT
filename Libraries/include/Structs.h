#pragma once

#ifndef STRUCTS

    #define STRUCTS
    #include "Inputs.h"
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

    struct Emission_law_parameters {

        /* Phenomenological emission model parameters */

        double Emission_scale;
        double Absorbtion_coeff;    
        double Emission_power_law;  // emission   ~ pow( redshity, EMISSION_POWER_LAW )
        double Source_f_power_law;  // absorbtion ~ pow( redshity, SOURCE_F_POWER_LAW + EMISSION_POWER_LAW )
        
    };

    struct Precomputed_e_pitch_angles {

        double sin_electron_pitch_angles[NUM_SAMPLES_TO_AVG]{};
        double one_over_sqrt_sin[NUM_SAMPLES_TO_AVG];
        double one_over_cbrt_sin[NUM_SAMPLES_TO_AVG];
        double one_over_sin_to_1_6[NUM_SAMPLES_TO_AVG];

    };

    struct Metric_type {

        double Metric[4][4];
        double Lapse_function;
        double Shift_function;

    };

    class Spacetime_Base_Class;
    class Optically_Thin_Toroidal_Model;
    class Novikov_Thorne_Model;

    struct Initial_conditions_type {

        double init_metric[4][4];
        double init_metric_Redshift_func;
        double init_metric_Shitft_func;

        double init_Pos[3];
        double init_Three_Momentum[3];

        Spacetime_Base_Class* Spacetimes[SPACETIME_NUMBER];
        Optically_Thin_Toroidal_Model* OTT_model;
        Novikov_Thorne_Model* NT_model;

    };

    struct Results_type {

        double Flux_NT[ORDER_NUM]{};
        double Redshift_NT[ORDER_NUM]{};

        double Intensity[ORDER_NUM]{};
        double Optical_Depth{};

        double Source_Coords[3][ORDER_NUM]{};
        double Three_Momentum[3][ORDER_NUM]{};

        double Image_Coords[2]{};

        std::vector<std::vector<double>> Ray_log;

        double Parameters[SPACETIME_NUMBER]{};

    };

    struct Trace_Results_type {

        double Intensity;
        double Optical_Depth;

        double Coordinates[3];
        double Three_Momentum[3];

        double Image_Coords[2];

    };

#endif