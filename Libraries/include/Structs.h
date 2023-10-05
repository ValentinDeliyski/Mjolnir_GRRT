#pragma once

#ifndef STRUCTS

    #define STRUCTS
    #include "Inputs.h"

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

#endif