#pragma once

#ifndef INPUTS

    #define INPUTS
    #define _USE_MATH_DEFINES

    #include <math.h>
    #include <string>
    #include "Enumerations.h"

    typedef double const Real;

    // ======================== Spacetime Inputs ======================== //

    const Spacetime_enums e_metric = Naked_Singularity; // Spacetime to be used

    Real MASS = 1.0f;
    Real SPIN = 0.5f;

    // Wormhole spacetime parameters //

    Real WH_REDSHIFT = 2.0f;
    Real WH_R_THROAT = MASS;

    const bool STOP_AT_THROAT = false;
 
    // Regular Black Hole spacetime parameters //

    Real RBH_PARAM = 0.5f;

    // Janis - Newman - Winicour Naked Singularity spacetime parameters //

    Real JNW_GAMMA		   = 0.48;
    Real JNW_R_SINGULARITY = 2 * MASS / JNW_GAMMA;

    // Gauss - Bonnet spacetime parameters //

    Real GAUSS_BONNET_GAMMA = 1.15;

    // Black Hole with Dark Matter Halo parameters //

    Real COMPACTNESS = 1e-4;
    Real M_HALO      = 1e4;  // Should be either one of [1e2, 1e4]
    Real A_0         = M_HALO / COMPACTNESS;

    // ======================== Observer Inputs ======================== //

    Real r_obs	   = 1e4;			     // Radial potision of the observer [ M ]
    Real theta_obs = 160.0 / 180 * M_PI; // Polar angle of the observer [ Rad ]
    Real phi_obs   = 0.0f;			     // Azimuthal angle of the observer ( not used ) [ Rad ]
    Real obs_cam_rotation_angle = -70.0f / 180 * M_PI - M_PI_4;

    // ======================== Emission Model Inputs ======================== //

    const Emission_model_enums     e_emission   = Synchotron_exact; // Emission model to be used
    const Disk_density_model_enums e_disk_model = Power_law;

    // Novikov - Thorne accretion disk parameters

    const bool Evaluate_NT_disk = true;

    Real r_in  = 4.5f;	// Inner accretion idsk radius [ M ]
    Real r_out = 34.0f;    // Outer accretion disk radius [ M ]

    // Exponential disk profile parameters //

    Real DISK_HEIGHT_SCALE = 1.0f / (10.0f / 3.0f); // disk density ~ exp( - cos(theta)^2 / DISK_HEIGHT_SCALE^2 / 2 )
    Real DISK_RADIAL_SCALE = 10.0f;               // disk density ~ exp( - r^2 / DISK_RADIAL_SCALE^2 / 2)

    // Power law disk profile parameters //

    Real DISK_OPENING_ANGLE = 1.0f / 10;  // disk density ~ exp( - ctan(theta)^2 / DISK_OPENING_ANGLE^2 / 2)
    Real DISK_CUTOFF_SCALE  = 0.4f;       // disk density ~ exp( - (r - R_Cutoff)^2 / DISK_CUTOFF_SCALE^2) if r < R_Cutoff
    Real R_Cutoff           = 5.0f;       // = r_isco[Inner] if = NULL or = r_isco[Outer] if < 0
    Real R_0                = 5.0f;

    // Phenomenological Synchotron emission parameters //
    
    Real EMISSION_POWER_LAW    = 0.0f;  // emission   ~ pow( redshity, EMISSION_POWER_LAW )
    Real SOURCE_F_POWER_LAW	   = 2.5;   // absorbtion ~ pow( redshity, SOURCE_F_POWER_LAW + EMISSION_POWER_LAW )
    Real DISK_ABSORBTION_COEFF = 1e-4;  // absorbtion ~ DISK_ABSORBTION_COEFF

    Real EMISSION_SCALE_PHENOMENOLOGICAL = 8e-19 * 0.6 / 1.964;

    // Exact Synchotron emission parameters //
    
    Real DISK_MAGNETIZATION    = 0.01;
    Real MAG_FIELD_GEOMETRY[3] = { 1, 0, 0 }; // Unused for now
    
    Real N_ELECTRON_EXACT_CGS = 5e+05;
    Real T_ELECTRON_EXACT_CGS = 7.2e+10;

    const int NUM_SAMPLES_TO_AVG = 50; // Number of samples used to average the emission function over the electron pitch angles

    // Hotspot paramteres //

    Real HOTSPOT_REL_SCALE  = 0.0f; // Hotspot density ~ HOTSPOT_REL_SCALE
    Real HOTSPOT_SCALE      = 1.0f; // Hotspot density ~ exp(-|r - r_c|^2 / HOTSPOT_SCALE&^2)
    Real HOTSPOT_R_COORD    = 6.0f;
    const int HOTSPOT_ANIMATION_NUMBER = 4;

    // ======================== Simulation Modes Inputs ======================== //

    const int Active_Sim_Mode = 1;
    const bool Truncate_files = true;

    // Simulation Mode 1 and 3 viewing window //

    Real V_angle_min = -atan(15 / r_obs);
    Real V_angle_max =  atan(15 / r_obs);

    Real H_angle_min = -atan(15 / r_obs);
    Real H_angle_max =  atan(15 / r_obs);

    const int RESOLUTION = 1024;                             // Linear size of the square pixel grid that makes up the image
    const int NUM_RAYS = RESOLUTION * RESOLUTION; // The size of the buffer to store the texture

    Real Scan_Step = (H_angle_max - H_angle_min) / (RESOLUTION - 1); // The angular step when iterating photons

    // Sim Mode 2 Input File Path //

    const std::string input_file_path = "C:\\Users\\Valentin\\Documents\\University stuff\\General Relativity\\Polarization\\Schwarzschild_Impact_parameters\\First_relativistic\\geodesic_data_20_deg_Sch_r6_50_photons.txt";

    // Sim Mode 4 Initial Conditions //

    Real X_INIT = 1.0f;
    Real Y_INIT = 1.0f;

    // ======================== Integrator Inputs ======================== //

    Real INIT_STEPSIZE     = 1e-5;  // > 0 otherwise not really important (unless you put the observer at r_obs > 1e6)
    Real INTEGRAL_ACCURACY = 1e-6;  // Used to compute the Flux integral for the Novikov-Thorne model - this value seems to be good
    Real RK45_ACCURACY     = 1e-9;  // 1e-9 Seems to be an opitimal tradeoff between accuracy and performace 
    Real SAFETY_1          = 0.8;   // Value between 0 and 1, used for scaling the integration step - between 0.8 and 0.9 is optimal
    Real SAFETY_2          = 1e-16; // Near zero positive number used to avoid division by 0 when calculating the integration step

#endif
