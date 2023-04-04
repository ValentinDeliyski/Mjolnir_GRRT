#pragma once

#ifndef INPUTS

    #define INPUTS
    #define _USE_MATH_DEFINES

    #include <math.h>
    #include <string>
    #include "Enumerations.h"

    typedef double const Real;

    // ======================== Spacetime Inputs ======================== //

    const Spacetime_enums e_metric = Gauss_Bonnet;	// Spacetime to be used

    Real MASS = 1.0;
    Real SPIN = 0.98;

    // Wormhole spacetime parameters //

    Real WH_REDSHIFT = 2.0;
    Real WH_R_THROAT = MASS;

    const bool STOP_AT_THROAT = false;
 
    // Regular Black Hole spacetime parameters //

    Real RBH_PARAM = 0.00;

    // Janis - Newman - Winicour Naked Singularity spacetime parameters //

    Real JNW_GAMMA		   = 0.48;
    Real JNW_R_SINGULARITY = 2 * MASS / JNW_GAMMA;

    // Gauss - Bonnet spacetime parameters //

    Real GAUSS_BONNET_GAMMA = 1.1495190525;

    // Black Hole with Dark Matter Halo parameters //

    Real COMPACTNESS = 1e-4;
    Real M_HALO      = 1e2;  // Should be either one of [1e2, 1e4]
    Real A_0         = M_HALO / COMPACTNESS;

    // ======================== Observer Inputs ======================== //

    Real r_obs	   = 1e4;			   // Radial potision of the observer [ M ]
    Real theta_obs = 60. / 180 * M_PI; // Polar angle of the observer [ Rad ]
    Real phi_obs   = 0;				   // Azimuthal angle of the observer ( not used ) [ Rad ]

    // ======================== Emission Model Inputs ======================== //

    const Emission_model_enums     e_emission   = Synchotron_phenomenological; // Emission model to be used
    const Disk_density_model_enums e_disk_model = Power_law;

    // Novikov - Thorne accretion disk parameters

    Real r_in  = 4.896029367;	// Inner accretion idsk radius [ M ]
    Real r_out = 35;            // Outer accretion disk radius [ M ]

    // Exponential disk profile parameters //

    Real DISK_HEIGHT_SCALE = 1. / (10. / 3); // disk density ~ exp( - cos(theta)^2 / DISK_HEIGHT_SCALE^2 / 2 )
    Real DISK_RADIAL_SCALE = 10;		      // disk density ~ exp( - r^2 / DISK_RADIAL_SCALE^2 / 2)

    // Power law disk profile parameters //

    Real DISK_OPENING_ANGLE = 0.1;	 // disk density ~ exp( - ctan(theta)^2 / DISK_OPENING_ANGLE^2 / 2)
    Real DISK_CUTOFF_SCALE  = 2;     // disk density ~ exp( - (r - r_ISCO)^2 / DISK_CUTOFF_SCALE^2) if r < r_ISCO
    Real R_0                = 1;

    // Phenomenological Synchotron emission parameters //
    
    Real EMISSION_POWER_LAW    = 0;   // emission   ~ pow( redshity, EMISSION_POWER_LAW )
    Real SOURCE_F_POWER_LAW	   = 2.5; // absorbtion ~ pow( redshity, SOURCE_F_POWER_LAW + EMISSION_POWER_LAW )
    Real DISK_ABSORBTION_COEFF = 0;   // absorbtion ~ DISK_ABSORBTION_COEFF

    Real EMISSION_SCALE_PHENOMENOLOGICAL = 3e-18;

    /* Exact Synchotron emission parameters */
    
    Real DISK_MAGNETIZATION    = 0.01;
    Real MAG_FIELD_GEOMETRY[3] = { 1, 0, 0 };

    Real N_ELECTRON_EXACT_CGS = 2e6;
    Real T_ELECTRON_EXACT_CGS = 1e11;

    /* Hotspot paramteres */

    Real HOTSPOT_REL_SCALE  = 0. / 2; // Hotspot density ~ HOTSPOT_REL_SCALE
    Real HOTSPOT_SCALE      = 1.;
    Real HOTSPOT_R_COORD    = 6.;
    const int HOTSPOT_ANIMATION_NUMBER = 8;

    // ======================== Simulation Modes Inputs ======================== //

    const int Active_Sim_Mode = 1;
    const bool Truncate_files = true;

    // Simulation Mode 1 and 3 viewing window //

    Real V_angle_min = -atan(30 / r_obs);
    Real V_angle_max =  atan(30 / r_obs);

    Real H_angle_min = -atan(30 / r_obs);
    Real H_angle_max =  atan(30 / r_obs);

    const int RESOLUTION = 512;                  // Linear size of the square pixel grid that makes up the image
    const int TEXTURE_BUFFER_SIZE = 3000 * 3000; // The size of the buffer to store the texture - must be >= RESOLUTION^2

    Real Scan_Step = (H_angle_max - H_angle_min) / (RESOLUTION - 1); // The angular step when iterating photons

    // Sim Mode 2 Input File Path //

    const std::string input_file_path = "C:\\Users\\Valentin\\Documents\\University stuff\\General Relativity\\Polarization\\Schwarzschild_Impact_parameters\\First_relativistic\\geodesic_data_20_deg_Sch_r6_50_photons.txt";

    // Sim Mode 4 Initial Conditions //

    Real X_INIT = 1.;
    Real Y_INIT = 1.;

    // ======================== Integrator Inputs ======================== //

    Real INIT_STEPSIZE     = 1e-5;  // > 0 otherwise not really important (unless you put the observer at r_obs > 1e6)
    Real INTEGRAL_ACCURACY = 5e-9;  // Used to compute the Flux integral for the Novikov-Thorne model - this value seems to be good
    Real RK45_ACCURACY     = 1e-8;  // 1e-8 Seems to be an opitimal tradeoff between accuracy and performace 
    Real SAFETY_1          = 0.8;   // Value between 0 and 1, used for scaling the integration step - between 0.8 and 0.9 is optimal
    Real SAFETY_2          = 1e-16; // Near zero positive number used to avoid division by 0 when calculating the integration step

#endif
