#pragma once

#ifndef INPUTS

    #define INPUTS
    #define _USE_MATH_DEFINES

    #include <math.h>
    #include <string>
    #include "Enumerations.h"

    typedef double const Real;

    // ======================== Spacetime Inputs ======================== //

    const Spacetime_enums e_metric = Kerr;	// Spacetime to be used

    Real MASS = 1.0;
    Real SPIN = 1e-5;

    // Wormhole spacetime parameters //

    Real WH_REDSHIFT = 2.0;
    Real WH_R_THROAT = MASS;

    // Regular Black Hole spacetime parameters //

    Real RBH_PARAM = 0.00;

    // Janis - Newman - Winicour Naked Singularity spacetime parameters //

    Real JNW_GAMMA		   = 0.48;
    Real JNW_R_SINGULARITY = 2 * MASS / JNW_GAMMA;


    // ======================== Observer Inputs ======================== //

    Real r_obs	   = 1e4;			   // Radial potision of the observer [ M ]
    Real theta_obs = 20. / 180 * M_PI; // Polar angle of the observer [ Rad ]
    Real phi_obs = 0;				   // Azimuthal angle of the observer ( not used ) [ Rad ]


    // ======================== Emission Model Inputs ======================== //

    const Emission_model_enums e_emission = Synchotron_phenomenological; // Emission model to be used

    // Novikov - Thorne accretion disk parameters

    Real r_in  = 3;	  // Inner accretion idsk radius [ M ]
    Real r_out = 35;  // Outer accretion disk radius [ M ]

    // Phenomenological Synchotron emission parameters //
    
    Real EMISSION_POWER_LAW    = 0;   // emission   ~ pow( redshity, EMISSION_POWER_LAW )
    Real SOURCE_F_POWER_LAW	   = 2.5; // absorbtion ~ pow( redshity, SOURCE_F_POWER_LAW + EMISSION_POWER_LAW )
    Real DISK_ABSORBTION_COEFF = 0;   // absorbtion ~ DISK_ABSORBTION_COEFF

    Real DISK_HEIGHT_SCALE = 1. / (100. / 3); // disk density ~ exp( - cos(theta)^2 / DISK_HEIGHT_SCALE^2 / 2 )
    Real DISK_RADIAL_SCALE = 10;		      // disk density ~ exp( - r^2 / DISK_RADIAL_SCALE^2 / 2)

    Real EMISSION_SCALE_PHENOMENOLOGICAL = 3e-18;

    /* Exact Synchotron emission parameters */
    
    Real DISK_OPENING_ANGLE = 1;			 // disk density ~ exp( - tan(theta)^2 / DISK_OPENING_ANGLE^2 / 2)
    Real DISK_CUTOFF_SCALE  = 1. / sqrt(10); // disk density ~ exp( - (r - r_ISCO)^2 / DISK_CUTOFF_SCALE^2) if r < r_ISCO

    Real DISK_MAGNETIZATION    = 0.01;
    Real MAG_FIELD_GEOMETRY[3] = { 1, 0, 0 };

    Real N_ELECTRON_EXACT_CGS = 2e6;
    Real T_ELECTRON_EXACT_CGS = 1e11;


    // ======================== Simulation Modes Inputs ======================== //

    const int Active_Sim_Mode = 2;

    // Simulation Mode 1 viewing window //

    Real V_angle_min = -atan(30 / r_obs); 
    Real V_angle_max =  atan(30 / r_obs);

    Real H_angle_min = -atan(30 / r_obs);
    Real H_angle_max =  atan(30 / r_obs);

    const int RESOLUTION = 2500;                 // Linear size of the square pixel grid that makes up the image
    const int TEXTURE_BUFFER_SIZE = 3000 * 3000; // The size of the buffer to store the texture - must be > RESOLUTION^2

    Real Scan_Step = (H_angle_max - H_angle_min) / (RESOLUTION - 1); // The angular step when iterating photons

    // Sim Mode 2 input file path //

    const std::string input_file_path = "C:\\Users\\Valentin\\Documents\\University stuff\\General Relativity\\Polarization\\Schwarzschild_Impact_parameters\\First_relativistic\\geodesic_data_20_deg_Sch_r6_50_photons.txt";

#endif
