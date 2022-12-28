#pragma once

#ifndef CONSTANTS

    #define CONSTANTS

    typedef const double Real;
    typedef const int  Const_int;
    typedef const bool Const_bool;

    Real INIT_STEPSIZE	   = 1e-5; // > 0 otherwise not really important (unless you put the observer at r_obs > 1e6)
    Real INTEGRAL_ACCURACY = 5e-9;  // Used to compute the Flux integral for the Novikov-Thorne model - this value seems to be good
    Real RK45_ACCURACY	   = 1e-7;  // 1e-7 Seems to be an opitimal tradeoff between accuracy and performace 
    Real SAFETY_1		   = 0.8;   // Value between 0 and 1, used for scaling the integration step - between 0.8 and 0.9 is optimal
    Real SAFETY_2		   = 1e-16; // Near zero positive number used to avoid division by 0 when calculating the integration step

    Const_int RK45_size			    = 7;	   // Number of integration sub-steps
    Const_int MAX_INTEGRATION_COUNT = 7600000; // Realistically the program will never reach this many integration steps, but I prefer this to an infinite loop

    Real Coeff_deriv[RK45_size][RK45_size - 1] =
    {
        {       0,               0,              0,             0,            0,           0   },
        {    1. / 5,             0,              0,             0,            0,           0   },
        {    3. / 40,         9. / 40,           0,             0,            0,           0   },
        {   44. / 45,       -56. / 15,       32. / 9,           0,            0,           0   },
        {19372. / 6561,  -25360. / 2187,  64448. / 6561,  -212. / 729,        0,           0   },
        { 9017. / 3168,    -355. / 33,    46732. / 5247,    49. / 176, -5103. / 18656,     0   },
        {   35. / 384,           0,         500. / 1113,   125. / 192, -2187. / 6784,  11. / 84}
    };

    Real Coeff_sol[RK45_size] = { 35. / 384, 0, 500. / 1113, 125. / 192, -2187. / 6784, 11. / 84, 0 };

    Real Coeff_test_sol[RK45_size] = { 5179. / 57600, 0, 7571. / 16695, 393. / 640, -92097. / 339200, 187. / 2100, 1. / 40 };

    /*

    Physical constants in SI

    */

    Real M_SUN_SI = 1.989e30;
    Real OBJECT_MASS_SI = 4e6 * M_SUN_SI;

    Real G_NEWTON_SI = 6.6743e-11;

    Real M_ELECTRON_SI = 9.1093837e-31;
    Real Q_ELECTRON_SI = 1.60217663e-19;

    Real M_PROTON_SI = 1.67262192e-27;

    Real C_LIGHT_SI         = 299792458;
    Real BOLTZMANN_CONST_SI = 1.380649e-23;

    Real PLANCK_CONSTANT_SI = 6.62607015e-34;

    Real OBS_FREQUENCY_SI = 230e9;

    /*

    Physical constants in CGS

    */

    Real METER_TO_CM = 100;
    Real MASS_TO_CM  = G_NEWTON_SI * OBJECT_MASS_SI / C_LIGHT_SI / C_LIGHT_SI * METER_TO_CM;

    Real M_ELECTRON_CGS = 9.1094e-28;
    Real Q_ELECTRON_CGS = 4.8032e-10;

    Real M_PROTON_CGS = 1.67262192e-24;

    Real C_LIGHT_CGS         = 2.99792458e10;
    Real BOLTZMANN_CONST_CGS = 1.380649e-16;

    Real PLANCK_CONSTANT_CGS = 6.626196e-27;

    Real OBS_FREQUENCY_CGS = 230e9;

    Real CGS_TO_JANSKY = 1e+23;

    /*

    Metric Parameters

    */

    Real MASS = 1.0;
    Real SPIN = 1.0;

    Real WH_REDSHIFT = 2.0;
    Real WH_R_THROAT = MASS;

    Real RBH_PARAM = 0.00;

    Real JNW_R_SINGULARITY = 4.5;
    Real JNW_GAMMA = 2 * MASS / JNW_R_SINGULARITY;

    /*
    
    Optically Thin Disk Paramters
    
    */

    /* Phenomenological model parameters */

    Real EMISSION_POWER_LAW = 0;
    Real SOURCE_F_POWER_LAW = 2.5;

    Real DISK_ABSORBTION_COEFF = 1e5;
    Real DISK_HEIGHT_SCALE = 100. / 3;

    Real EMISSION_SCALE_PHENOMENOLOGICAL = 3e-18;

    /* Exact model parameters */

    Real DISK_OPENING_ANGLE = 1;
    Real DISK_CUTOFF_SCALE  = 10;

    Real DISK_MAGNETIZATION = 0.01;
    Real MAG_FIELD_GEOMETRY[3] = {1, 0, 0};

    Real N_ELECTRON_EXACT_CGS = 2e6;
    Real T_ELECTRON_EXACT_CGS = 1e11;

    Const_int RESOLUTION = 3000 * 3000;

#endif