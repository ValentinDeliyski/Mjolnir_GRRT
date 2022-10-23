#pragma once

#ifndef CONSTANTS

    #define CONSTANTS

    #define _USE_MATH_DEFINES

    typedef const double Const_Float;
    typedef const int  Const_int;
    typedef const bool Const_bool;

    Const_Float INIT_STEPSIZE	   = 1e-5;  // > 0 otherwise not really important (unless you put the observer at r_obs > 1e6)
    Const_Float INTEGRAL_ACCURACY  = 5e-9;  // Used to compute the Flux integral for the Novikov-Thorne model - this value seems to be good
    Const_Float RK45_ACCURACY	   = 1e-8;  // 1e-7 Seems to be an opitimal tradeoff between accuracy and performace 
    Const_Float SAFETY_1		   = 0.8;   // Value between 0 and 1, used for scaling the integration step - between 0.8 and 0.9 is optimal
    Const_Float SAFETY_2		   = 1e-16; // Near zero positive number used to avoid division by 0 when calculating the integration step

    Const_int RK45_size			    = 7;	   // Number of integration sub-steps
    Const_int MAX_INTEGRATION_COUNT = 7600000; // Realistically the program will never reach this many integration steps, but I prefer this to an infinite loop

    Const_Float Coeff_deriv[RK45_size][RK45_size-1] =
    {
        {       0,               0,              0,             0,            0,           0   },
        {    1. / 5,             0,              0,             0,            0,           0   },
        {    3. / 40,         9. / 40,           0,             0,            0,           0   },
        {   44. / 45,       -56. / 15,       32. / 9,           0,            0,           0   },
        {19372. / 6561,  -25360. / 2187,  64448. / 6561,  -212. / 729,        0,           0   },
        { 9017. / 3168,    -355. / 33,    46732. / 5247,    49. / 176, -5103. / 18656,     0   },
        {   35. / 384,           0,         500. / 1113,   125. / 192, -2187. / 6784,  11. / 84}
    };

    Const_Float Coeff_sol[RK45_size] = { 35. / 384, 0, 500. / 1113, 125. / 192, -2187. / 6784, 11. / 84, 0 };

    Const_Float Coeff_test_sol[RK45_size] = { 5179. / 57600, 0, 7571. / 16695, 393. / 640, -92097. / 339200, 187. / 2100, 1. / 40 };

    /*

    Physical constants in SI

    */

    Const_Float M_ELECTRON_SI = 9.1093837e-31;
    Const_Float Q_ELECTRON_SI = 1.60217663e-19;
    Const_Float T_ELECTRON_SI = 1e11;
    Const_Float N_ELECTRON_SI = 1e12;

    Const_Float M_PROTON_SI = 1.67262192e-27;

    Const_Float C_LIGHT_SI         = 299792458;
    Const_Float BOLTZMANN_CONST_SI = 1.380649e-23;

    Const_Float OBS_FREQUENCY_SI = 230e9;

    /*

    Physical constants in CGS

    */

    Const_Float M_ELECTRON_CGS = 9.1094e-28;
    Const_Float Q_ELECTRON_CGS = 4.8032e-10;
    Const_Float T_ELECTRON_CGS = 1e11;
    Const_Float N_ELECTRON_CGS = 2e6;

    Const_Float M_PROTON_CGS = 1.67262192e-24;

    Const_Float C_LIGHT_CGS         = 2.99792458e10;
    Const_Float BOLTZMANN_CONST_CGS = 1.380649e-16;

    Const_Float OBS_FREQUENCY_CGS = 230e9;

    Const_Float CGS_TO_JANSKY = 1e+23;

    /*

    Metric Parameters

    */

    Const_Float MASS = 1.0;
    Const_Float SPIN = 0.0001;

    Const_Float WH_REDSHIFT = 1.0;
    Const_Float WH_R_THROAT = MASS;

    Const_Float RBH_PARAM = 0.00;

    Const_Float JNW_R_SINGULARITY = 3.0;
    Const_Float JNW_GAMMA = 2 * MASS / JNW_R_SINGULARITY;

    /*
    
    Optically Thin Disk Paramters
    
    */

    Const_Float DISK_ALPHA = 1;

    Const_Float DISK_HEIGHT_SCALE = 0.1;

    Const_Float DISK_RAD_CUTOFF = 3.0 * MASS;
    Const_Float DISK_OMEGA = MASS / 3.5;

    Const_Float DISK_MAGNETIZATION = 0.01;
    Const_Float MAG_FIELD_GEOMETRY[3] = {1, 0, 0};

#endif