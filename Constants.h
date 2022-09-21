#pragma once

#ifndef CONSTANTS

    #define CONSTANTS

    #define _USE_MATH_DEFINES

    typedef double const RK45;
    typedef int const RK45_int;


    RK45 INIT_STEPSIZE	   = 1e-5;  // > 0 otherwise not really important
    RK45 INTEGRAL_ACCURACY = 5e-6;  // Used to compute the Flux integral for the Novikov-Thorne model - this value seems to be good
    RK45 RK45_ACCURACY	   = 1e-7;  // 1e-9 Seems to be an opitimal tradeoff between accuracy and performace 
    RK45 SAFETY_1		   = 0.8;   // Value between 0 and 1, used for scaling the integration step - between 0.8 and 0.9 is optimal
    RK45 SAFETY_2		   = 1e-16; // Near zero positive number used to avoid division by 0 when calculating the integration step

    RK45_int RK45_size			   = 7;		  // Number of integration sub-steps
    RK45_int MAX_INTEGRATION_COUNT = 7600000; // Realistically the program will never reach this many integration steps, but I prefer this to an infinite loop

    RK45 Coeff_deriv[RK45_size][RK45_size-1] =
    {
        {       0,               0,              0,             0,            0,           0   },
        {    1. / 5,             0,              0,             0,            0,           0   },
        {    3. / 40,         9. / 40,           0,             0,            0,           0   },
        {   44. / 45,       -56. / 15,       32. / 9,           0,            0,           0   },
        {19372. / 6561,  -25360. / 2187,  64448. / 6561,  -212. / 727,        0,           0   },
        { 9017. / 3168,    -355. / 33,    46732. / 5247,    49. / 176, -5103. / 18656,     0   },
        {   35. / 384,           0,         500. / 1113,   125. / 192, -2187. / 6784,  11. / 84}
    };

    RK45 Coeff_sol[RK45_size] = { 35. / 384, 0, 500. / 1113, 125. / 192, -2187. / 6784, 11. / 84, 0 };

    RK45 Coeff_test_sol[RK45_size] = { 5179. / 57600, 0, 7571. / 16695, 393. / 640, -92097. / 339200, 187. / 2100, 1. / 40 };

    /*

    Physical constants in SI

    */

    double const M_ELECTRON_SI = 9.1093837e-31;
    double const Q_ELECTRON_SI = 1.60217663e-19;
    double const T_ELECTRON_SI = 1e11;
    double const N_ELECTRON_SI = 1e12;

    double const M_PROTON_SI = 1.67262192e-27;

    double const C_LIGHT_SI         = 299792458;
    double const BOLTZMANN_CONST_SI = 1.380649e-23;

    double const OBS_FREQUENCY_SI = 230e9;

    /*

    Physical constants in CGS

    */

    double const M_ELECTRON_CGS = 9.1094e-28;
    double const Q_ELECTRON_CGS = 4.8032e-10;
    double const T_ELECTRON_CGS = 1e11;
    double const N_ELECTRON_CGS = 2e6;

    double const M_PROTON_CGS = 1.67262192e-24;

    double const C_LIGHT_CGS         = 2.99792458e10;
    double const BOLTZMANN_CONST_CGS = 1.380649e-16;

    double const OBS_FREQUENCY_CGS = 230e9;

    double const CGS_JANSKY = 1e+23;

#endif