#pragma once

#ifndef CONSTANTS

    #define CONSTANTS
    
    #include "Inputs.h"
    #include <string>

    const int RK45_size			    = 7;   // Number of integration sub-steps
    const int MAX_INTEGRATION_COUNT = 1e9; // Realistically the program will never reach this many integration steps, but I prefer this to an infinite loop

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
    Real OBJECT_MASS_SI = 6.2e9 * M_SUN_SI;

    Real G_NEWTON_SI = 6.6743e-11;

    Real M_ELECTRON_SI = 9.1093837e-31;
    Real Q_ELECTRON_SI = 1.60217663e-19;

    Real M_PROTON_SI = 1.67262192e-27;

    Real C_LIGHT_SI         = 299792458;
    Real BOLTZMANN_CONST_SI = 1.380649e-23;

    Real PLANCK_CONSTANT_SI = 6.62607015e-34;

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
    
    Shader File Paths
    
    */

    const auto vert_shader_path = "C:\\Users\\Valentin\\Documents\\Repos\\Gravitational_Lenser\\Libraries\\shaders\\default.vert";
    const auto frag_shader_path = "C:\\Users\\Valentin\\Documents\\Repos\\Gravitational_Lenser\\Libraries\\shaders\\default.frag";

#endif