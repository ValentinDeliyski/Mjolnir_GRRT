#pragma once

#ifndef INPUTS

    #define INPUTS
    #define _USE_MATH_DEFINES

    #include <math.h>
    #include <string>
    #include "Enumerations.h"

    typedef double const Real;

    // ======================== Simulation Modes Inputs ======================== //

    const int Active_Sim_Mode = 1;
    const bool Truncate_files = true;

    // Sim Mode 2 Configuration //

    const int PARAM_SWEEP_NUMBER = 1;
    Real INIT_PARAM_VALUE        = 0.001;
    Real FINAL_PARAM_VALUE       = 1;

    const Metric_Parameter_Selector PARAM_TYPE = WH_Redshift;

    // Sim Mode 4 Initial Conditions //

    Real X_INIT = 1.0f;
    Real Y_INIT = 1.0f;

#endif
