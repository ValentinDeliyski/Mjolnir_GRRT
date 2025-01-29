#pragma once

#define _USE_MATH_DEFINES

#include "Enumerations.h"
#include "Constants.h"
#include "Structs.h"
#include "General_math_functions.h"

void Implicit_Trapezoid_Radiative_Transfer(double const emission_functions[STOKES_PARAM_NUM],
                                           double const absorbtion_functions[STOKES_PARAM_NUM],
                                           double const faradey_functions[STOKES_PARAM_NUM],
                                           double const step,
                                           double Intensity[STOKES_PARAM_NUM]);

void Analytic_Radiative_Transfer(double const emission_functions[STOKES_PARAM_NUM],
                                 double const absorbtion_functions[STOKES_PARAM_NUM],
                                 double const faradey_functions[STOKES_PARAM_NUM],
                                 double const step,
                                 double Intensity[STOKES_PARAM_NUM]);

void RK4_Radiative_Transfer(Transfer_functions_type* const Transfer_functions,
                            double const step,
                            double Stokes_Vector[STOKES_PARAM_NUM]);