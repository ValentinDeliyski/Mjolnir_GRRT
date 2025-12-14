#pragma once

#define _USE_MATH_DEFINES

#include "Enumerations.h"
#include "Emission_Models.h"
#include "Constants.h"
#include "Spacetimes.h"
#include "Structs.h"
#include "General_math_functions.h"

void Implicit_Trapezoid_Radiative_Transfer(double* const emission_functions,
                                           double* const absorbtion_functions,
                                           double* const faradey_functions,
                                           double const step,
                                           double* const Intensity);

void Analytic_Radiative_Transfer(double* const emission_functions,
                                 double* const absorbtion_functions,
                                 double* const faradey_functions,
                                 double const step,
                                 double* const Intensity);

void RK5_radiative_transfer(double* const State_Vector_Global,
                            const Simulation_Context_type* p_Sim_Context,
                            double* const Stokes_Vector);