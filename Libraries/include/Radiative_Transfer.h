#pragma once

#define _USE_MATH_DEFINES

#include "Enumerations.h"
#include "Constants.h"
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
