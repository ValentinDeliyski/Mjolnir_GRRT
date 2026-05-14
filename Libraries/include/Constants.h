#pragma once
#include <complex>
#include <numbers>
typedef double const Real;

const int Nyström_size = 6;

Real Nyström_Deriv_coeffs[Nyström_size][Nyström_size] = { {   0.,      0.,         0.,       0.,   0., 0.},
                                                          {1. / 3,     0.,         0.,       0.,   0., 0.},
                                                          {4. / 25, 6. / 25,       0.,       0.,   0., 0.},
                                                          {1. / 4,    -3.,     15. / 4,      0.,   0., 0.},
                                                          {2. / 27, 10. / 9,  -50. / 81,  8. / 81, 0., 0.},
                                                          {2. / 25, 12. / 25,   2. / 15,  8. / 75, 0., 0.} };

Real Nyström_Coeff_sol[Nyström_size] = { 23. / 192 , 0., 125. / 192, 0., -27. / 64, 125. / 192 };
Real Nyström_Coeff_param[Nyström_size] = { 0 , 1. / 3, 2. / 5, 1., 2. / 3, 4. / 5 };

constexpr double Minkowski_Metric[4][4] = { {-1., 0., 0., 0.},
                                            { 0., 1., 0., 0.},
                                            { 0., 0., 1., 0.},
                                            { 0., 0., 0., 1.} };

constexpr std::complex<double> complex_i = { 0.0, 1.0 };

/*

Physical constants in SI

*/

constexpr double M_SUN_SI = 1.989e30;

constexpr double G_NEWTON_SI = 6.6743e-11;

constexpr double M_PROTON_SI = 1.6726219e-27;
constexpr double M_ELECTRON_SI = 9.1093837e-31;
constexpr double Q_ELECTRON_SI = 1.60217663e-19;

constexpr double C_LIGHT_SI         = 299792458;
constexpr double BOLTZMANN_CONST_SI = 1.380649e-23;

constexpr double PLANCK_CONSTANT_SI = 6.62607015e-34;

/*

Physical constants in CGS

*/

constexpr double METER_TO_CM = 100;

constexpr double M_ELECTRON_CGS = 9.1094e-28;
constexpr double Q_ELECTRON_CGS = 4.8032e-10;

constexpr double M_PROTON_CGS = 1.67262192e-24;

constexpr double C_LIGHT_CGS         = 2.99792458e10;
constexpr double BOLTZMANN_CONST_CGS = 1.380649e-16;

constexpr double PLANCK_CONSTANT_CGS = 6.626196e-27;

constexpr double OBS_FREQUENCY_CGS = 230e9;

constexpr double CGS_TO_JANSKY = 1e+23;

constexpr double MASS_TO_CM = M_SUN_SI * G_NEWTON_SI / C_LIGHT_SI / C_LIGHT_SI * METER_TO_CM;

constexpr double MIN_INTENSITY_THRESHOLD = 1e-40;