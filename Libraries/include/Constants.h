#pragma once
#include <complex>

typedef double const Real;

const int Nyström_size = 6;

Real Nyström_Deriv_coeffs[Nyström_size][Nyström_size] = { {   0.,      0.,         0.,       0.,   0., 0.},
                                                          {1. / 3,     0.,         0.,       0.,   0., 0.},
                                                          {4. / 25, 6. / 25,       0.,       0.,   0., 0.},
                                                          {1. / 4,    -3.,     15. / 4,      0.,   0., 0.},
                                                          {2. / 27, 10. / 9,  -50. / 81,  8. / 81, 0., 0.},
                                                          {2. / 25, 12. / 25,   2. / 15,  8. / 75, 0., 0.} };

Real Nyström_Coeff_sol[Nyström_size] = { 23. / 192 , 0., 125. / 192, 0., -27. / 64, 125. / 192 };

const int RK78_size	= 13;   // Number of integration sub-steps

const int ESDIRK54_size = 7;

Real RK78_Coeff_sol[RK78_size]      = { 41. / 840, 0, 0, 0, 0, 34. / 105, 9. / 35, 9. / 35, 9. / 280, 9. / 280, 41. / 840,     0,         0};
Real RK78_Coeff_test_sol[RK78_size] = {     0,     0, 0, 0, 0, 34. / 105, 9. / 35, 9. / 35, 9. / 280, 9. / 280,     0,     41. / 840, 41. / 840 };

const std::complex<double> complex_i = { 0.0, 1.0 };

/*

Physical constants in SI

*/

Real M_SUN_SI = 1.989e30;

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

Real M_ELECTRON_CGS = 9.1094e-28;
Real Q_ELECTRON_CGS = 4.8032e-10;

Real M_PROTON_CGS = 1.67262192e-24;

Real C_LIGHT_CGS         = 2.99792458e10;
Real BOLTZMANN_CONST_CGS = 1.380649e-16;

Real PLANCK_CONSTANT_CGS = 6.626196e-27;

Real OBS_FREQUENCY_CGS = 230e9;

Real CGS_TO_JANSKY = 1e+23;

/* This constant is used in the dimentionless radiative transfer equations as a density scale in units of [g/cm^3]. */
Real Global_density_scale = 1.0e6; 