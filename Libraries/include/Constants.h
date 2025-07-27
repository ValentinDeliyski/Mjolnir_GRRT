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

// The coefficients are from https://ntrs.nasa.gov/api/citations/19680027281/downloads/19680027281.pdf, table X
Real RK78_Coeff_deriv[RK78_size][RK78_size - 1] =
{
    {       0,         0,         0,          0,            0,            0,          0,          0,        0,         0,     0,  0  },
    {    2. / 27,      0,         0,          0,            0,            0,          0,          0,        0,         0,     0,  0  },
    {    1. / 36,   1. / 12,      0,          0,            0,            0,          0,          0,        0,         0,     0,  0  },
    {    1. / 24,      0,      1. / 8,        0,            0,            0,          0,          0,        0,         0,     0,  0  },
    {    5. / 12,      0,    -25. / 16,   25. / 16,         0,            0,          0,          0,        0,         0,     0,  0  },
    {    1. / 20,      0,         0,       1. / 4,       1. / 5,          0,          0,          0,        0,         0,     0,  0  },
    {  -25. / 108,     0,         0,     125. / 108,   -65. / 27,    125. / 54,       0,          0,        0,         0,     0,  0  },
    {   31. / 300,     0,         0,          0,        61. / 225,    -2. / 9,    13. / 900,      0,        0,         0,     0,  0  },
    {      2.,         0,         0,     -53. / 6,     704. / 45,   -107. / 9,    67. / 90,       3.,       0,         0,     0,  0  },
    {  -91. / 108,     0,         0,      23. / 108,  -976. / 135,   311. / 54,  -19. / 60,   17. / 6,  -1. / 12,      0,     0,  0  },
    { 2383. / 4100,    0,         0,    -341. / 164,  4496. / 1025, -301. / 82, 2133. / 4100, 45. / 82, 45. / 164, 18. / 41,  0,  0  },
    {    3. / 205,     0,         0,          0,            0,        -6. / 41,   -3. / 205,  -3. / 41,  3. / 41,   6. / 41,  0,  0  },
    {-1777. / 4100,    0,         0,    -341. / 164,  4496. / 1025, -289. / 82, 2193. / 4100, 51. / 82, 33. / 164, 12. / 41,  0,  1. }     
};  

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