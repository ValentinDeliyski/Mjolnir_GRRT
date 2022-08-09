#pragma once
#define _USE_MATH_DEFINES

typedef double const RK45;

RK45 INIT_STEPSIZE	   = 1e-5;
RK45 INTEGRAL_ACCURACY = 5e-6;
RK45 RK45_ACCURACY	   = 1e-10;
RK45 SAFETY			   = 1e-16;

RK45 Coeff_deriv[7][6] =
{
    {       0,               0,              0,             0,            0,           0   },
    {    1. / 5,             0,              0,             0,            0,           0   },
    {    3. / 40,         9. / 40,           0,             0,            0,           0   },
    {   44. / 45,       -56. / 15,       32. / 9,           0,            0,           0   },
    {19372. / 6561,  -25360. / 2187,  64448. / 6561,  -212. / 727,        0,           0   },
    { 9017. / 3168,    -355. / 33,    46732. / 5247,    49. / 176, -5103. / 18656,     0   },
    {   35. / 384,           0,         500. / 1113,   125. / 192, -2187. / 6784,  11. / 84}
};

RK45 Coeff_sol[7] = { 35. / 384, 0, 500. / 1113, 125. / 192, -2187. / 6784, 11. / 84, 0 };

RK45 Coeff_test_sol[7] = { 5179. / 57600, 0, 7571. / 16695, 393. / 640, -92097. / 339200, 187. / 2100, 1. / 40 };

typedef enum tag_Spacetimes {

	Kerr		   = 0,
	Reg_Black_Hole = 1,
	Wormhole	   = 2

}Spacetimes;

typedef enum tag_Coord_enums {

	e_r		  = 0,
	e_theta   = 1,
	e_phi	  = 2,
	e_phi_FD  = 3,
	e_p_theta = 4,
	e_p_r	  = 5

}Coord_enums;

typedef enum tag_XYZ_enums {

	x = 0,
	y = 1,
	z = 2

}XYZ_enums;
