#pragma once

#define _USE_MATH_DEFINES

#include <string>

typedef double const RK45;
typedef int const RK45_int;

RK45 INIT_STEPSIZE	   = 1e-5;  // > 0 otherwise not really important
RK45 INTEGRAL_ACCURACY = 5e-6;  // Used to compute the Flux integral for the Novikov-Thorne model - this value seems to be good
RK45 RK45_ACCURACY	   = 1e-7;  // 1e-9 Seems to be an opitimal tradeoff between accuracy and performace 
RK45 SAFETY_1		   = 0.8;   // Value between 0 and 1, used for scaling the integration step - between 0.8 and 0.9 is optimal
RK45 SAFETY_2		   = 1e-16; // Near zero positive number used to avoid division by 0 when calculating the integration step

RK45 TOROIDAL_DISK_BOUNDY_DENSITY = 8e-6;

RK45_int RK45_size			   = 7;		 // Number of integration sub-steps
RK45_int MAX_INTEGRATION_COUNT = 7600000; // Realistically the program will never reach this many integration steps, but I prefer this to an infinite loop

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

	e_r		       = 0,
	e_theta        = 1,
	e_phi	       = 2,
	e_phi_FD       = 3,
	e_p_theta      = 4,
	e_p_r		   = 5,
	e_Coord_Number = 6

}Coord_enums;

typedef enum tag_XYZ_enums {

	x = 0,
	y = 1,
	z = 2

}XYZ_enums;

typedef enum tag_Image_Orders {

	direct	  = 0,
	first     = 1,
	second	  = 2,
	third     = 3,
	ORDER_NUM

}Order_enums;

typedef enum tag_Return_Values {

	OK    = 0,
	ERROR = 255

}Return_Value_enums;

std::string Return_Value_String[] = {

	"OK",
	"ERROR"

};

typedef enum tag_Disk_Model {

	Novikov_Thorne			= 0,
	Optically_Thin_Toroidal = 1 

}Disk_Models;

typedef enum tag_Disk_Intersection {

	Outside_Disk  = 0,
	Inside_Disk   = 1,
	Entering_Disk = 2,
	Exiting_Disk  = 3

}Disk_Intersection;

std::string File_Names[] = {

	"Kerr_Data0",
	"Kerr_Data1",
	"Kerr_Data2",
	"Kerr_Data3",

	"Kerr_Momentum_Data0",
	"Kerr_Momentum_Data1",
	"Kerr_Momentum_Data2",
	"Kerr_Momentum_Data3",

	"Wormhole_Data0",
	"Wormhole_Data1",
	"Wormhole_Data2",
	"Wormhole_Data3",

	"Wormhole_Momentum_Data0",
	"Wormhole_Momentum_Data1",
	"Wormhole_Momentum_Data2",
	"Wormhole_Momentum_Data3",

	"RBH_Data0",
	"RBH_Data1",
	"RBH_Data2",
	"RBH_Data3",

	"RBH_Momentum_Data0",
	"RBH_Momentum_Data1",
	"RBH_Momentum_Data2",
	"RBH_Momentum_Data3"

};