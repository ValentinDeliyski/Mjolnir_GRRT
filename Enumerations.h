#pragma once

#ifndef ENUMS

	#define ENUMS

	typedef enum tag_Spacetime_enums {

		Kerr = 0,
		Reg_Black_Hole = 1,
		Wormhole = 2,
		Naked_Singularity = 3,
		SPACETIME_NUMBER = 4

	}e_Spacetimes;

	typedef enum tag_State_enums {

		e_r = 0,
		e_theta = 1,
		e_phi = 2,
		e_phi_FD = 3,
		e_p_theta = 4,
		e_p_r = 5,
		e_Intensity = 6,
		e_Optical_Depth = 7,
		e_State_Number = 8

	}State_enums;

	typedef enum tag_Spacetime_coords {

		e_t_coord     = 0,
		e_r_coord     = 1,
		e_theta_coord = 2,
		e_phi_coord   = 3

	}Spacetime_coords;

	typedef enum tag_XYZ_enums {

		x = 0,
		y = 1,
		z = 2

	}XYZ_enums;

	typedef enum tag_Image_Orders {

		direct = 0,
		first = 1,
		second = 2,
		third = 3,
		ORDER_NUM = 4,
		Total = 5,

	}Order_enums;

	typedef enum tag_Return_Values {

		OK = 0,
		ERROR = 255

	}Return_Value_enums;

	typedef enum tag_Disk_Model {

		Novikov_Thorne = 0,
		Optically_Thin_Toroidal = 1,
		DISK_MODEL_NUM = 2

	}Disk_Models;

	typedef enum tag_Disk_Intersection {

		Outside_Disk = 0,
		Inside_Disk = 1,
		Entering_Disk = 2,
		Exiting_Disk = 3

	}Disk_Intersection;

	typedef enum tag_Orbit_Orientation {

		Prograde   = 0,
		Retrograde = 1,

	}Orbit_Orientation;

	typedef enum tag_Spacetime_Parameters {

		/* Wormhole */

		e_WH_Redshift_parameter = 0,
		e_WH_Throat_radius	  = 1,

		/* Regular Black Hole */

		e_RBH_parameter = 2,

		/* JNW Naked Singularity */

		e_JNW_r_Singularity = 3,
		e_JNW_Gamma		  = 4,

		PARAMETER_NUM = 5

	}Spacetime_Parameters;


#endif
