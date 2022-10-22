#pragma once

#ifndef LENSING

    #define LENSING

    #include "Enumerations.h"
    #include "Spacetimes.h"
    #include "Disk_Models.h"


	Return_Value_enums RK45(double State_Vector[], double Derivatives[], double* step, double J, bool* continue_integration,
							c_Observer Observer_class, Optically_Thin_Toroidal_Model OTT_Model, std::vector<c_Spacetime_Base*> Spacetimes);

	Return_Value_enums Lens(double initial_conditions[], std::ofstream data[], std::ofstream momentum_data[]);

	typedef struct tag_results {

		double Flux_NT[ORDER_NUM];
		double Redshift_NT[ORDER_NUM];

		double Intensity[ORDER_NUM];
		double Optical_Depth;

		double Source_Coords[3][ORDER_NUM];
		double Three_Momentum[3][ORDER_NUM];

		double Image_Coords[2];

		double Parameters[SPACETIME_NUMBER];

	}results;

#endif 