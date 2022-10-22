#pragma once

#ifndef LENSING

    #define LENSING

    #include "Enumerations.h"
    #include "Spacetimes.h"
    #include "Disk_Models.h"

    Return_Value_enums Lens(double initial_conditions[], bool lens_from_file, std::ofstream data[], std::ofstream momentum_data[], 
                            c_Observer Observer_class, Novikov_Thorne_Model NT_Model,
                            Optically_Thin_Toroidal_Model OTT_Model, std::vector<c_Spacetime_Base*> VECTOR);

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