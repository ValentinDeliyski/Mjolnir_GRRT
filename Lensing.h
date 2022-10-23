#pragma once

#ifndef LENSING

    #define LENSING

    #include "Enumerations.h"
    #include "Spacetimes.h"
    #include "Disk_Models.h"

	Return_Value_enums RK45(double State_Vector[], double Derivatives[], double* step, double J, bool* continue_integration);

	Return_Value_enums Lens(Initial_conditions_type* p_Initial_Conditions, std::ofstream data[], std::ofstream momentum_data[]);

#endif 