#pragma once

#ifndef LENSING

    #define LENSING

    #include "Enumerations.h"
    #include "Spacetimes.h"
    #include "Disk_Models.h"

    void get_Radiative_Transfer(double State_Vector[], double Derivatives[], int iteration, double J);

    void RK45(double State_Vector[], double Derivatives[], double* step, double J, bool* continue_integration);

	void Lens(Initial_conditions_type* p_Initial_Conditions, std::ofstream data[], std::ofstream momentum_data[]);

#endif 