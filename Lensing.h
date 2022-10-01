#pragma once

#ifndef LENSING

    #define LENSING

    #include "Enumerations.h"
    #include "Spacetimes.h"
    #include "Disk_Models.h"

    Return_Value_enums Lens(double initial_conditions[], bool lens_from_file, std::ofstream data[], std::ofstream momentum_data[], 
                            c_Observer Observer_class, Disk_Models e_Disk_Model, Novikov_Thorne_Model NT_Model,
                            Optically_Thin_Toroidal_Model OTT_Model, std::vector<c_Spacetime_Base*> VECTOR);

#endif 