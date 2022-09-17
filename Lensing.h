#pragma once

#ifndef LENSING

    #define LENSING

    #include "Enumerations.h"
    #include "Spacetimes.h"
    #include "Disk_Models.h"

    Return_Value_enums Lens(double initial_conditions[], double M, double metric_parameter, double a, double r_throat, double r_in, double r_out, bool lens_from_file,
                            std::ofstream data[], std::ofstream momentum_data[], e_Spacetimes e_metric, c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class,
                            Disk_Models e_Disk_Model, Novikov_Thorne_Model NT_Model, Optically_Thin_Toroidal_Model OTT_Model);

#endif 