/****************************************************************************************************
|                                                                                                   |
|                          ---------  Gravitational Ray Tracer  ---------                           | 
|                                                                                                   |
|    * @Version: 3.2                                                                                |
|    * @Author: Valentin Deliyski                                                                   |
|    * @Description: This program numeriaclly integrates the equations of motion                    |
|    for null geodesics and ratiative transfer in a curved spacetime and projects                   |
|    them onto an observer's screen to construct relativistic images of accretion disks             |
|                                                                                                   |
|    * @Supported Spacetimes:                                                                       | 
|        ** Kerr Black Holes                                                                        |
|        ** Static Regular Black Holes                                                              |
|        ** Rotating Traversable Wormholes                                                          |
|        ** Janis - Newman - Winicour Naked Singularities                                           |
|                                                                                                   |
|    * @Supported Disk Models                                                                       |
|        ** Novikov-Thorne                                                                          |
|        ** Generic Optically Thin Disk With Arbitrary Density, Emission and Absorbtion Profiles    |
|                                                                                                   |
****************************************************************************************************/

#define _USE_MATH_DEFINES

#include <string>
#include <iostream>
#include <iomanip> 
#include <fstream>
#include <cmath>
#include <string>
#include <vector>

#include "Constants.h"
#include "Spacetimes.h"
#include "Enumerations.h"
#include "IO_files.h"

#include "Disk_Models.h"
#include "General_functions.h"

#include "Lensing.h"

e_Spacetimes e_metric = Naked_Singularity;

/*

Define classes that hold the spacetime properites

*/

std::vector<c_Spacetime_Base*> Spacetimes = {

    new derived_Kerr_class(),
    new derived_RBH_class(),
    new derived_Wormhole_class(),
    new derived_JNW_class()
};

/*

Define the Observer class

*/

extern Const_Float r_obs = 1e4;
extern Const_Float theta_obs = 20. / 180 * M_PI;
Const_Float phi_obs = 0;

c_Observer Observer_class(r_obs, theta_obs, phi_obs);

/*

Define the Optically Thin Disk Class

*/

Optically_Thin_Toroidal_Model OTT_Model(DISK_ALPHA, DISK_HEIGHT_SCALE, DISK_RAD_CUTOFF, DISK_OMEGA, DISK_MAGNETIZATION, MAG_FIELD_GEOMETRY);

/*

Define the Novikov-Thorne Disk Class

*/

Const_Float r_in = 0.9*Spacetimes[e_metric]->get_ISCO(Prograde);
Const_Float r_out = 20 * r_in;

Novikov_Thorne_Model NT_Model(r_in, r_out);

/*

Define some global boolians

*/

extern Const_bool lens_from_file = true;
extern Const_bool truncate       = true;

int main() {

    /*
    
    Create/Open the logging files

    */
        
    std::ofstream data[4], momentum_data[4];

        open_output_files(data, momentum_data);
  
    /*
    
    Get the metric at the observer to feed into the initial conditions functions

    */

    Initial_conditions_type s_Initial_Conditions{};
    s_Initial_Conditions.init_Pos[e_r] = r_obs;
    s_Initial_Conditions.init_Pos[e_theta] = theta_obs;
    s_Initial_Conditions.init_Pos[e_phi] = phi_obs;

    double metric[4][4]{}, N_obs, omega_obs;

        Spacetimes[e_metric]->get_metric(s_Initial_Conditions.init_metric, &N_obs, &omega_obs, r_obs, theta_obs);
    
        s_Initial_Conditions.init_metric_Redshift_func = N_obs;
        s_Initial_Conditions.init_metric_Shitft_func = omega_obs;

    /*
    
    Im not even sure if these status checks do anything...
    
    */

    Return_Value_enums Integration_status = OK;

    print_ASCII_art();

    std::cout << "Observer Radial Position [GM/c^2] = " << r_obs << '\n';
    std::cout << "Observer Inclination [deg]        = " << int(theta_obs / M_PI * 180) << '\n';

    if (lens_from_file) {

        /*
        
        Read the initial conditions from file

        */

        double J_data[500]{}, p_theta_data[500]{};
        int Data_number = 0;

        get_geodesic_data(J_data, p_theta_data, &Data_number);

        if (Integration_status == OK) {

            for (int photon = 0; photon <= Data_number; photon += 1) {

                /*

                This function polulates the initial momentum inside the s_Initial_Conditions struct

                */

                Spacetimes[e_metric]->get_initial_conditions_from_file(&s_Initial_Conditions, J_data, p_theta_data, photon);

                Integration_status = Lens(&s_Initial_Conditions, data, momentum_data);
    
                print_progress(photon, Data_number, lens_from_file);
            }

            std::cout << '\n';
        }
    }
    else{

        /*
        
        Setup a viewing window for the observer and loop trough it

        */

        double V_angle_min = -10 / r_obs;
        double V_angle_max = 10 / r_obs;

        double H_angle_min = -15 / r_obs;
        double H_angle_max = 15 / r_obs;

        double Scan_Step = 2 * H_angle_max / 1000;
   
        int progress = 0;
        
        int V_num = floor(log10f((V_angle_max - V_angle_min) * 10000) + 1);

        if (Integration_status == OK) {

            for (double V_angle = 0; V_angle >= V_angle_min; V_angle -= Scan_Step) {

                print_progress(progress, int((V_angle_max - V_angle_min) / Scan_Step), lens_from_file);

                progress += 1;

                for (double H_angle = H_angle_min; H_angle <= H_angle_max; H_angle += Scan_Step) {

                    /*
                    
                    This function polulates the initial momentum inside the s_Initial_Conditions struct
                    
                    */

                    get_intitial_conditions_from_angles(&s_Initial_Conditions, V_angle, H_angle);

                    Integration_status = Lens(&s_Initial_Conditions, data, momentum_data);

                }

            }
        }
    }            

    close_output_files(data, momentum_data);

    std::string Return_Value_String[] = {

        "OK",
        "ERROR"

    };

    std::cout << "Program Status = " << Return_Value_String[Integration_status];

    return Integration_status;
}