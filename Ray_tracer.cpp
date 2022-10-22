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

e_Spacetimes e_metric = Kerr;

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
extern Const_Float theta_obs = 85. / 180 * M_PI;
Const_Float phi_obs = 0;

c_Observer Observer_class(r_obs, theta_obs, phi_obs);

/*

Define the Optically Thin Disk Class

*/

Optically_Thin_Toroidal_Model OTT_Model(DISK_ALPHA, DISK_HEIGHT_SCALE, DISK_RAD_CUTOFF, DISK_OMEGA, DISK_MAGNETIZATION, MAG_FIELD_GEOMETRY);

/*

Define the Novikov-Thorne Disk Class

*/

Const_Float r_in = Spacetimes[e_metric]->get_ISCO(Prograde);
Const_Float r_out = 20 * r_in;

Novikov_Thorne_Model NT_Model(r_in, r_out);

/*

Define some global boolians

*/

extern Const_bool lens_from_file = false;
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

    double metric[4][4]{}, N_obs, omega_obs;

        Spacetimes[e_metric]->get_metric(metric, &N_obs, &omega_obs, r_obs, theta_obs);

    double J, p_theta_0, p_r_0;

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
                    
                    Feed those initial conditions to the lenser

                    */

                Spacetimes[e_metric]->get_initial_conditions_from_file(&J, J_data, &p_theta_0, p_theta_data, &p_r_0, photon, r_obs,
                                                                       theta_obs, metric, N_obs, omega_obs);

                double initial_conditions[6] = { r_obs, theta_obs, phi_obs, J, p_theta_0, p_r_0 };

                Integration_status = Lens(initial_conditions, data, momentum_data);
    
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

                    get_intitial_conditions_from_angles(&J, &p_theta_0, &p_r_0, metric, V_angle, H_angle);

                    double initial_conditions[6] = { r_obs, theta_obs, phi_obs, J, p_theta_0, p_r_0 };

                    Integration_status = Lens(initial_conditions, data, momentum_data);

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