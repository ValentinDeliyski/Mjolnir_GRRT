/***************************************************************************************************
|                                                                                                  |
|                      ---------  Gravitational Ray Tracer  ---------                              | 
|                                                                                                  |
|    * @Version: 3.1                                                                               |
|    * @Author: Valentin Deliyski                                                                  |
|    * @Description: This program numeriaclly integrates the equations of motion                   |
|    for null geodesics in a curved spacetime and and projects them onto an observer's screen      |
|                                                                                                  |
|    * @Supported Spacetimes:                                                                      | 
|        ** Kerr Black Holes                                                                       |
|        ** Static Regular Black Holes                                                             |
|        ** Rotating Traversable Wormholes                                                         |
|                                                                                                  |
|    * @Supported Disk Models                                                                      |
|        ** Novikov-Thorne                                                                         |
|        ** Generic Optically Thin Disk With Arbitrary Density, Emission and Absorbtion Profiles   |
|                                                                                                  |
***************************************************************************************************/

#define _USE_MATH_DEFINES

#include <string>

#include <iostream>
#include <iomanip> 
#include <fstream>
#include <cmath>

#include "Enumerations.h"
#include "Constants.h"
#include "Spacetimes.h"
#include "IO_files.h"

#include "Disk_Models.h"
#include "General_functions.h"

#include "Lensing.h"

e_Spacetimes e_metric = Kerr;
Disk_Models e_Disk_Model = Optically_Thin_Toroidal;

int main() {

    bool lens_from_file = false;
    bool truncate       = false;

    double r_obs, theta_obs, phi_obs;

        r_obs = 10'000;
        theta_obs = 85. / 180 * M_PI;
        phi_obs = 0;

    double M, a, r_throat, metric_parameter;

        M = 1.0;
        metric_parameter = 0.0;
        a = 0.98;

    /*
    Define classes that hold the spacetime properites
    */

    c_Kerr Kerr_class(a);
    c_RBH RBH_class(metric_parameter);
    c_Wormhole Wormhole_class(metric_parameter,a);

        r_throat = Wormhole_class.get_r_throat();

    /*
    Create/Open the logging files
    */
        
    std::ofstream data[4], momentum_data[4];

        open_output_files(e_metric, data, momentum_data, truncate);

    /*
    Get the ISCO orbits from the spacetime classes and set the inner Novikov-Thorne disk radius
    */

    double r_in, r_out, r_ISCO;

        switch (e_metric) {

            case Kerr:

                r_ISCO = Kerr_class.get_ISCO();
                    
                r_in  = 6;
                r_out = 50;

                break;

            case Reg_Black_Hole:

                r_ISCO = RBH_class.get_ISCO();

                r_in  = 4.5;
                r_out = 50;

                break;

            case Wormhole:

                r_ISCO = Wormhole_class.get_ISCO();

                r_in  = 1;
                r_out = 50 * r_ISCO;

                break;

            default:

                std::cout << "Wrong metric!" << '\n';

                return ERROR;

        }


    Novikov_Thorne_Model NT_Model(r_in, r_out);

    /*
    Set the Optically Thin Toroidal Disk model parameters
    */

    double disk_alpha, disk_height_scale, disk_rad_cutoff, disk_omega, disk_magnetization;

        disk_alpha = 3;
        disk_height_scale = 0.1;
        disk_rad_cutoff = 4 * M;
        disk_omega = sqrt(1. / 12) * M;
        disk_magnetization = 0.01;

    Optically_Thin_Toroidal_Model OTT_Model(disk_alpha, disk_height_scale, disk_rad_cutoff, disk_omega);
    
    /*
    Get the metric at the observer to feed into the initial conditions functions
    */

    double metric[4][4], N_obs, omega_obs;

        get_metric(e_metric, metric, &N_obs, &omega_obs, r_obs, theta_obs,
                   Kerr_class, RBH_class, Wormhole_class);

    double J, p_theta_0, p_r_0;

    Return_Value_enums Integration_status = OK;

    print_ASCII_art();

    std::cout << "Observer Radial Position [GM/c^2] = " << r_obs << '\n';
    std::cout << "Observer Inclination [deg] = "   << int(theta_obs / M_PI * 180) << '\n';

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

                get_initial_conditions_from_file(e_metric, &J, J_data, &p_theta_0, p_theta_data, &p_r_0, photon, r_obs, 
                                                 theta_obs, metric, N_obs, omega_obs,
                                                 Kerr_class, RBH_class, Wormhole_class);

                double initial_conditions[6] = { r_obs, theta_obs, phi_obs, J, p_theta_0, p_r_0 };

                Integration_status = Lens(initial_conditions, M, metric_parameter, a, r_throat, r_in, r_out,
                                          lens_from_file, data, momentum_data, e_metric, Kerr_class, RBH_class, Wormhole_class,
                                          e_Disk_Model, NT_Model, OTT_Model);

    
                print_progress(photon, Data_number, lens_from_file);
            }

            std::cout << '\n';
        }
    }
    else{

        /*
        Setup a viewing window for the observer and loop trough it
        */

        double V_angle_min = -0.0015;
        double V_angle_max = 0.0025;

        double H_angle_min = -0.0055;
        double H_angle_max = 0.0055;

        double Scan_Step = 5e-6;
   
        int progress = 0;
        
        int V_num = floor(log10f((V_angle_max - V_angle_min) * 10000) + 1);

        if (Integration_status == OK) {

            for (double V_angle = 0; V_angle >= -V_angle_max; V_angle -= Scan_Step) {

                print_progress(progress, int((V_angle_max - V_angle_min) / Scan_Step), lens_from_file);

                progress += 1;

                for (double H_angle = H_angle_min; H_angle <= H_angle_max; H_angle += Scan_Step) {

                    get_intitial_conditions_from_angles(&J, &p_theta_0, &p_r_0, metric, V_angle, H_angle);

                    double initial_conditions[6] = { r_obs, theta_obs, phi_obs, J, p_theta_0, p_r_0 };

                    Integration_status = Lens(initial_conditions, M, metric_parameter, a, r_throat, r_in, r_out,
                                              lens_from_file, data, momentum_data, e_metric, Kerr_class, RBH_class, Wormhole_class,
                                              e_Disk_Model, NT_Model, OTT_Model);

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