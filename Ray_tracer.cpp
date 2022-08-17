/***************************************************************************************************
|                                                                                                  |
|                          ---------  Numerical Ray Tracer  ---------                              | 
|                                                                                                  |
|    * @Version: 3.0                                                                               |
|    * @Author: Valentin Deliyski                                                                  |
|    * @Description: This program numeriaclly integrates the equations of motion                   |
|    for null geodesics in a curved spacetime and and projects them onto an observer's screen      |
|                                                                                                  |
|    * @Supported Spacetimes:                                                                      | 
|        ** Kerr Black Holes                                                                       |
|        ** Static Regular Black Holes                                                             |
|        ** Rotating Traversable Wormholes                                                         |
|                                                                                                  |
***************************************************************************************************/

#include "Enums_Constants.h"

#include <iostream>
#include <iomanip> 
#include <fstream>
#include <math.h>

#include "Spacetimes.h"
#include "IO_files.h"
#include "General_functions.h"
#include "Lensing.h"

Spacetimes e_metric = Kerr;

int main() {

    bool lens_from_file = false;

    double max_error     = RK45_ACCURACY;
    double safety_factor = SAFETY;

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
    c_Wormhole Wormhole_class(metric_parameter);

        r_throat = Wormhole_class.get_r_throat();

    /*
    Create/Open the logging files
    */
        
    std::ofstream data[4], momentum_data[4];

        open_output_files(e_metric, data, momentum_data);

    /*
    Get the ISCO orbits from the spacetime classes and set the inner disc radius
    */

    double r_in, r_out, r_ISCO;

        switch (e_metric) {

            case Kerr:

                r_ISCO = Kerr_class.get_r_ISCO();
                    
                r_in  = 4.5;
                r_out = 50;

                break;

            case Reg_Black_Hole:

                r_ISCO = RBH_class.get_r_ISCO();

                r_in  = 3.2;
                r_out = 50;

                break;

            case Wormhole:

                r_ISCO = Wormhole_class.get_r_ISCO();

                r_in  = r_ISCO;
                r_out = 50 * r_ISCO;

                break;

            default:

                std::cout << "Wrong metric!" << '\n';

                return -1;

            }
    
    /*
    Get the metric at the observer to feed into the initial conditions functions
    */

    double metric[4][4], N_obs, omega_obs;

        get_metric(e_metric, metric, &N_obs, &omega_obs, M, r_throat, a, metric_parameter, r_obs, theta_obs,
                   Kerr_class, RBH_class, Wormhole_class);

    double J, p_theta_0, p_r_0;

    if (lens_from_file) {

        /*
        Read the initial conditions from file
        */

        double J_data[500]{}, p_theta_data[500]{};
        int Data_number = 0;

            get_geodesic_data(J_data, p_theta_data, &Data_number);

        for (int photon = 0; photon <= Data_number; photon += 1) {

           /*
           Feed those initial conditions to the lenser
           */

            get_initial_conditions_from_file(e_metric, &J, J_data, &p_theta_0, p_theta_data, &p_r_0,
                                             photon, r_obs, theta_obs, metric, N_obs, omega_obs, M ,a, metric_parameter,
                                             Kerr_class, RBH_class, Wormhole_class);

            double initial_conditions[6] = { r_obs, theta_obs, phi_obs, J, p_theta_0, p_r_0 };

            Lens(initial_conditions, M, metric_parameter, a, r_throat, Coeff_deriv, Coeff_sol, Coeff_test_sol,
                 max_error, safety_factor, r_in, r_out, lens_from_file,data, momentum_data, e_metric,
                 Kerr_class, RBH_class, Wormhole_class);

        }

    }
    else{

        /*
        Setup a viewing window for the observer and loop trough it
        */

        double V_angle_max = 0.0017;
        double H_angle_max = 0.006;

        for (double V_angle = -0.0015; V_angle <= 0.0025; V_angle += 5e-6) {

            std::cout << std::fixed << std::setprecision(6) << V_angle << " ";

            for (double H_angle = -0.0055; H_angle <= 0.0055; H_angle += 5e-6) {

                get_intitial_conditions_from_angles(&J, &p_theta_0, &p_r_0, metric, V_angle, H_angle);

                double initial_conditions[6] = { r_obs, theta_obs, phi_obs, J, p_theta_0, p_r_0 };

                Lens(initial_conditions, M, metric_parameter, a, r_throat, Coeff_deriv, Coeff_sol, Coeff_test_sol,
                     max_error, safety_factor, r_in, r_out, lens_from_file,data, momentum_data, e_metric,
                     Kerr_class, RBH_class, Wormhole_class);

            }

        }

    }            

    close_output_files(data, momentum_data);

    return 0;
}