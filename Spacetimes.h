#pragma once

#ifndef SPACETIMES
    
    #define SPACETIMES

    #include <iostream>
    #include <vector>
    #include "Enumerations.h"

    typedef class tag_Spacetime_Base_Class {

        public:

            virtual double get_ISCO(Orbit_Orientation Orientation) {

                std::cout << "Using Base Spacetime Class - Something Broke!" << '\n';

                return ERROR;

            };

            virtual double get_Photon_Sphere(Orbit_Orientation Orientation) {

                std::cout << "Using Base Spacetime Class - Something Broke!" << '\n';

                return ERROR;

            };

            /* Metric and its derivatives */

            virtual int get_metric(double metric[4][4], double* N_metric, 
                                   double* omega_metric, double r, double theta) {

                std::cout << "Using Base Spacetime Class - Something Broke!" << '\n';

                return ERROR;

            };

            virtual int get_dr_metric(double metric[4][4], double* N_metric,
                                      double* omega_metric, double r, double theta) {

                std::cout << "Using Base Spacetime Class - Something Broke!" << '\n';

                return ERROR;

            };

            virtual int get_d2r_metric(double metric[4][4], double* N_metric,
                                       double* omega_metric, double r, double theta) {

                std::cout << "Using Base Spacetime Class - Something Broke!" << '\n';

                return ERROR;

            };

            /* Initial conditions derived from images */

            virtual int get_initial_conditions_from_file(double* J, double J_data[], double* p_theta, double p_theta_data[], double* p_r,
                                                         int photon, double r_obs, double theta_obs, double metric[4][4], double N, double omega) {

                std::cout << "Using Base Spacetime Class - Something Broke!" << '\n';

                return ERROR;

            };

            /* Equations of motion */

            virtual int get_EOM(double inter_State_vector[], double J, double Derivatives[], int iteration) {

                std::cout << "Using Base Spacetime Class - Something Broke!" << '\n';
                
                return ERROR;
            
            };

            /* Integration Termination Conditions */

            virtual bool terminate_integration(double State_vector[], double Derivatives[]) { 

                std::cout << "Using Base Spacetime Class - Something Broke!" << '\n';

                return true; 
            
            };

    }c_Spacetime_Base;

    typedef class derived_Kerr_class : public c_Spacetime_Base{

        public:

            double get_ISCO(Orbit_Orientation Orientation);
            double get_Photon_Sphere(Orbit_Orientation Orientation);

            /* Metric and its derivatives */

            int get_metric(double metric[4][4], double* N_metric,
                           double* omega_metric, double r, double theta);

            int get_dr_metric(double metric[4][4], double* N_metric,
                              double* omega_metric, double r, double theta);

            int get_d2r_metric(double metric[4][4], double* N_metric,
                               double* omega_metric, double r, double theta);

            ///* Initial conditions derived from images */

            int get_initial_conditions_from_file(double* J, double J_data[], double* p_theta, double p_theta_data[], double* p_r,
                                                 int photon, double r_obs, double theta_obs, double metric[4][4], double N, double omega);

            /* Equations of motion */

            int get_EOM(double inter_State_vector[], double J, double Derivatives[], int iteration);

            /* Integration Termination Conditions */

            virtual bool terminate_integration(double State_vector[], double Derivatives[]);


    }Kerr_class;

    typedef class derived_Wormhole_class : public c_Spacetime_Base {

    public:

        double get_ISCO(Orbit_Orientation Orientation);
        double get_Photon_Sphere(Orbit_Orientation Orientation);

        /* Metric and its derivatives */

        int get_metric(double metric[4][4], double* N_metric,
                       double* omega_metric, double r, double theta);

        int get_dr_metric(double metric[4][4], double* N_metric,
                          double* omega_metric, double r, double theta);

        int get_d2r_metric(double metric[4][4], double* N_metric,
                           double* omega_metric, double r, double theta);

        /* Initial conditions derived from images */

        int get_initial_conditions_from_file(double* J, double J_data[], double* p_theta, double p_theta_data[], double* p_r,
                                             int photon, double r_obs, double theta_obs, double metric[4][4], double N, double omega);

        /* Equations of motion */

        int get_EOM(double inter_State_vector[], double J, double Derivatives[], int iteration);

        /* Integration Termination Conditions */

        bool terminate_integration(double State_vector[], double Derivatives[]);


    }Wormhole_class;

    typedef class derived_RBH_class : public c_Spacetime_Base {

    public:

        double get_ISCO(Orbit_Orientation Orientation);
        double get_Photon_Sphere(Orbit_Orientation Orientation);

        /* Metric and its derivatives */

        int get_metric(double metric[4][4], double* N_metric,
            double* omega_metric, double r, double theta);

        int get_dr_metric(double metric[4][4], double* N_metric,
                          double* omega_metric, double r, double theta);

        int get_d2r_metric(double metric[4][4], double* N_metric,
                           double* omega_metric, double r, double theta);

        /* Initial conditions derived from images */

        int get_initial_conditions_from_file(double* J, double J_data[], double* p_theta, double p_theta_data[], double* p_r,
            int photon, double r_obs, double theta_obs, double metric[4][4], double N, double omega);

        /* Equations of motion */

        int get_EOM(double inter_State_vector[], double J, double Derivatives[], int iteration);

        /* Integration Termination Conditions */

        bool terminate_integration(double State_vector[], double Derivatives[]);


    }RBH_class;

    typedef class derived_JNW_class : public c_Spacetime_Base {

    public:

        double get_ISCO(Orbit_Orientation Orientation);
        double get_Photon_Sphere(Orbit_Orientation Orientation);

        /* Metric and its derivatives */

        int get_metric(double metric[4][4], double* N_metric,
                       double* omega_metric, double r, double theta);

        int get_dr_metric(double metric[4][4], double* N_metric,
                          double* omega_metric, double r, double theta);

        int get_d2r_metric(double metric[4][4], double* N_metric,
                           double* omega_metric, double r, double theta);

        /* Initial conditions derived from images */

        int get_initial_conditions_from_file(double* J, double J_data[], double* p_theta, double p_theta_data[], double* p_r,
                                             int photon, double r_obs, double theta_obs, double metric[4][4], double N, double omega);

        /* Equations of motion */

        int get_EOM(double inter_State_vector[], double J, double Derivatives[], int iteration);

        /* Integration Termination Conditions */

        bool terminate_integration(double State_vector[], double Derivatives[]);


    }JNW_class;

    typedef class tag_observer {

    private:

        double r_obs;
        double theta_obs;
        double phi_obs;

    public:

        tag_observer(double r, double theta, double phi);

        double get_r_obs();
        double get_theta_obs();
        double get_phi_obs();

        int get_obs_velocity(double Obs_velocity[4], std::vector<c_Spacetime_Base*> Spacetimes);

    }c_Observer;

#endif
