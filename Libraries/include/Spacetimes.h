#pragma once

#ifndef SPACETIMES
    
    #define SPACETIMES

    #include <iostream>
    #include <vector>
    #include "Enumerations.h"

    class Spacetime_Base_Class {

        public:

            virtual double get_ISCO(Orbit_Orientation Orientation) {

                std::cout << "Using Base Spacetime Class - Something Broke!" << '\n';

                return ERROR;

            };

            virtual double* get_Photon_Sphere(Orbit_Orientation Orientation) {

                std::cout << "Using Base Spacetime Class - Something Broke!" << '\n';

                return NULL;

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

            virtual int get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon) {

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

    };

    class Kerr_class : public Spacetime_Base_Class {

        public:

            double get_ISCO(Orbit_Orientation Orientation) override;
            double* get_Photon_Sphere(Orbit_Orientation Orientation) override;

            /* Metric and its derivatives */

            int get_metric(double metric[4][4], double* N_metric,
                           double* omega_metric, double r, double theta) override;

            int get_dr_metric(double dr_metric[4][4], double* dr_N_metric,
                              double* dr_omega_metric, double r, double theta);

            int get_d2r_metric(double d2r_metric[4][4], double* d2r_N_metric,
                               double* d2r_omega_metric, double r, double theta);

            ///* Initial conditions derived from images */

            int get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon);

            /* Equations of motion */

            int get_EOM(double inter_State_vector[], double J, double Derivatives[], int iteration);

            /* Integration Termination Conditions */

            virtual bool terminate_integration(double State_vector[], double Derivatives[]);


    };

    class Wormhole_class : public Spacetime_Base_Class {

    public:

        double get_ISCO(Orbit_Orientation Orientation);
        double* get_Photon_Sphere(Orbit_Orientation Orientation);

        /* Metric and its derivatives */

        int get_metric(double metric[4][4], double* N_metric,
                       double* omega_metric, double r, double theta);

        int get_dr_metric(double dr_metric[4][4], double* dr_N_metric,
                          double* dr_omega_metric, double r, double theta);

        int get_d2r_metric(double d2r_metric[4][4], double* d2r_N_metric,
                           double* d2r_omega_metric, double r, double theta);

        /* Initial conditions derived from images */

        int get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon);

        /* Equations of motion */

        int get_EOM(double inter_State_vector[], double J, double Derivatives[], int iteration);

        /* Integration Termination Conditions */

        bool terminate_integration(double State_vector[], double Derivatives[]);


    };

    class RBH_class : public Spacetime_Base_Class {

    public:

        double get_ISCO(Orbit_Orientation Orientation);
        double* get_Photon_Sphere(Orbit_Orientation Orientation);

        /* Metric and its derivatives */

        int get_metric(double metric[4][4], double* N_metric,
                       double* omega_metric, double r, double theta);

        int get_dr_metric(double dr_metric[4][4], double* dr_N_metric,
                          double* dr_omega_metric, double r, double theta);

        int get_d2r_metric(double d2r_metric[4][4], double* d2r_N_metric,
                           double* d2r_omega_metric, double r, double theta);

        /* Initial conditions derived from images */

        int get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon);

        /* Equations of motion */

        int get_EOM(double inter_State_vector[], double J, double Derivatives[], int iteration);

        /* Integration Termination Conditions */

        bool terminate_integration(double State_vector[], double Derivatives[]);


    };

    class JNW_class : public Spacetime_Base_Class {

    public:

        double get_ISCO(Orbit_Orientation Orientation);
        double* get_Photon_Sphere(Orbit_Orientation Orientation);

        /* Metric and its derivatives */

        int get_metric(double metric[4][4], double* N_metric,
                       double* omega_metric, double r, double theta);

        int get_dr_metric(double dr_metric[4][4], double* dr_N_metric,
                          double* dr_omega_metric, double r, double theta);

        int get_d2r_metric(double d2r_metric[4][4], double* d2r_N_metric,
                           double* d2r_omega_metric, double r, double theta);

        /* Initial conditions derived from images */

        int get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon);

        /* Equations of motion */

        int get_EOM(double inter_State_vector[], double J, double Derivatives[], int iteration);

        /* Integration Termination Conditions */

        bool terminate_integration(double State_vector[], double Derivatives[]);


    };

    class Gauss_Bonnet_class : public Spacetime_Base_Class {

    public:

        double get_ISCO(Orbit_Orientation Orientation);
        double* get_Photon_Sphere(Orbit_Orientation Orientation);

        /* Metric and its derivatives */

        int get_metric(double metric[4][4], double* N_metric,
                       double* omega_metric, double r, double theta);

        int get_dr_metric(double dr_metric[4][4], double* dr_N_metric,
                          double* dr_omega_metric, double r, double theta);

        int get_d2r_metric(double d2r_metric[4][4], double* d2r_N_metric,
                           double* d2r_omega_metric, double r, double theta);

        /* Initial conditions derived from images */

        int get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon);

        /* Equations of motion */

        int get_EOM(double inter_State_vector[], double J, double Derivatives[], int iteration);

        /* Integration Termination Conditions */

        bool terminate_integration(double State_vector[], double Derivatives[]);

    };

    class Black_Hole_w_Dark_Matter_Halo_class : public Spacetime_Base_Class {

    public:

        /* Metric and its derivatives */

        int get_metric(double metric[4][4], double* N_metric,
            double* omega_metric, double r, double theta);

        int get_dr_metric(double dr_metric[4][4], double* dr_N_metric,
            double* dr_omega_metric, double r, double theta);


        /* Initial conditions derived from images */

 
        /* Equations of motion */

        int get_EOM(double inter_State_vector[], double J, double Derivatives[], int iteration);

        /* Integration Termination Conditions */

        bool terminate_integration(double State_vector[], double Derivatives[]);

    };

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

        int get_obs_velocity(double Obs_velocity[4], std::vector<Spacetime_Base_Class*> Spacetimes);

    }c_Observer;

#endif
