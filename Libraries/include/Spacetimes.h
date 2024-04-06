#pragma once

#ifndef SPACETIMES
    
    #define SPACETIMES

    #include <iostream>
    #include "Structs.h"

    class Spacetime_Base_Class {

        public:

            virtual double* get_ISCO() {

                std::cout << "Using Base Spacetime Class - Something Broke!" << '\n';

                return NULL;

            };

            virtual double* get_Photon_Sphere() {

                std::cout << "Using Base Spacetime Class - Something Broke!" << '\n';

                return NULL;

            };

            /* Metric and its derivatives */

            virtual Metric_type get_metric(double State_vector[]) {

                std::cout << "Using Base Spacetime Class - Something Broke!" << '\n';

                return {};

            };

            virtual Metric_type get_dr_metric(double State_vector[]) {

                std::cout << "Using Base Spacetime Class - Something Broke!" << '\n';

                return {};

            };


            virtual Metric_type get_dtheta_metric(double State_vector[]) {

                std::cout << "Using Base Spacetime Class - Something Broke!" << '\n';

                return {};

            };

            virtual Metric_type get_d2r_metric(double State_vector[]) {

                std::cout << "Using Base Spacetime Class - Something Broke!" << '\n';

                return {};

            };

            /* Initial conditions derived from images */

            virtual int get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon) {

                std::cout << "Using Base Spacetime Class - Something Broke!" << '\n';

                return ERROR;

            };

            /* Equations of motion */

            virtual int get_EOM(double inter_State_vector[], double Derivatives[]) {

                std::cout << "Using Base Spacetime Class - Something Broke!" << '\n';
                
                return ERROR;
            
            };

            /* Integration Termination Conditions */

            virtual bool terminate_integration(double State_vector[], double Derivatives[]) { 

                std::cout << "Using Base Spacetime Class - Something Broke!" << '\n';

                return true; 
            
            };

            virtual bool load_parameters(Metric_Parameters_type Metric_Parameters) {
            
                std::cout << "Using Base Spacetime Class - Something Broke!" << '\n';

                return true;
            
            };

            virtual Metric_Parameters_type get_parameters() {
            
                std::cout << "Using Base Spacetime Class - Something Broke!" << '\n';

                Metric_Parameters_type Error{};

                return Error;
            
            };

            virtual void update_parameters(double Param_value, Metric_Parameter_Selector Parameter) {
            
                std::cout << "Using Base Spacetime Class - Something Broke!" << '\n';
            
            };

    };

    class Kerr_class : public Spacetime_Base_Class {

        private:

            Metric_type s_Metric;
            Metric_type s_dr_Metric;
            Metric_type s_dtheta_Metric{};
            Metric_type s_d2r_Metric;
            
            double Mass = 1.0f;
            double Spin_Param{};

        public:

            double* get_ISCO() override;
            double* get_Photon_Sphere() override;

            /* Metric and its derivatives */

            int update_metric(double State_vector[]);
            Metric_type get_metric(double State_vector[]) override;

            int update_dr_metric(double State_vector[]);
            Metric_type get_dr_metric(double State_vector[]) override;

            int update_dtheta_metric(double State_vector[]);
            Metric_type get_dtheta_metric(double State_vector[]) override;

            int update_d2r_metric(double State_vector[]);
            Metric_type get_d2r_metric(double State_vector[]) override;

            /* Initial conditions derived from images */

            int get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon) override;

            /* Equations of motion */

            int get_EOM(double inter_State_vector[], double Derivatives[]) override;

            /* Integration Termination Conditions */

            bool terminate_integration(double State_vector[], double Derivatives[]) override;

            bool load_parameters(Metric_Parameters_type Metric_Parameters);

            Metric_Parameters_type get_parameters();

            void update_parameters(double Param_value, Metric_Parameter_Selector Parameter);

    };

    class Wormhole_class : public Spacetime_Base_Class {

    private:

        Metric_type s_Metric{};
        Metric_type s_dr_Metric{};
        Metric_type s_dtheta_Metric{};
        Metric_type s_d2r_Metric{};

        double Mass = 1.0f;
        double R_Throat = this->Mass;

        double Spin_Param{};
        double Redshift_Param{};
        bool Stop_at_Throat;

    public:


        double* get_ISCO();
        double* get_Photon_Sphere();

        /* Metric and its derivatives */

        int update_metric(double State_vector[]);
        Metric_type get_metric(double State_vector[]) override;

        int update_dr_metric(double State_vector[]);
        Metric_type get_dr_metric(double State_vector[]) override;

        int update_dtheta_metric(double State_vector[]);
        Metric_type get_dtheta_metric(double State_vector[]) override;

        int update_d2r_metric(double State_vector[]);
        Metric_type get_d2r_metric(double State_vector[]) override;

        /* Initial conditions derived from images */

        int get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon) override;

        /* Equations of motion */

        int get_EOM(double inter_State_vector[], double Derivatives[]) override;

        /* Integration Termination Conditions */

        bool terminate_integration(double State_vector[], double Derivatives[]) override;

        bool load_parameters(Metric_Parameters_type Metric_Parameters);

        Metric_Parameters_type get_parameters();

        void update_parameters(double Param_value, Metric_Parameter_Selector Parameter);

    };

    class RBH_class : public Spacetime_Base_Class {

    private:

        Metric_type s_Metric{};
        Metric_type s_dr_Metric{};
        Metric_type s_dtheta_Metric{};
        Metric_type s_d2r_Metric{};

        double Parameter;

    public:

        double* get_ISCO();
        double* get_Photon_Sphere();

        /* Metric and its derivatives */

        int update_metric(double State_vector[]);
        Metric_type get_metric(double State_vector[]) override;

        int update_dr_metric(double State_vector[]);
        Metric_type get_dr_metric(double State_vector[]) override;

        int update_dtheta_metric(double State_vector[]);
        Metric_type get_dtheta_metric(double State_vector[]) override;

        int update_d2r_metric(double State_vector[]);
        Metric_type get_d2r_metric(double State_vector[]) override;

        /* Initial conditions derived from images */

        int get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon) override;

        /* Equations of motion */

        int get_EOM(double inter_State_vector[], double Derivatives[]) override;

        /* Integration Termination Conditions */

        bool terminate_integration(double State_vector[], double Derivatives[]) override;

        bool load_parameters(Metric_Parameters_type Metric_Parameters);

        Metric_Parameters_type get_parameters();

        void update_parameters(double Param_value, Metric_Parameter_Selector Parameter);

    };

    class JNW_class : public Spacetime_Base_Class {

    private:

        Metric_type s_Metric{};
        Metric_type s_dr_Metric{};
        Metric_type s_dtheta_Metric{};
        Metric_type s_d2r_Metric{};

        double Gamma;

    public:

        double* get_ISCO();
        double* get_Photon_Sphere();

        /* Metric and its derivatives */

        int update_metric(double State_vector[]);
        Metric_type get_metric(double State_vector[]) override;

        int update_dr_metric(double State_vector[]);
        Metric_type get_dr_metric(double State_vector[]) override;

        int update_dtheta_metric(double State_vector[]);
        Metric_type get_dtheta_metric(double State_vector[]) override;

        int update_d2r_metric(double State_vector[]);
        Metric_type get_d2r_metric(double State_vector[]) override;

        /* Initial conditions derived from images */

        int get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon) override;

        /* Equations of motion */

        int get_EOM(double inter_State_vector[], double Derivatives[]) override;

        /* Integration Termination Conditions */

        bool terminate_integration(double State_vector[], double Derivatives[]) override;

        bool load_parameters(Metric_Parameters_type Metric_Parameters);

        Metric_Parameters_type get_parameters();

        void update_parameters(double Param_value, Metric_Parameter_Selector Parameter);

    };

    class Gauss_Bonnet_class : public Spacetime_Base_Class {

    private:

        Metric_type s_Metric{};
        Metric_type s_dr_Metric{};
        Metric_type s_dtheta_Metric{};
        Metric_type s_d2r_Metric{};

        double Gamma;

    public:

        double* get_ISCO() ;
        double* get_Photon_Sphere() ;

        /* Metric and its derivatives */

        int update_metric(double State_vector[]);
        Metric_type get_metric(double State_vector[]) override;

        int update_dr_metric(double State_vector[]);
        Metric_type get_dr_metric(double State_vector[]) override;

        int update_dtheta_metric(double State_vector[]);
        Metric_type get_dtheta_metric(double State_vector[]) override;

        int update_d2r_metric(double State_vector[]);
        Metric_type get_d2r_metric(double State_vector[]) override;

        /* Initial conditions derived from images */

        int get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon) override;

        /* Equations of motion */

        int get_EOM(double inter_State_vector[], double Derivatives[]) override;

        /* Integration Termination Conditions */

        bool terminate_integration(double State_vector[], double Derivatives[]) override;

        bool load_parameters(Metric_Parameters_type Metric_Parameters);

        Metric_Parameters_type get_parameters();

        void update_parameters(double Param_value, Metric_Parameter_Selector Parameter);

    };

    class Black_Hole_w_Dark_Matter_Halo_class : public Spacetime_Base_Class {

    private:

        Metric_type s_Metric{};
        Metric_type s_dr_Metric{};
        Metric_type s_dtheta_Metric{};
        Metric_type s_d2r_Metric{};

        double Compactness;
        double Halo_Mass;

    public:

        double* get_ISCO();

        /* Metric and its derivatives */

        int update_metric(double State_vector[]);
        Metric_type get_metric(double State_vector[]) override;

        int update_dr_metric(double State_vector[]);
        Metric_type get_dr_metric(double State_vector[]) override;

        int update_dtheta_metric(double State_vector[]);
        Metric_type get_dtheta_metric(double State_vector[]) override;

        /* Initial conditions derived from images */

        int get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon) override;
 
        /* Equations of motion */

        int get_EOM(double inter_State_vector[], double Derivatives[]) override;

        /* Integration Termination Conditions */

        bool terminate_integration(double State_vector[], double Derivatives[]) override;

        bool load_parameters(Metric_Parameters_type Metric_Parameters);

        Metric_Parameters_type get_parameters();

        void update_parameters(double Param_value, Metric_Parameter_Selector Parameter);

    };

    class Observer_class {

    private:

        double r_obs;
        double theta_obs;
        double phi_obs;
        double obs_velocity[4];

    public:

        Observer_class(double r, double theta, double phi);

        double get_r_obs();
        double get_theta_obs();
        double get_phi_obs();

        int get_obs_velocity(double Obs_velocity[4]);

    };

#endif
