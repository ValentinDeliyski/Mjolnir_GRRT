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

            virtual Metric_type get_metric(const double* const State_Vector) {

                std::cout << "Using Base Spacetime Class - Something Broke!" << '\n';

                return {};

            };

            virtual Metric_type get_dr_metric(const double* const State_Vector) {

                std::cout << "Using Base Spacetime Class - Something Broke!" << '\n';

                return {};

            };


            virtual Metric_type get_dtheta_metric(const double* const State_Vector) {

                std::cout << "Using Base Spacetime Class - Something Broke!" << '\n';

                return {};

            };

            virtual Metric_type get_d2r_metric(const double* const State_Vector) {

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

            virtual bool load_parameters(Metric_parameters_type Metric_Parameters) {
            
                std::cout << "Using Base Spacetime Class - Something Broke!" << '\n';

                return true;
            
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
            
            double Mass = 1.0;
            double Spin_Param{};

        public:

            double* get_ISCO() override;
            double* get_Photon_Sphere() override;

            /* Metric and its derivatives */

            Metric_type get_metric(const double* const State_Vector) override;
            Metric_type get_dr_metric(const double* const State_Vector) override;
            Metric_type get_dtheta_metric(const double* const State_Vector) override;
            Metric_type get_d2r_metric(const double* const State_Vector) override;

            /* Initial conditions derived from images */

            int get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon) override;

            /* Equations of motion */

            int get_EOM(double inter_State_vector[], double Derivatives[]) override;

            /* Integration Termination Conditions */

            bool terminate_integration(double State_vector[], double Derivatives[]) override;

            bool load_parameters(Metric_parameters_type Metric_Parameters);

            void update_parameters(double Param_value, Metric_Parameter_Selector Parameter);

    };

    class Wormhole_class : public Spacetime_Base_Class {

    private:

        Metric_type s_Metric{};
        Metric_type s_dr_Metric{};
        Metric_type s_dtheta_Metric{};
        Metric_type s_d2r_Metric{};

        double Mass = 1.0;
        double R_Throat = this->Mass;

        double Spin_Param;
        double Redshift_Param;
        bool Stop_at_Throat;

    public:


        double* get_ISCO();
        double* get_Photon_Sphere();

        /* Metric and its derivatives */

        Metric_type get_metric(const double* const State_Vector) override;
        Metric_type get_dr_metric(const double* const State_Vector) override;
        Metric_type get_dtheta_metric(const double* const State_Vector) override;
        Metric_type get_d2r_metric(const double* const State_Vector) override;

        /* Initial conditions derived from images */

        int get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon) override;

        /* Equations of motion */

        int get_EOM(double inter_State_vector[], double Derivatives[]) override;

        /* Integration Termination Conditions */

        bool terminate_integration(double State_vector[], double Derivatives[]) override;

        bool load_parameters(Metric_parameters_type Metric_Parameters);

        void update_parameters(double Param_value, Metric_Parameter_Selector Parameter);

    };

    class RBH_class : public Spacetime_Base_Class {

    private:

        Metric_type s_Metric{};
        Metric_type s_dr_Metric{};
        Metric_type s_dtheta_Metric{};
        Metric_type s_d2r_Metric{};

        double Mass = 1.0;
        double Parameter;

    public:

        double* get_ISCO();
        double* get_Photon_Sphere();

        /* Metric and its derivatives */

        Metric_type get_metric(const double* const State_Vector) override;
        Metric_type get_dr_metric(const double* const State_Vector) override;
        Metric_type get_dtheta_metric(const double* const State_Vector) override;
        Metric_type get_d2r_metric(const double* const State_Vector) override;

        /* Initial conditions derived from images */

        int get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon) override;

        /* Equations of motion */

        int get_EOM(double inter_State_vector[], double Derivatives[]) override;

        /* Integration Termination Conditions */

        bool terminate_integration(double State_vector[], double Derivatives[]) override;

        bool load_parameters(Metric_parameters_type Metric_Parameters);

        void update_parameters(double Param_value, Metric_Parameter_Selector Parameter);

    };

    class JNW_class : public Spacetime_Base_Class {

    private:

        Metric_type s_Metric{};
        Metric_type s_dr_Metric{};
        Metric_type s_dtheta_Metric{};
        Metric_type s_d2r_Metric{};

        double Mass = 1.0;
        double Gamma;

    public:

        double* get_ISCO();
        double* get_Photon_Sphere();

        /* Metric and its derivatives */

        Metric_type get_metric(const double* const State_Vector) override;
        Metric_type get_dr_metric(const double* const State_Vector) override;
        Metric_type get_dtheta_metric(const double* const State_Vector) override;
        Metric_type get_d2r_metric(const double* const State_Vector) override;

        /* Initial conditions derived from images */

        int get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon) override;

        /* Equations of motion */

        int get_EOM(double inter_State_vector[], double Derivatives[]) override;

        /* Integration Termination Conditions */

        bool terminate_integration(double State_vector[], double Derivatives[]) override;

        bool load_parameters(Metric_parameters_type Metric_Parameters);

        void update_parameters(double Param_value, Metric_Parameter_Selector Parameter);

    };

    class Gauss_Bonnet_class : public Spacetime_Base_Class {

    private:

        Metric_type s_Metric{};
        Metric_type s_dr_Metric{};
        Metric_type s_dtheta_Metric{};
        Metric_type s_d2r_Metric{};

        double Mass = 1.0;
        double Gamma;

    public:

        double* get_ISCO() ;
        double* get_Photon_Sphere() ;

        /* Metric and its derivatives */

        Metric_type get_metric(const double* const State_Vector) override;
        Metric_type get_dr_metric(const double* const State_Vector) override;
        Metric_type get_dtheta_metric(const double* const State_Vector) override;
        Metric_type get_d2r_metric(const double* const State_Vector) override;

        /* Initial conditions derived from images */

        int get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon) override;

        /* Equations of motion */

        int get_EOM(double inter_State_vector[], double Derivatives[]) override;

        /* Integration Termination Conditions */

        bool terminate_integration(double State_vector[], double Derivatives[]) override;

        bool load_parameters(Metric_parameters_type Metric_Parameters);

        void update_parameters(double Param_value, Metric_Parameter_Selector Parameter);

    };

    class Black_Hole_w_Dark_Matter_Halo_class : public Spacetime_Base_Class {

    private:

        Metric_type s_Metric{};
        Metric_type s_dr_Metric{};
        Metric_type s_dtheta_Metric{};
        Metric_type s_d2r_Metric{};

        double Mass = 1.0;
        double Compactness;
        double Halo_Mass;

    public:

        double* get_ISCO();

        /* Metric and its derivatives */

        Metric_type get_metric(const double* const State_Vector) override;
        Metric_type get_dr_metric(const double* const State_Vector) override;
        Metric_type get_dtheta_metric(const double* const State_Vector) override;

        /* Initial conditions derived from images */

        int get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon) override;
 
        /* Equations of motion */

        int get_EOM(double inter_State_vector[], double Derivatives[]) override;

        /* Integration Termination Conditions */

        bool terminate_integration(double State_vector[], double Derivatives[]) override;

        bool load_parameters(Metric_parameters_type Metric_Parameters);

        void update_parameters(double Param_value, Metric_Parameter_Selector Parameter);

    };

    class Observer_class {

    private:

        Observer_parameters_type obs_params;
        double obs_velocity[4];

    public:

        Observer_class(Simulation_Context_type* p_Sim_Context);

        Observer_parameters_type get_parameters();

        double* get_obs_velocity();

    };

#endif
