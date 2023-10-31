#pragma once

#ifndef SPACETIMES
    
    #define SPACETIMES

    #include <iostream>
    #include <vector>
    #include "Enumerations.h"
    #include "Structs.h"

    class Spacetime_Base_Class {

        private:

            Metric_type s_Metric;
            Metric_type s_dr_Metric;
            Metric_type s_d2r_Metric;
            bool eval_bitmask[3]{};
            bool ignore_flag{};

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

            virtual int get_EOM(double inter_State_vector[], double Derivatives[], int iteration) {

                std::cout << "Using Base Spacetime Class - Something Broke!" << '\n';
                
                return ERROR;
            
            };

            /* Integration Termination Conditions */

            virtual bool terminate_integration(double State_vector[], double Derivatives[]) { 

                std::cout << "Using Base Spacetime Class - Something Broke!" << '\n';

                return true; 
            
            };

            virtual void set_ignore_flag(bool flag);

            virtual void reset_eval_bitmask();

    };

    class Kerr_class : public Spacetime_Base_Class {

        private:

            Metric_type s_Metric;
            Metric_type s_dr_Metric;
            Metric_type s_d2r_Metric;
            bool eval_bitmask[3]{};
            bool ignore_flag{};

        public:

            double* get_ISCO() override;
            double* get_Photon_Sphere() override;

            /* Metric and its derivatives */

            int update_metric(double State_vector[]);
            Metric_type get_metric(double State_vector[]);

            int update_dr_metric(double State_vector[]);
            Metric_type get_dr_metric(double State_vector[]);

            int update_d2r_metric(double State_vector[]);
            Metric_type get_d2r_metric(double State_vector[]);

            /* Initial conditions derived from images */

            int get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon);

            /* Equations of motion */

            int get_EOM(double inter_State_vector[], double Derivatives[], int iteration);

            /* Integration Termination Conditions */

            bool terminate_integration(double State_vector[], double Derivatives[]);
            void set_ignore_flag(bool flag);
            void reset_eval_bitmask();

    };

    class Wormhole_class : public Spacetime_Base_Class {

    private:

        Metric_type s_Metric{};
        Metric_type s_dr_Metric{};
        Metric_type s_d2r_Metric{};
        bool eval_bitmask[3]{};
        bool ignore_flag{};

    public:


        double* get_ISCO();
        double* get_Photon_Sphere();

        /* Metric and its derivatives */

        int update_metric(double State_vector[]);
        Metric_type get_metric(double State_vector[]);

        int update_dr_metric(double State_vector[]);
        Metric_type get_dr_metric(double State_vector[]);

        int update_d2r_metric(double State_vector[]);
        Metric_type get_d2r_metric(double State_vector[]);

        /* Initial conditions derived from images */

        int get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon);

        /* Equations of motion */

        int get_EOM(double inter_State_vector[], double Derivatives[], int iteration);

        /* Integration Termination Conditions */

        bool terminate_integration(double State_vector[], double Derivatives[]);
        void set_ignore_flag(bool flag);
        void reset_eval_bitmask();

    };

    class RBH_class : public Spacetime_Base_Class {

    private:

        Metric_type s_Metric{};
        Metric_type s_dr_Metric{};
        Metric_type s_d2r_Metric{};
        bool eval_bitmask[3]{};
        bool ignore_flag{};

    public:

        double* get_ISCO();
        double* get_Photon_Sphere();

        /* Metric and its derivatives */

        int update_metric(double State_vector[]);
        Metric_type get_metric(double State_vector[]);

        int update_dr_metric(double State_vector[]);
        Metric_type get_dr_metric(double State_vector[]);

        int update_d2r_metric(double State_vector[]);
        Metric_type get_d2r_metric(double State_vector[]);

        /* Initial conditions derived from images */

        int get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon);

        /* Equations of motion */

        int get_EOM(double inter_State_vector[], double Derivatives[], int iteration);

        /* Integration Termination Conditions */

        bool terminate_integration(double State_vector[], double Derivatives[]);
        void set_ignore_flag(bool flag);
        void reset_eval_bitmask();

    };

    class JNW_class : public Spacetime_Base_Class {

    private:

        Metric_type s_Metric{};
        Metric_type s_dr_Metric{};
        Metric_type s_d2r_Metric{};
        bool eval_bitmask[3] = {false, false, false};
        bool ignore_flag{};

    public:

        double* get_ISCO();
        double* get_Photon_Sphere();

        /* Metric and its derivatives */

        int update_metric(double State_vector[]);
        Metric_type get_metric(double State_vector[]);

        int update_dr_metric(double State_vector[]);
        Metric_type get_dr_metric(double State_vector[]);

        int update_d2r_metric(double State_vector[]);
        Metric_type get_d2r_metric(double State_vector[]);

        /* Initial conditions derived from images */

        int get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon);

        /* Equations of motion */

        int get_EOM(double inter_State_vector[], double Derivatives[], int iteration);

        /* Integration Termination Conditions */

        bool terminate_integration(double State_vector[], double Derivatives[]);

        void set_ignore_flag(bool flag);
        void reset_eval_bitmask();

    };

    class Gauss_Bonnet_class : public Spacetime_Base_Class {

    private:

        Metric_type s_Metric{};
        Metric_type s_dr_Metric{};
        Metric_type s_d2r_Metric{};
        bool eval_bitmask[3]{};
        bool ignore_flag{};

    public:

        double* get_ISCO() ;
        double* get_Photon_Sphere() ;

        /* Metric and its derivatives */

        int update_metric(double State_vector[]);
        Metric_type get_metric(double State_vector[]);

        int update_dr_metric(double State_vector[]);
        Metric_type get_dr_metric(double State_vector[]);

        int update_d2r_metric(double State_vector[]);
        Metric_type get_d2r_metric(double State_vector[]);

        /* Initial conditions derived from images */

        int get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon);

        /* Equations of motion */

        int get_EOM(double inter_State_vector[], double Derivatives[], int iteration);

        /* Integration Termination Conditions */

        bool terminate_integration(double State_vector[], double Derivatives[]);
        void set_ignore_flag(bool flag);
        void reset_eval_bitmask();
    };

    class Black_Hole_w_Dark_Matter_Halo_class : public Spacetime_Base_Class {

    private:

        Metric_type s_Metric{};
        Metric_type s_dr_Metric{};
        Metric_type s_d2r_Metric{};
        bool eval_bitmask[3]{};
        bool ignore_flag{};

    public:

        double* get_ISCO();

        /* Metric and its derivatives */

        int update_metric(double State_vector[]);
        Metric_type get_metric(double State_vector[]);

        int update_dr_metric(double State_vector[]);
        Metric_type get_dr_metric(double State_vector[]);

        /* Initial conditions derived from images */

        int get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon);
 
        /* Equations of motion */

        int get_EOM(double inter_State_vector[], double Derivatives[], int iteration);

        /* Integration Termination Conditions */

        bool terminate_integration(double State_vector[], double Derivatives[]);
        void set_ignore_flag(bool flag);
        void reset_eval_bitmask();

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
