#pragma once

#ifndef DISK_MODELS

    #define DISK_MODELS

    #include "Structs.h"

    class Spacetime_Base_Class;

    class Novikov_Thorne_Model {

        private:

            double r_in;
            double r_out;

        public:

            Novikov_Thorne_Model(double x, double y, Spacetime_Base_Class* Spacetimes[]);

            double get_r_in();
            double get_r_out();

            double Keplerian_angular_velocity(double r, Spacetime_Base_Class* Spacetimes[]);

            double dr_Keplerian_angular_velocity(double r, Spacetime_Base_Class* Spacetimes[]);

            double Redshift(double J, double State_Vector[], double r_obs, double theta_obs, Spacetime_Base_Class* Spacetimes[]);

            double disk_Energy(double r, Spacetime_Base_Class* Spacetimes[]);

            double disk_Angular_Momentum(double r, Spacetime_Base_Class* Spacetimes[]);

            double Flux_integrand(double r, Spacetime_Base_Class* Spacetimes[]);

            double solve_Flux_integral(double lower_bound, double upper_bound, double tolerance, Spacetime_Base_Class* Spacetimes[]);

            double get_flux(double r, Spacetime_Base_Class* Spacetimes[]);

    };

    class Optically_Thin_Toroidal_Model {

        private:
            Disk_model_parameters s_Disk_Parameters{};
            Emission_law_parameters s_Emission_Parameters{};
            Precomputed_e_pitch_angles s_Precomputed_e_pitch_angles{};
            double Disk_velocity[4]{};
            double Disk_Temperature{};
            double Disk_density_profile{};
            double Disk_hotspot_density{};

            int update_disk_velocity(double State_vector[], Initial_conditions_type* s_Initial_Conditions);
            int update_disk_temperature(double State_vector[]);
            int update_disk_density_profile(double State_vector[]);
            int update_disk_hotspot(double State_vector[]);

            void get_synchotron_emission_fit_function(Sync_emission_fit_functions e_fit_functions,
                                                         double Emission_functions[4],
                                                         double X,
                                                         double X_1_2,
                                                         double X_1_3);

        public:

            Optically_Thin_Toroidal_Model(Disk_model_parameters* p_Disk_Parameters, Emission_law_parameters* p_Emission_Parameters);

            double  get_disk_temperature(double State_vector[]);
            double* get_disk_velocity(double State_vector[], Initial_conditions_type* s_Initial_Conditions);
            double  get_disk_density_profile(double State_vector[]);
            double  get_disk_hotspot(double State_Vector[]);

            double get_magnetic_field(double B_field[3], double State_vector[]);

            void get_emission_function_synchotron_exact(double State_vector[], 
                                                        Initial_conditions_type* s_Initial_Conditions, 
                                                        double Emission_fucntions[4]);

            void get_emission_function_synchotron_phenomenological(double State_vector[], Initial_conditions_type* s_Initial_Conditions, double Emission_functions[4]);

            void get_absorbtion_function(double Emission_Functions[4], double State_vector[], double redshift, double Frequency, double Absorbtion_functions[4]);

            double get_electron_pitch_angle(double State_vector[], double B_field_local[], Initial_conditions_type* s_Initial_Conditions);

            void precompute_electron_pitch_angles();

            int load_parameters(Disk_model_parameters* p_Disk_Parameters, Emission_law_parameters* p_Emission_Parameters);

            Disk_model_parameters get_disk_params();

    };

#endif

