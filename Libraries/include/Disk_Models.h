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

            Novikov_Thorne_Model(double x, double y, Spacetime_Base_Class* Spacetime);

            double get_r_in();
            double get_r_out();

            double Keplerian_angular_velocity(double r, Spacetime_Base_Class* Spacetime);

            double dr_Keplerian_angular_velocity(double r, Spacetime_Base_Class* Spacetime);

            double Redshift(double J, double State_Vector[], double r_obs, double theta_obs, Spacetime_Base_Class* Spacetime);

            double disk_Energy(double r, Spacetime_Base_Class* Spacetime);

            double disk_Angular_Momentum(double r, Spacetime_Base_Class* Spacetime);

            double Flux_integrand(double r, Spacetime_Base_Class* Spacetime);

            double solve_Flux_integral(double lower_bound, double upper_bound, double tolerance, Spacetime_Base_Class* Spacetime);

            double get_flux(double r, Spacetime_Base_Class* Spacetime);

    };

    class Generic_Optically_Thin_Model {

        private:
            Disk_model_parameters s_Disk_Parameters{};
            Emission_law_parameters s_Emission_Parameters{};
            Precomputed_e_pitch_angles s_Precomputed_e_pitch_angles{};
            double Disk_velocity[4]{};
            double Disk_Temperature{};
            double Disk_density_profile{};
            double Hotspot_density_profile{};


            // ====================== Thermally Distributed Synchotron Fit Functions ====================== //

            void get_thermal_synchotron_emission_fit_functions(Thermal_Syncotron_fit_selector e_fit_functions,
                                                              double Emission_functions[4],
                                                              double X,
                                                              double sqrt_X,
                                                              double cbrt_X);

            void get_thermal_synchotron_faradey_fit_functions(double X,
                                                              double X_to_1_point_2,
                                                              double X_frac, 
                                                              double faradey_fucntions[STOKES_PARAM_NUM]);

            void get_thermal_synchotron_fit_functions(double Emission_fucntions[STOKES_PARAM_NUM],
                                                      double Faradey_functions[STOKES_PARAM_NUM],
                                                      Thermal_emission_f_arguments* Emission_args,
                                                      Thermal_faradey_f_arguments* Faradey_args);


            // ====================== Kappa Distributed Synchotron Fit Functions ====================== //

            void get_kappa_synchotron_emission_fit_functions(double Emission_functions[STOKES_PARAM_NUM],
                                                             double X,
                                                             double sqrt_X,
                                                             double cbrt_X,
                                                             double X_to_7_over_20,
                                                             double kappa,
                                                             double sin_emission_angle,
                                                             double T_electron_dim);

            void get_kappa_synchotron_absorbtion_fit_functions(double Absorbtion_functions[STOKES_PARAM_NUM],
                                                               double X,
                                                               double sqrt_X,
                                                               double cbrt_X,
                                                               double X_to_7_over_20,
                                                               double kappa,
                                                               double sin_emission_angle,
                                                               double T_electron_dim);

            // ====================== Power Distributed Synchotron Fit Functions ====================== //


        public:

            // ====================== Accretion Disk State Functions ====================== //

            double  get_disk_temperature(double State_vector[]);

            double* get_disk_velocity(double State_Vector[], Simulation_Context_type* p_Sim_Context);

            double  get_disk_density_profile(double State_vector[]);

            double  get_disk_hotspot_density_profile(double State_Vector[]);

            double get_magnetic_field(double B_field[4], 
                                      double State_Vector[],
                                      Simulation_Context_type* p_Sim_Context);

            // ======================= Radiative Fransfer Functions ======================= //

            void get_radiative_transfer_functions(double State_Vector[e_State_Number],
                                                  Simulation_Context_type* p_Sim_Context,
                                                  double Emission_functions[STOKES_PARAM_NUM],
                                                  double Absorbtion_functions[STOKES_PARAM_NUM],
                                                  double Faradey_functions[STOKES_PARAM_NUM]);

            void get_thermal_synchotron_transfer_functions(double State_vector[],
                                                           Simulation_Context_type* p_Sim_Context,
                                                           double Emission_functions[STOKES_PARAM_NUM],
                                                           double Faradey_functions[STOKES_PARAM_NUM],
                                                           double Absorbtion_functions[STOKES_PARAM_NUM]);

            void get_phenomenological_synchotron_functions(double State_Vector[],
                                                           Simulation_Context_type* p_Sim_Context, 
                                                           double Emission_functions[STOKES_PARAM_NUM],
                                                           double Faradey_functions[STOKES_PARAM_NUM],
                                                           double Absorbtion_functions[STOKES_PARAM_NUM]);

            // ======================= Electron Pitch Angle Functions ======================= //

            double get_electron_pitch_angle(double B_field_local[4], 
                                            double State_Vector[], 
                                            Simulation_Context_type* p_Sim_Context);

            // =============================== Misc Functions =============================== //

            void precompute_electron_pitch_angles();

            int load_parameters(Disk_model_parameters* p_Disk_Parameters, 
                                Emission_law_parameters* p_Emission_Parameters);

    };

#endif

