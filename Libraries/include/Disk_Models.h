#pragma once

#ifndef DISK_MODELS

	#define DISK_MODELS

    #include "Enumerations.h"
    #include "Spacetimes.h"

    #include <iostream>
    #include <cmath>
    #include <vector>

    extern Emission_model_enums e_emission;

    typedef class tag_Novikov_Thorne_Model {

        private:

            double r_in;
            double r_out;

        public:

            tag_Novikov_Thorne_Model(double x, double y);

            double get_r_in();
            double get_r_out();

            double Keplerian_angular_velocity(double r, std::vector<c_Spacetime_Base*> Spacetimes);

            double dr_Keplerian_angular_velocity(double r, std::vector<c_Spacetime_Base*> Spacetimes);

            double Redshift(double J, double State_Vector[], double r_obs, double theta_obs, std::vector<c_Spacetime_Base*> Spacetimes);

            double disk_Energy(double r, std::vector<c_Spacetime_Base*> Spacetimes);

            double disk_Angular_Momentum(double r, std::vector<c_Spacetime_Base*> Spacetimes);

            double Flux_integrand(double r, std::vector<c_Spacetime_Base*> Spacetimes);

            double solve_Flux_integral(double lower_bound, double upper_bound, double tolerance, std::vector<c_Spacetime_Base*> Spacetimes);

            double get_flux(double r, std::vector<c_Spacetime_Base*> Spacetimes);

    }Novikov_Thorne_Model;

    typedef class tag_Optically_Thin_Toroidal_Model {

        public:

            tag_Optically_Thin_Toroidal_Model();
     
            int get_disk_velocity(double Disk_velocity[], double State_vector[], std::vector<c_Spacetime_Base*> Spacetimes);

            double get_disk_density(double State_vector[]);

            double get_magnetic_field(double B_field[3], double State_vector[]);

            double get_emission_fucntion_synchotron_exact(double State_vector[], double J, std::vector<c_Spacetime_Base*> Spacetimes);

            double get_emission_fucntion_synchotron_phenomenological(double State_vector[], double J, std::vector<c_Spacetime_Base*> Spacetimes);

            double get_absorbtion_fucntion(double Emission_Function, double State_vector[], double J, double Frequency, double Temperature);

            double get_electron_pitch_angle(double State_vector[], double B_field_local[], std::vector<c_Spacetime_Base*> Spacetimes);

            double get_disk_temperature(double State_vector[]);

    }Optically_Thin_Toroidal_Model;

#endif

