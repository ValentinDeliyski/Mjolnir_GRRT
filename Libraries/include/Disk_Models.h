#pragma once

#ifndef DISK_MODELS

    #define DISK_MODELS

    #include "Spacetimes.h"

    #include <iostream>
    #include <vector>

    class Novikov_Thorne_Model {

        private:

            double r_in;
            double r_out;

        public:

            Novikov_Thorne_Model(double x, double y);

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

        public:

            double HOTSPOT_R_COORD;
            double HOTSPOT_PHI_COORD;

            Optically_Thin_Toroidal_Model(double hotspot_radial_position ,double hotspot_phi_angle);
     
            int get_disk_velocity(double Disk_velocity[], double State_vector[], Spacetime_Base_Class* Spacetimes[]);

            double get_disk_hotspot(double State_Vector[]);

            double get_disk_density_profile(double State_vector[]);

            double get_magnetic_field(double B_field[3], double State_vector[]);

            double get_emission_function_synchotron_exact(double State_vector[], Spacetime_Base_Class* Spacetimes[]);

            double get_emission_function_synchotron_phenomenological(double State_vector[], Spacetime_Base_Class* Spacetimes[]);

            double get_absorbtion_function(double Emission_Function, double State_vector[], double redshift, double Frequency, double Temperature);

            double get_electron_pitch_angle(double State_vector[], double B_field_local[], Spacetime_Base_Class* Spacetimes[]);

            double get_disk_temperature(double State_vector[]);

    };

#endif

