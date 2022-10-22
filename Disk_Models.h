#pragma once

#ifndef DISK_MODELS

    #include <vector>

	#define DISK_MODELS

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

        private:

            double DISK_ALPHA;
            double DISK_HEIGHT_SCALE;
            double DISK_RAD_CUTOFF;
            double DISK_OMEGA;
            double DISK_MAGNETIZATION;
            double MAG_FIELD_GEOMETRY[3];

        public:

            tag_Optically_Thin_Toroidal_Model(double alpha, double height_scale, double rad_cutoff, double omega,
                                              double magnetization, double mag_field[3]);

            double get_disk_alpha();
            double get_disk_height_scale();
            double get_disk_rad_cutoff();
            double get_disk_omega();
            double get_disk_magnetization();

            int get_disk_velocity(double Disk_velocity[], double State_vector[], std::vector<c_Spacetime_Base*> Spacetimes);

            double get_disk_density(double State_vector[]);

            double get_magnetic_field(double B_field[3], double State_vector[]);

    }Optically_Thin_Toroidal_Model;

#endif

