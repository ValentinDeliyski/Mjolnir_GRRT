#pragma once

#ifndef DISK_MODELS

	#define DISK_MODELS

    typedef class tag_Novikov_Thorne_Model {

        private:

            double r_in;
            double r_out;

        public:

            tag_Novikov_Thorne_Model(double x, double y);

            double get_r_in();
            double get_r_out();

            double Keplerian_angular_velocity(e_Spacetimes e_metric, double r,
                c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class, c_JNW_Naked_Singularity JNW_class);

            double dr_Keplerian_angular_velocity(e_Spacetimes e_metric, double r, double Kepler,
                c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class, c_JNW_Naked_Singularity JNW_class);

            double Redshift(e_Spacetimes e_metric, double J, double State_Vector[], double r_obs, double theta_obs,
                c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class, c_JNW_Naked_Singularity JNW_class);

            double disk_Energy(e_Spacetimes e_metric, double r,
                c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class, c_JNW_Naked_Singularity JNW_class);

            double disk_Angular_Momentum(e_Spacetimes e_metric, double r,
                c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class, c_JNW_Naked_Singularity JNW_class);

            double Flux_integrand(e_Spacetimes e_metric, double r,
                c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class, c_JNW_Naked_Singularity JNW_class);

            double solve_Flux_integral(e_Spacetimes e_metric, double lower_bound, double upper_bound, double tolerance,
                c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class, c_JNW_Naked_Singularity JNW_class);

            double get_flux(e_Spacetimes e_metric, double r, double r_in,
                c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class, c_JNW_Naked_Singularity JNW_class);

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

            int get_disk_velocity(double Disk_velocity[], double State_vector[], e_Spacetimes e_metric,
                                  c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class, c_JNW_Naked_Singularity JNW_class);

            double get_disk_density(double State_vector[]);

            double get_magnetic_field(double B_field[3], double State_vector[]);

    }Optically_Thin_Toroidal_Model;

#endif

