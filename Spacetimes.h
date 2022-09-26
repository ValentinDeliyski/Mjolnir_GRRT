#pragma once

#ifndef SPACETIMES
    
    #define SPACETIMES

    #include "Enumerations.h"

    typedef class tag_Kerr {

	    private:

		    double a;
		    double M;

		    double r_ph_retrograde;
		    double r_ph_prograde;
		    double r_horizon;
		    double r_ISCO;

	    public:

            tag_Kerr(double x);

            double get_spin();
            double get_ISCO();
            double get_r_horizon();;
            double get_r_ph_prograde();
            double get_r_ph_retrograde();

            int metric(double metric[4][4], double* N_metric, double* omega_metric, double r, double theta);

            int metric_first_derivatives(class tag_Kerr Kerr_class, double dr_metric[4][4], double* dr_N, double* dr_omega,
                double r, double theta);

            int metric_second_derivatives(class tag_Kerr Kerr_class, double d2r_metric[4][4], double* d2r_N, double* d2r_omega,
                double r, double theta);

            int intitial_conditions_from_file(double* J, double J_data[], double* p_theta, double p_theta_data[], double* p_r,
                int photon, double r_obs, double theta_obs, double metric[4][4]);

            int EOM(double inter_State_vector[], double J, double Derivatives[], int iteration);


    }c_Kerr;

    typedef class tag_Wormhole {

	    private:

		    double alpha_metric;
		    double r_throat;
            double a;
		    double M;

		    double r_ISCO;
		    double r_ph;

	    public:

            tag_Wormhole(double x, double spin);

            double get_metric_parameter();
            double get_r_throat();
            double get_ISCO();
            double get_r_ph();
            double get_spin();

            int metric(double metric[4][4], double* N_metric, double* omega, double l, double theta);

            int metric_first_derivatives(class tag_Wormhole Wormhole_class, double dr_metric[4][4], double* dr_N, double* dr_omega,
                double l, double theta);

            int metric_second_derivatives(class tag_Wormhole Wormhole_class, double d2r_metric[4][4], double* d2r_N, double* d2r_omega,
                double l, double theta);

            int intitial_conditions_from_file(double* J, double J_data[], double* p_theta, double p_theta_data[], double* p_r,
                int photon, double r_obs, double theta_obs, double metric[4][4], double N, double omega);

            int EOM(double inter_State_vector[], double J, double Derivatives[], int iteration);

    }c_Wormhole;

    typedef class tag_Regular_Black_Hole {

	    private:

		    double metric_parameter;
		    double M;

		    double r_horizon;
		    double r_ISCO;
		    double r_ph;

	    public:

            tag_Regular_Black_Hole(double x);

            double get_metric_parameter();
            double get_r_horizon();
            double get_ISCO();
            double get_r_ph();

            int metric(double metric[4][4], double* N_metric, double* omega_metric,
                double r, double theta);

            int metric_first_derivatives(class tag_Regular_Black_Hole RBH_class, double dr_metric[4][4], double* dr_N,
                double* dr_omega, double r, double theta);

            int metric_second_derivatives(class tag_Regular_Black_Hole RBH_class, double d2r_metric[4][4], double* d2r_N,
                double* d2r_omega, double r, double theta);

            int intitial_conditions_from_file(double* J, double J_data[], double* p_theta, double p_theta_data[], double* p_r,
                int photon, double r_obs, double theta_obs, double metric[4][4]);

            int EOM(double inter_State_vector[], double J, double Derivatives[], int iteration);

    }c_RBH;

    typedef class tag_JNW_Naked_Singularity {

    private:

        double M;
        double GAMMA;
        double r_SINGULARITY;

        double r_ph;
        double r_ISCO_outer;
        double r_ISCO_inner;

    public:

        tag_JNW_Naked_Singularity(double gamma, double b);

        double get_r_ph();
        double get_gamma();
        double get_r_ISCO_outer();
        double get_r_ISCO_inner();
        double get_r_singularity();

        int metric(double metric[4][4], double* N_metric, double* omega_metric, double r, double theta);

        int metric_first_derivatives(class tag_JNW_Naked_Singularity JNW_class, double d2r_metric[4][4],
            double* d2r_N, double* d2r_omega, double r, double theta);

        int metric_second_derivatives(class tag_JNW_Naked_Singularity JNW_class, double d2r_metric[4][4],
            double* d2r_N, double* d2r_omega, double r, double theta);

        int intitial_conditions_from_file(double* J, double J_data[], double* p_theta, double p_theta_data[], double* p_r,
            int photon, double r_obs, double theta_obs, double metric[4][4]);

        int EOM(double inter_State_vector[], double J, double Derivatives[], int iteration);

    }c_JNW_Naked_Singularity;

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

        int get_obs_velocity(double Obs_velocity[4],
                             e_Spacetimes e_metric, c_Kerr Kerr_class, c_Wormhole Wormhole_class, c_RBH RBH_class, c_JNW_Naked_Singularity JNW_class);

    }c_Observer;

   
#endif
