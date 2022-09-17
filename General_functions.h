#pragma once

#ifndef GENERAL_FUNCTIONS

    #define GENERAL_FUNCTIONS

    #include "Enumerations.h"
    #include "Spacetimes.h"
    #include "Disk_Models.h"

    int Rorate_to_obs_plane(double theta_obs, double phi_obs, double Image_point[3], double rotated_Image_point[3]);

    double my_max(double vector[]);

    bool crossed_equatior(double State_vector[], double Old_State_Vector[]);

    double dot_product(double vector_1[3], double vector_2[3]);

    int invert_metric(double inv_metric[4][4], double metric[4][4]);

    double get_ISCO(e_Spacetimes e_metric, c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class);

    int get_metric(e_Spacetimes e_metric, double metric[4][4], double* N_metric, double* omega, double r, double theta,
        c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class);

    double get_metric_det(double metric[4][4]);

    int get_metric_fist_derivatives(e_Spacetimes e_metric, double dr_metric[4][4], double* dr_N_metric, double* dr_omega,
        double r, double theta, c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class);

    int get_metric_second_derivatives(e_Spacetimes e_metric, double d2r_metric[4][4], double* d2r_N_metric, double* d2r_omega,
        double r, double theta, c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class);

    int get_intitial_conditions_from_angles(double* J, double* p_theta, double* p_r, double metric[4][4],
        double V_angle, double H_angle);

    int get_initial_conditions_from_file(e_Spacetimes e_metric, double* J, double J_data[], double* p_theta, double p_theta_data[],
        double* p_r, int photon, double r_obs, double theta_obs, double metric[4][4],
        double N_metric, double omega_metric,
        c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class);

    double Redshift(e_Spacetimes e_metric, Disk_Models Disk_Model, double J, double State_Vector[], double r_obs, double theta_obs,
        c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class);

    int Loretnz_boost();

    int get_Radiative_Transfer(double State_Vector[], double Derivatives[], int iteration, Optically_Thin_Toroidal_Model OTT_Model,
        double redshift);

    int get_EOM(e_Spacetimes e_metric, double inter_State_vector[7 * 6], double J, double Derivatives[7 * 6], int iteration,
        c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class);

    void print_ASCII_art();

    void print_progress(int current, int max, bool lens_from_file);

    Disk_Intersection Disk_event(Disk_Models e_Disk_Model, double State_Vector[], double Old_State_Vector[],
        Novikov_Thorne_Model NT_Model, Optically_Thin_Toroidal_Model OTT_Model);

#endif