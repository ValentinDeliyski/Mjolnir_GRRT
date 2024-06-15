#pragma once

#ifndef GENERAL_FUNCTIONS

    #define _USE_MATH_DEFINES

    #define GENERAL_FUNCTIONS

    struct Initial_conditions_type;
    struct Metric_type;

    double get_planck_function_CGS(double Frequency, double Temperature);

    void invert_metric(double inv_metric[4][4], double metric[4][4]);

    double get_metric_det(double metric[4][4]);

    double get_eq_induced_metric_det(double metric[4][4]);

    void get_intitial_conditions_from_angles(Initial_conditions_type* p_Initial_Conditions, double V_angle, double H_angle);

    void get_impact_parameters(Initial_conditions_type* p_Initial_Conditions, double Image_coords[2]);

    double Redshift(double State_Vector[], double U_source[]);

    void get_ZAMO_tetrad(double e_t[4], double e_r[4], double e_theta[4], double e_phi[4], double metric[4][4]);

    void Contravariant_coord_to_ZAMO(double metric[4][4], double Contravariant_Vector[4], double ZAMO_Vector[4]);

    void Lorentz_boost_matrix(double Boost_matrix[4][4], double U_source[4]);

    int Increment_theta_turning_points(double State_Vector[], double Old_State[]);

    int compute_image_order(int N_theta_turning_points, Initial_conditions_type* p_Initial_Conditions);

    void get_connection_coefficients(Metric_type s_Metric, Metric_type s_dr_metric, Metric_type s_dtheta_metric, double Connectrion_Coeffs[4][4][4]);

#endif 