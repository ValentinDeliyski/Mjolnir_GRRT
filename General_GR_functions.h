#pragma once

#ifndef GENERAL_FUNCTIONS

    #define GENERAL_FUNCTIONS

    #include "Enumerations.h"
    #include "Spacetimes.h"
    #include "Disk_Models.h"
    #include <vector>

    double get_planck_function_CGS(double Frequency, double Temperature);

    int invert_metric(double inv_metric[4][4], double metric[4][4]);

    double get_metric_det(double metric[4][4]);

    int get_intitial_conditions_from_angles(Initial_conditions_type* p_Initial_Conditions, double V_angle, double H_angle);

    double Redshift(double J, double State_Vector[], double U_source[]);

    double get_photon_t_momentum(double State_vector[], double J, double metric[4][4]);

    int get_ZAMO_tetrad(double e_t[4], double e_r[4], double e_theta[4], double e_phi[4], double metric[4][4]);

    int Contravariant_coord_to_ZAMO(double metric[4][4], double Contravariant_Vector[4], double ZAMO_Vector[4]);

    int Lorentz_boost_matrix(double Boost_matrix[4][4], double U_source[4], double metric[4][4]);

    int get_Radiative_Transfer(double State_Vector[], double Derivatives[], int iteration, double J);

#endif