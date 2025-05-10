#pragma once
#define _USE_MATH_DEFINES
#include "Structs.h"
#include <iostream>
#include <cmath>

class Spacetime_Base_Class {

public:

    virtual double* get_ISCO() {

        std::cout << "Using Base Spacetime Class - Something Broke!" << '\n';

        return NULL;

    };

    virtual double* get_Photon_Sphere() {

        std::cout << "Using Base Spacetime Class - Something Broke!" << '\n';

        return NULL;

    };

    /* Metric and its derivatives */

    virtual Metric_type get_metric(const double* const State_Vector) const  {

        std::cout << "Using Base Spacetime Class - Something Broke!" << '\n';

        return {};

    };

    virtual Metric_type get_dr_metric(const double* const State_Vector) const  {

        std::cout << "Using Base Spacetime Class - Something Broke!" << '\n';

        return {};

    };

    virtual Metric_type get_dtheta_metric(const double* const State_Vector) const  {

        std::cout << "Using Base Spacetime Class - Something Broke!" << '\n';

        return {};

    };

    virtual Metric_type get_d2r_metric(const double* const State_Vector) const  {

        std::cout << "Using Base Spacetime Class - Something Broke!" << '\n';

        return {};

    };

    /* Initial conditions derived from images */

    virtual int get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon) {

        std::cout << "Using Base Spacetime Class - Something Broke!" << '\n';

        return ERROR;

    };

    /* Equations of motion */

    virtual void get_EOM(double State_Vector[], double Derivatives[]) const {

        /* ----------- Temporary matrix, used to store intermediate calculations ----------- */
        double temp_matrix[4][4]{};

        Metric_type Metric = this->get_metric(State_Vector);
        Metric_type dr_Metric = this->get_dr_metric(State_Vector);
        Metric_type dtheta_Metric = this->get_dtheta_metric(State_Vector);

        double inv_metric[4][4]{};
        invert_metric(inv_metric, Metric.Metric);

        double dr_inv_metric[4][4]{};
        matrix_matrix_multiply(dr_Metric.Metric, inv_metric, temp_matrix);
        matrix_matrix_multiply(inv_metric, temp_matrix, dr_inv_metric);

        double dtheta_inv_metric[4][4]{};
        matrix_matrix_multiply(dtheta_Metric.Metric, inv_metric, temp_matrix);
        matrix_matrix_multiply(inv_metric, temp_matrix, dtheta_inv_metric);

        /* ----------- There is a minus sign infront of the whole expression for the derivative of an inverse of a matrix ----------- */
        for (int left_idx = 0; left_idx <= 3; left_idx++) {

            for (int right_idx = 0; right_idx <= 3; right_idx++) {

                dr_inv_metric[left_idx][right_idx] *= -1;
                dtheta_inv_metric[left_idx][right_idx] *= -1;

            }

        }

        for (int right_idx = 0; right_idx <= 3; right_idx++) {

            *(Derivatives + e_t) += inv_metric[e_t][right_idx] * State_Vector[right_idx + 4];
            *(Derivatives + e_r) += inv_metric[e_r][right_idx] * State_Vector[right_idx + 4];
            *(Derivatives + e_theta) += inv_metric[e_theta][right_idx] * State_Vector[right_idx + 4];
            *(Derivatives + e_phi) += inv_metric[e_phi][right_idx] * State_Vector[right_idx + 4];

            for (int left_idx = 0; left_idx <= 3; left_idx++) {

                *(Derivatives + e_p_r) += -1. / 2 * dr_inv_metric[left_idx][right_idx] * State_Vector[left_idx + 4] * State_Vector[right_idx + 4];
                *(Derivatives + e_p_theta) += -1. / 2 * dtheta_inv_metric[left_idx][right_idx] * State_Vector[left_idx + 4] * State_Vector[right_idx + 4];

            }

        }

        *(Derivatives + e_p_t) = 0.0;
        *(Derivatives + e_p_phi) = 0.0;
    
    };

    /* Integration Termination Conditions */

    virtual bool terminate_integration(double State_vector[], double Derivatives[]) { 

        std::cout << "Using Base Spacetime Class - Something Broke!" << '\n';

        return true; 
    
    };

    virtual Return_Values load_parameters(const Metric_parameters_type* const Metric_Parameters) {
    
        std::cout << "Using Base Spacetime Class - Something Broke!" << '\n';

        return ERROR;
    
    };

};

class Kerr_class : public Spacetime_Base_Class {

private:

    double Mass = 1.0;
    double Spin_Param;

public:

    double* get_ISCO() override;
    double* get_Photon_Sphere() override;

    /* Metric and its derivatives */

    Metric_type get_metric(const double* const State_Vector) const override;
    Metric_type get_dr_metric(const double* const State_Vector) const override;
    Metric_type get_dtheta_metric(const double* const State_Vector) const override;
    Metric_type get_d2r_metric(const double* const State_Vector) const override;

    /* Initial conditions derived from images */

    int get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon) override;

    /* Equations of motion */

    void get_EOM(double inter_State_vector[], double Derivatives[]) const override;

    /* Integration Termination Conditions */

    bool terminate_integration(double State_vector[], double Derivatives[]) override;

    Return_Values load_parameters(const Metric_parameters_type* const Metric_Parameters) override;
     
};

class Wormhole_class : public Spacetime_Base_Class {

private:

    double Mass = 1.0;
    double R_Throat = this->Mass;
    double Spin_Param;
    double Redshift_Param;

    bool Stop_at_Throat;

public:


    double* get_ISCO();
    double* get_Photon_Sphere();

    /* Metric and its derivatives */

    Metric_type get_metric(const double* const State_Vector) const override;
    Metric_type get_dr_metric(const double* const State_Vector) const override;
    Metric_type get_dtheta_metric(const double* const State_Vector) const override;
    Metric_type get_d2r_metric(const double* const State_Vector) const override;

    /* Initial conditions derived from images */

    int get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon) override;

    /* Equations of motion */

    void get_EOM(double inter_State_vector[], double Derivatives[]) const override;

    /* Integration Termination Conditions */

    bool terminate_integration(double State_vector[], double Derivatives[]) override;

    Return_Values load_parameters(const Metric_parameters_type* const Metric_Parameters) override;

};

class RBH_class : public Spacetime_Base_Class {

private:

    double Mass = 1.0;
    double Parameter;

public:

    double* get_ISCO();
    double* get_Photon_Sphere();

    /* Metric and its derivatives */

    Metric_type get_metric(const double* const State_Vector) const override;
    Metric_type get_dr_metric(const double* const State_Vector) const override;
    Metric_type get_dtheta_metric(const double* const State_Vector) const override;
    Metric_type get_d2r_metric(const double* const State_Vector) const override;

    /* Initial conditions derived from images */

    int get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon) override;

    /* Equations of motion */

    void get_EOM(double inter_State_vector[], double Derivatives[]) const override;

    /* Integration Termination Conditions */

    bool terminate_integration(double State_vector[], double Derivatives[]) override;

    Return_Values load_parameters(const Metric_parameters_type* const Metric_Parameters) override;

};

class JNW_class : public Spacetime_Base_Class {

private:

    double Mass = 1.0;
    double Gamma;

public:

    double* get_ISCO();
    double* get_Photon_Sphere();

    /* Metric and its derivatives */

    Metric_type get_metric(const double* const State_Vector) const override;
    Metric_type get_dr_metric(const double* const State_Vector) const override;
    Metric_type get_dtheta_metric(const double* const State_Vector) const override;
    Metric_type get_d2r_metric(const double* const State_Vector) const override;

    /* Initial conditions derived from images */

    int get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon) override;

    /* Equations of motion */

    void get_EOM(double inter_State_vector[], double Derivatives[]) const override;

    /* Integration Termination Conditions */

    bool terminate_integration(double State_vector[], double Derivatives[]) override;

    Return_Values load_parameters(const Metric_parameters_type* const Metric_Parameters) override;

};

class Gauss_Bonnet_class : public Spacetime_Base_Class {

private:;

    double Mass = 1.0;
    double Gamma;

public:

    double* get_ISCO() ;
    double* get_Photon_Sphere() ;

    /* Metric and its derivatives */

    Metric_type get_metric(const double* const State_Vector) const override;
    Metric_type get_dr_metric(const double* const State_Vector) const override;
    Metric_type get_dtheta_metric(const double* const State_Vector) const override;
    Metric_type get_d2r_metric(const double* const State_Vector) const override;

    /* Initial conditions derived from images */

    int get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon) override;

    /* Equations of motion */

    void get_EOM(double inter_State_vector[], double Derivatives[]) const override;

    /* Integration Termination Conditions */

    bool terminate_integration(double State_vector[], double Derivatives[]) override;

    Return_Values load_parameters(const Metric_parameters_type* const Metric_Parameters) override;

};

class Black_Hole_w_Dark_Matter_Halo_class : public Spacetime_Base_Class {

private:

    double Mass = 1.0;
    double Compactness;
    double Halo_Mass;

public:

    double* get_ISCO();

    /* Metric and its derivatives */

    Metric_type get_metric(const double* const State_Vector) const override;
    Metric_type get_dr_metric(const double* const State_Vector) const override;
    Metric_type get_dtheta_metric(const double* const State_Vector) const override;

    /* Initial conditions derived from images */

    int get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon) override;

    /* Equations of motion */

    void get_EOM(double inter_State_vector[], double Derivatives[]) const override;

    /* Integration Termination Conditions */

    bool terminate_integration(double State_vector[], double Derivatives[]) override;

    Return_Values load_parameters(const Metric_parameters_type* const Metric_Parameters) override;

};

class Numerical_metric : public Spacetime_Base_Class {

private:

    Numerical_metric_params_type Parameters;

    inline void get_control_point_matrix(const double* const Control_vector, const int r_idx, const int theta_idx, double Control_matrix[4][4]) const;

    Metric_type comute_Minkowski_metric(const double* const State_Vector) const;
    Metric_type comute_dr_Minkowski_metric(const double* const State_Vector) const;
    Metric_type comute_d2r_Minkowski_metric(const double* const State_Vector) const;
    Metric_type comute_dtheta_Minkowski_metric(const double* const State_Vector) const;

    void get_polynomial_basis_vector(const double natural_parameter, double* const Polynomial_basis_vector) const;
    void get_derivative_polynomial_basis_vector(const double natural_parameter, double* const Polynomial_basis_vector) const;
    void get_second_derivative_polynomial_basis_vector(const double natural_parameter, double* const Polynomial_basis_vector) const;

    double evaluate_single_spline(const double Control_point_matrix[4][4], const double Radial_natural_parameter, const double Theta_natural_parameter, Derivative_selector_enums Derivative_selector) const;
    Metric_type evaluate_all_splines(const double* const State_Vector, Derivative_selector_enums Derivative_selector) const;
    Metric_type compute_metric_components_from_spline(const double* const State_vector, Derivative_selector_enums Derivative_selector) const;

    double compactify_radial_coordiante(const double r) const;
    double uncompactify_radial_coordinate(const double x);

public:

    /* Metric and its derivatives */

    Metric_type get_metric(const double* const State_Vector) const override;
    Metric_type get_dr_metric(const double* const State_Vector) const override;
    Metric_type get_dtheta_metric(const double* const State_Vector) const override;
    Metric_type get_d2r_metric(const double* const State_Vector) const override;

    /* Initial conditions derived from images */

    int get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon) override;

    /* Equations of motion */

    void get_EOM(double inter_State_vector[], double Derivatives[]) const override;

    /* Integration Termination Conditions */

    bool terminate_integration(double State_vector[], double Derivatives[]) override;

    Return_Values load_parameters(const Metric_parameters_type* const Metric_Parameters) override;

};

class Observer_class {

private:

    Observer_parameters_type obs_params;
    double obs_velocity[4];

public:

    Observer_class(Simulation_Context_type* p_Sim_Context);

    Observer_parameters_type get_parameters() const;

    const double* get_obs_velocity() const;

};