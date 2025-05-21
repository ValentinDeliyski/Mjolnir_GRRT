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

        std::cout << "Using Base Spacetime Class - Something Broke!'\n'";
    
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

    Numerical_metric_potentials_type evaluate_all_splines(const double* const State_Vector, int radial_grid_idx, int theta_grid_idx, Derivative_selector_enums Derivative_selector) const;
    Numerical_metric_potentials_type compute_metric_components_from_spline(const double* const State_vector, int radial_grid_idx, int theta_grid_idx, Derivative_selector_enums Derivative_selector) const;

    double compactify_radial_coordiante(const double r) const;
    double uncompactify_radial_coordinate(const double x);

    Metric_type get_metric(const double* const State_Vector, int radial_grid_idx, int theta_grid_idx) const;
    Metric_type get_dr_metric(const double* const State_Vector, int radial_grid_idx, int theta_grid_idx) const;
    Metric_type get_dtheta_metric(const double* const State_Vector, int radial_grid_idx, int theta_grid_idx) const;
    Metric_type get_d2r_metric(const double* const State_Vector, int radial_grid_idx, int theta_grid_idx) const;
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