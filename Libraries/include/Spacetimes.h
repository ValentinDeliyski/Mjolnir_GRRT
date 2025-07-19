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

    virtual void get_EOM(const double* const State_vector, double* const Derivatives) const {

        std::cout << "Using Base Spacetime Class - Something Broke!'\n'";
    
    };

    /* Integration Termination Conditions */

    virtual bool terminate_integration(const double* const State_vector) {

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
    double Horizon_radius;
    double Scattering_radius;
    double Min_distance_to_singular_point;

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

    void get_EOM(const double* const State_vector, double* const Derivatives) const override;

    /* Integration Termination Conditions */

    bool terminate_integration(const double* const State_vector) override;

    Return_Values load_parameters(const Metric_parameters_type* const Metric_Parameters) override;
     
};

class Minkowski_class : public Spacetime_Base_Class {

private:

    double Scattering_radius;

public:

    /* Metric and its derivatives */

    Metric_type get_metric(const double* const State_Vector) const override;
    Metric_type get_dr_metric(const double* const State_Vector) const override;
    Metric_type get_dtheta_metric(const double* const State_Vector) const override;
    Metric_type get_d2r_metric(const double* const State_Vector) const override;

    /* Initial conditions derived from images */

    int get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon) override;

    /* Equations of motion */

    void get_EOM(const double* const State_vector, double* const Derivatives) const override;

    /* Integration Termination Conditions */

    bool terminate_integration(const double* const State_vector) override;

    Return_Values load_parameters(const Metric_parameters_type* const Metric_Parameters) override;

};

class Wormhole_class : public Spacetime_Base_Class {

private:

    double Mass = 1.0;
    double R_Throat = this->Mass;
    double Spin_Param;
    double Redshift_Param;

    bool Stop_at_Throat;

    double Scattering_radius;
    double Min_distance_to_throat;

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

    void get_EOM(const double* const State_vector, double* const Derivatives) const override;

    /* Integration Termination Conditions */

    bool terminate_integration(const double* const State_vector) override;

    Return_Values load_parameters(const Metric_parameters_type* const Metric_Parameters) override;

};

class RBH_class : public Spacetime_Base_Class {

private:

    double Mass = 1.0;
    double Parameter;

    double Scattering_radius;
    double Min_distance_to_singular_point;
    double Horizon_radius;

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

    void get_EOM(const double* const State_vector, double* const Derivatives) const override;

    /* Integration Termination Conditions */

    bool terminate_integration(const double* const State_vector) override;

    Return_Values load_parameters(const Metric_parameters_type* const Metric_Parameters) override;

};

class JNW_class : public Spacetime_Base_Class {

private:

    double Mass = 1.0;
    double Gamma;

    double Scattering_radius;
    double Min_distance_to_singular_point;
    double Horizon_radius;

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

    void get_EOM(const double* const State_vector, double* const Derivatives) const override;

    /* Integration Termination Conditions */

    bool terminate_integration(const double* const State_vector) override;

    Return_Values load_parameters(const Metric_parameters_type* const Metric_Parameters) override;

};

class Gauss_Bonnet_class : public Spacetime_Base_Class {

private:;

    double Mass = 1.0;
    double Gamma;

    double Scattering_radius;
    double Min_distance_to_singular_point;
    double Horizon_radius;

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

    void get_EOM(const double* const State_vector, double* const Derivatives) const override;

    /* Integration Termination Conditions */

    bool terminate_integration(const double* const State_vector) override;

    Return_Values load_parameters(const Metric_parameters_type* const Metric_Parameters) override;

};

class Black_Hole_w_Dark_Matter_Halo_class : public Spacetime_Base_Class {

private:

    double Mass = 1.0;
    double Compactness;
    double Halo_Mass;

    double Scattering_radius;
    double Min_distance_to_singular_point;

public:

    double* get_ISCO();

    /* Metric and its derivatives */

    Metric_type get_metric(const double* const State_Vector) const override;
    Metric_type get_dr_metric(const double* const State_Vector) const override;
    Metric_type get_dtheta_metric(const double* const State_Vector) const override;

    /* Initial conditions derived from images */

    int get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon) override;

    /* Equations of motion */

    void get_EOM(const double* const State_vector, double* const Derivatives) const override;

    /* Integration Termination Conditions */

    bool terminate_integration(const double* const State_vector) override;

    Return_Values load_parameters(const Metric_parameters_type* const Metric_Parameters) override;

};

class Numerical_metric : public Spacetime_Base_Class {

private:

    Numerical_metric_params_type Parameters;

    double Scattering_radius;
    double Min_distance_to_singular_point;

    inline void get_control_point_matrix(const double* const Control_vector, const int r_idx, const int theta_idx, double Control_matrix[4][4]) const;

    void get_polynomial_basis_vector(const double natural_parameter, double* const Polynomial_basis_vector) const;
    void get_derivative_polynomial_basis_vector(const double natural_parameter, double* const Polynomial_basis_vector) const;
    void get_second_derivative_polynomial_basis_vector(const double natural_parameter, double* const Polynomial_basis_vector) const;

    double evaluate_single_spline(const double Control_point_matrix[4][4], const double Radial_natural_parameter, const double Theta_natural_parameter, Derivative_selector_enums Derivative_selector) const;

    Numerical_metric_potentials_type evaluate_all_splines(const double* const State_Vector, int radial_grid_idx, int theta_grid_idx, Derivative_selector_enums Derivative_selector) const;
    Numerical_metric_potentials_type compute_metric_components_from_spline(const double* const State_vector, int radial_grid_idx, int theta_grid_idx, Derivative_selector_enums Derivative_selector) const;

    double compactify_radial_coordiante(const double r) const;

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

    void get_EOM(const double* const State_vector, double* const Derivatives) const override;

    /* Integration Termination Conditions */

    bool terminate_integration(const double* const State_vector) override;

    Return_Values load_parameters(const Metric_parameters_type* const Metric_Parameters) override;

};

class Observer_class {

private:

    Observer_parameters_type obs_params;
    double obs_velocity[4];
    double fiducial_obs_velocity[4];

public:

    Observer_class(Simulation_Context_type* p_Sim_Context);

    Observer_parameters_type get_parameters() const;

    const double* const get_obs_velocity() const;

    double* const get_fiducial_obs_velocity(Metric_type* p_Metric);

};