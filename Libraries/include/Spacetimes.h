#pragma once
#define _USE_MATH_DEFINES
#include "Structs.h"
#include <complex>
#include <iostream>
#include <cmath>
#include <format>

class Spacetime_Base_Class {

public:

    virtual double* get_ISCO() { throw std::runtime_error("Using Base Spacetime Class. Something Broke in get_ISCO!"); };

    virtual double* get_Photon_Sphere() { throw std::runtime_error("Using Base Spacetime Class. Something Broke in get_Photon_sphere!"); };

    /* --------------------------------------------------- Metric and its derivatives --------------------------------------------------- */

    virtual Metric_type get_metric(const double* const State_Vector) const  { throw std::runtime_error("Using Base Spacetime Class. Something Broke in get_metric!"); };

    virtual Metric_type get_dr_metric(const double* const State_Vector) const  { throw std::runtime_error("Using Base Spacetime Class. Something Broke in get_dr_metric!"); };

    virtual Metric_type get_dtheta_metric(const double* const State_Vector) const  { throw std::runtime_error("Using Base Spacetime Class. Something Broke in get_dtheta_metric!"); };

    virtual Metric_type get_d2r_metric(const double* const State_Vector) const  { throw std::runtime_error("Using Base Spacetime Class. Something Broke in get_d2r_metric!"); };

    /* ------------------------------------------------------ Equations of motion ------------------------------------------------------ */

    virtual void get_EOM(const double* const State_vector, double* const Derivatives) const { throw std::runtime_error("Using Base Spacetime Class. Something Broke in get_EOM!"); };

    /* ---------------------------------------------- Integration Termination Conditions ----------------------------------------------- */

    virtual bool terminate_integration(const double* const State_vector) { throw std::runtime_error("Using Base Spacetime Class. Something Broke in terminate_integration!"); };

};

class Kerr_class : public Spacetime_Base_Class {

private:

    double Mass = 1.0;
    double Spin_Param;
    double Horizon_radius;
    double Scattering_radius;
    double Min_distance_to_singular_point;

public:

    Kerr_class(const Metric_parameters_type* const Metric_Parameters);

    double* get_ISCO() override;
    double* get_Photon_Sphere() override;

    /* Metric and its derivatives */

    Metric_type get_metric(const double* const State_Vector) const override;
    Metric_type get_dr_metric(const double* const State_Vector) const override;
    Metric_type get_dtheta_metric(const double* const State_Vector) const override;
    Metric_type get_d2r_metric(const double* const State_Vector) const override;

    /* Equations of motion */

    void get_EOM(const double* const State_vector, double* const Derivatives) const override;

    /* Integration Termination Conditions */

    bool terminate_integration(const double* const State_vector) override;
   
};

class Minkowski_class : public Spacetime_Base_Class {

private:

    double Scattering_radius;

public:

    Minkowski_class(const Metric_parameters_type* const p_Metric_Parameters);

    /* Metric and its derivatives */

    Metric_type get_metric(const double* const State_Vector) const override;
    Metric_type get_dr_metric(const double* const State_Vector) const override;
    Metric_type get_dtheta_metric(const double* const State_Vector) const override;
    Metric_type get_d2r_metric(const double* const State_Vector) const override;

    /* Equations of motion */

    void get_EOM(const double* const State_vector, double* const Derivatives) const override;

    /* Integration Termination Conditions */

    bool terminate_integration(const double* const State_vector) override;

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

    Wormhole_class(const Metric_parameters_type* const p_Metric_Parameters);

    double* get_ISCO();
    double* get_Photon_Sphere();

    /* Metric and its derivatives */

    Metric_type get_metric(const double* const State_Vector) const override;
    Metric_type get_dr_metric(const double* const State_Vector) const override;
    Metric_type get_dtheta_metric(const double* const State_Vector) const override;
    Metric_type get_d2r_metric(const double* const State_Vector) const override;

    /* Equations of motion */

    void get_EOM(const double* const State_vector, double* const Derivatives) const override;

    /* Integration Termination Conditions */

    bool terminate_integration(const double* const State_vector) override;

};

class RBH_class : public Spacetime_Base_Class {

private:

    double Mass = 1.0;
    double Parameter;

    double Scattering_radius;
    double Min_distance_to_singular_point;
    double Horizon_radius;

public:

    RBH_class(const Metric_parameters_type* const p_Metric_Parameters);

    double* get_ISCO();
    double* get_Photon_Sphere();

    /* Metric and its derivatives */

    Metric_type get_metric(const double* const State_Vector) const override;
    Metric_type get_dr_metric(const double* const State_Vector) const override;
    Metric_type get_dtheta_metric(const double* const State_Vector) const override;
    Metric_type get_d2r_metric(const double* const State_Vector) const override;

    /* Equations of motion */

    void get_EOM(const double* const State_vector, double* const Derivatives) const override;

    /* Integration Termination Conditions */

    bool terminate_integration(const double* const State_vector) override;

};

class JNW_class : public Spacetime_Base_Class {

private:

    double Mass = 1.0;
    double Gamma;

    double Scattering_radius;
    double Min_distance_to_singular_point;
    double Horizon_radius;

public:

    JNW_class(const Metric_parameters_type* const p_Metric_Parameters);

    double* get_ISCO();
    double* get_Photon_Sphere();

    /* Metric and its derivatives */

    Metric_type get_metric(const double* const State_Vector) const override;
    Metric_type get_dr_metric(const double* const State_Vector) const override;
    Metric_type get_dtheta_metric(const double* const State_Vector) const override;
    Metric_type get_d2r_metric(const double* const State_Vector) const override;

    /* Equations of motion */

    void get_EOM(const double* const State_vector, double* const Derivatives) const override;

    /* Integration Termination Conditions */

    bool terminate_integration(const double* const State_vector) override;

};

class Gauss_Bonnet_class : public Spacetime_Base_Class {

private:

    double Mass = 1.0;
    double Gamma;

    double Scattering_radius;
    double Min_distance_to_singular_point;
    double Horizon_radius;

public:

    Gauss_Bonnet_class(const Metric_parameters_type* const p_Metric_Parameters);

    double* get_ISCO() ;
    double* get_Photon_Sphere() ;

    /* Metric and its derivatives */

    Metric_type get_metric(const double* const State_Vector) const override;
    Metric_type get_dr_metric(const double* const State_Vector) const override;
    Metric_type get_dtheta_metric(const double* const State_Vector) const override;
    Metric_type get_d2r_metric(const double* const State_Vector) const override;

    /* Equations of motion */

    void get_EOM(const double* const State_vector, double* const Derivatives) const override;

    /* Integration Termination Conditions */

    bool terminate_integration(const double* const State_vector) override;

};

class Black_Hole_w_Dark_Matter_Halo_class : public Spacetime_Base_Class {

private:

    double Mass = 1.0;
    double Compactness;
    double Halo_Mass;

    double Scattering_radius;
    double Min_distance_to_singular_point;

public:

    Black_Hole_w_Dark_Matter_Halo_class(const Metric_parameters_type* const p_Metric_Parameters);

    double* get_ISCO();

    /* Metric and its derivatives */

    Metric_type get_metric(const double* const State_Vector) const override;
    Metric_type get_dr_metric(const double* const State_Vector) const override;
    Metric_type get_dtheta_metric(const double* const State_Vector) const override;

    /* Equations of motion */

    void get_EOM(const double* const State_vector, double* const Derivatives) const override;

    /* Integration Termination Conditions */

    bool terminate_integration(const double* const State_vector) override;

};

class Numerical_metric : public Spacetime_Base_Class {

private:
    
    Numerical_metric_params_type Parameters;

    double Scattering_radius;
    double Min_distance_to_singular_point;

    /* ------------------------------ Functions that evaluate the cubic B-spline of the metric potentials ------------------------------ */

    void get_control_point_matrix(const double* const Control_vector, const long long r_idx, const long long theta_idx, double Control_matrix[4][4]) const;
    void get_polynomial_basis_vector(const double natural_parameter, double* const Polynomial_basis_vector) const;
    void get_derivative_polynomial_basis_vector(const double natural_parameter, double* const Polynomial_basis_vector) const;
    void get_second_derivative_polynomial_basis_vector(const double natural_parameter, double* const Polynomial_basis_vector) const;

    double evaluate_single_spline(const double Control_point_matrix[4][4], const double Radial_natural_parameter, const double Theta_natural_parameter, Derivative_selector_enums Derivative_selector) const;

    Numerical_metric_potentials_type evaluate_all_splines(const double radial_natural_param, const double theta_natural_param, long long radial_grid_idx, long long theta_grid_idx, Derivative_selector_enums Derivative_selector) const;
    Numerical_metric_potentials_type compute_metric_components_from_spline(Spline_arguments_type s_Splnie_args, Derivative_selector_enums Derivative_selector) const;

    double compactify_radial_coordiante(const double r) const;

    /* --------------------------------------------------- Metric and its derivatives -------------------------------------------------- */

    Metric_type get_metric(const double* const State_Vector, long long radial_grid_idx, long long theta_grid_idx) const;
    Metric_type get_dr_metric(const double* const State_Vector, long long radial_grid_idx, long long theta_grid_idx) const;
    Metric_type get_dtheta_metric(const double* const State_Vector, long long radial_grid_idx, long long theta_grid_idx) const;
    Metric_type get_d2r_metric(const double* const State_Vector, long long radial_grid_idx, long long theta_grid_idx) const;

public:

    Numerical_metric(const Metric_parameters_type* const p_Metric_Parameters);

    /* --------------------------------------------- Metric and its derivatives (wrappers) --------------------------------------------- */

    Metric_type get_metric(const double* const State_Vector) const override;
    Metric_type get_dr_metric(const double* const State_Vector) const override;
    Metric_type get_dtheta_metric(const double* const State_Vector) const override;
    Metric_type get_d2r_metric(const double* const State_Vector) const override;

    /* ------------------------------------------------------ Equations of motion ------------------------------------------------------ */

    void get_EOM(const double* const State_vector, double* const Derivatives) const override;

    /* ---------------------------------------------- Integration Termination Conditions ----------------------------------------------- */

    bool terminate_integration(const double* const State_vector) override;

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