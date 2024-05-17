#pragma once

#define _USE_MATH_DEFINES

#include "General_math_functions.h"
#include "General_GR_functions.h"
#include "Disk_Models.h"
#include "Spacetimes.h"
#include "Constants.h"

#include <vector>

/***************************************************
|                                                  |
| Novikov-Thorne Model Class Functions Definitions |
|                                                  |
***************************************************/

Novikov_Thorne_Model::Novikov_Thorne_Model(double x, double y, Spacetime_Base_Class* Spacetimes[]) {

    r_in = x;

    if (x == NULL) {

        r_in = Spacetimes[e_metric]->get_ISCO()[Outer];

    }

    r_out = y;

};

double Novikov_Thorne_Model::Keplerian_angular_velocity(double r, Spacetime_Base_Class* Spacetimes[]) {

    double State_Vector[2] = { r, M_PI_2 };

    Metric_type s_dr_Metric = Spacetimes[e_metric]->get_dr_metric(State_Vector);

    return (-s_dr_Metric.Metric[0][3] + sqrt(s_dr_Metric.Metric[0][3] * s_dr_Metric.Metric[0][3] - s_dr_Metric.Metric[0][0] * s_dr_Metric.Metric[3][3])) / s_dr_Metric.Metric[3][3];

}

double Novikov_Thorne_Model::dr_Keplerian_angular_velocity(double r, Spacetime_Base_Class* Spacetimes[]) {

    double State_Vector[2] = { r, M_PI_2 };

    Metric_type s_dr_Metric = Spacetimes[e_metric]->get_dr_metric(State_Vector);
    Metric_type s_d2r_Metric = Spacetimes[e_metric]->get_d2r_metric(State_Vector);

    double root = sqrt(s_dr_Metric.Metric[0][3] * s_dr_Metric.Metric[0][3] - s_dr_Metric.Metric[0][0] * s_dr_Metric.Metric[3][3]);

    double Kepler = this->Keplerian_angular_velocity(r, Spacetimes);

    return  -Kepler / s_dr_Metric.Metric[3][3] * s_d2r_Metric.Metric[3][3] + (-s_d2r_Metric.Metric[0][3]
        + 1.0 / root / 2 * (2 * s_dr_Metric.Metric[0][3] * s_d2r_Metric.Metric[0][3] - s_dr_Metric.Metric[0][0] * s_d2r_Metric.Metric[3][3]
            - s_d2r_Metric.Metric[0][0] * s_dr_Metric.Metric[3][3])) / s_dr_Metric.Metric[3][3];

}

double Novikov_Thorne_Model::Redshift(double J, double State_Vector[], double r_obs, double theta_obs,
    Spacetime_Base_Class* Spacetimes[]) {

    double& r_source = State_Vector[e_r];
    double& theta_source = State_Vector[e_theta];

    /*
    Get the observer 4-velocity
    */

    double State_Vector_obs[2] = { r_obs, theta_obs };

    Metric_type s_Metric_obs = Spacetimes[e_metric]->get_metric(State_Vector_obs);

    double U_obs[4] = { 1.0 / s_Metric_obs.Lapse_function, 0 ,0 , s_Metric_obs.Shift_function / s_Metric_obs.Lapse_function };

    /*
    Get the source 4-velocity
    */

    Metric_type s_Metric_source = Spacetimes[e_metric]->get_metric(State_Vector);

    double Kepler = this->Keplerian_angular_velocity(r_source, Spacetimes);

    double gamma = 1 / sqrt(-s_Metric_source.Metric[0][0] - 2 * s_Metric_source.Metric[0][3] * Kepler - s_Metric_source.Metric[3][3] * Kepler * Kepler);

    if (isnan(gamma) || isinf(gamma)) {

        exit(ERROR);
    }

    double U_source[4] = { gamma, 0, 0, gamma * Kepler };

    return  (-U_obs[0] + U_obs[3] * J) / (-U_source[0] + U_source[3] * J);

}

double Novikov_Thorne_Model::disk_Energy(double r, Spacetime_Base_Class* Spacetimes[]) {

    double State_Vector[2] = { r, M_PI_2 };

    Metric_type s_Metric_source = Spacetimes[e_metric]->get_metric(State_Vector);

    double Kepler = this->Keplerian_angular_velocity(r, Spacetimes);

    double root = sqrt(-s_Metric_source.Metric[0][0] - 2 * s_Metric_source.Metric[0][3] * Kepler - s_Metric_source.Metric[3][3] * Kepler * Kepler);

    return  -(s_Metric_source.Metric[0][0] + s_Metric_source.Metric[0][3] * Kepler) / root;

}

double Novikov_Thorne_Model::disk_Angular_Momentum(double r, Spacetime_Base_Class* Spacetimes[]) {

    double State_Vector[2] = { r, M_PI_2 };

    Metric_type s_Metric_source = Spacetimes[e_metric]->get_metric(State_Vector);

    double Kepler = this->Keplerian_angular_velocity(r, Spacetimes);

    double root = sqrt(-s_Metric_source.Metric[0][0] - 2 * s_Metric_source.Metric[0][3] * Kepler - s_Metric_source.Metric[3][3] * Kepler * Kepler);

    return  (s_Metric_source.Metric[3][3] * Kepler + s_Metric_source.Metric[0][3]) / root;

}

double Novikov_Thorne_Model::Flux_integrand(double r, Spacetime_Base_Class* Spacetimes[]) {

    double State_Vector[2] = { r, M_PI_2 };

    Metric_type s_Metric = Spacetimes[e_metric]->get_metric(State_Vector);
    Metric_type s_dr_Metric = Spacetimes[e_metric]->get_dr_metric(State_Vector);

    double Kepler = this->Keplerian_angular_velocity(r, Spacetimes);
    double dr_Kepler = this->dr_Keplerian_angular_velocity(r, Spacetimes);

    double root = sqrt(-s_Metric.Metric[0][0] - 2 * s_Metric.Metric[0][3] * Kepler - s_Metric.Metric[3][3] * Kepler * Kepler);
    double dr_root = (-s_dr_Metric.Metric[0][0] - 2 * (s_dr_Metric.Metric[0][3] * Kepler + s_Metric.Metric[0][3] * dr_Kepler)
        - s_dr_Metric.Metric[3][3] * Kepler * Kepler - 2 * s_Metric.Metric[3][3] * Kepler * dr_Kepler);

    double E = this->disk_Energy(r, Spacetimes);
    double L = this->disk_Angular_Momentum(r, Spacetimes);

    double dr_L = (s_dr_Metric.Metric[3][3] * Kepler + s_Metric.Metric[3][3] * dr_Kepler + s_dr_Metric.Metric[0][3]) / root - L / root / root / 2 * dr_root;

    return (E - Kepler * L) * dr_L;

}

double Novikov_Thorne_Model::solve_Flux_integral(double lower_bound, double upper_bound, double tolerance, Spacetime_Base_Class* Spacetimes[]) {

    double mid_point = (lower_bound + upper_bound) / 2;
    double left_mid_point = (lower_bound + mid_point) / 2;
    double right_mid_point = (mid_point + upper_bound) / 2;

    double F_lower_bound = this->Flux_integrand(lower_bound, Spacetimes);
    double F_mid_point = this->Flux_integrand(mid_point, Spacetimes);
    double F_upper_bound = this->Flux_integrand(upper_bound, Spacetimes);

    double F_left_mid = this->Flux_integrand(left_mid_point, Spacetimes);
    double F_right_mid = this->Flux_integrand(right_mid_point, Spacetimes);

    double S_left = (mid_point - lower_bound) / 6 * (F_lower_bound + 4 * F_left_mid + F_mid_point);
    double S_right = (upper_bound - mid_point) / 6 * (F_mid_point + 4 * F_right_mid + F_upper_bound);

    double S_2 = S_left + S_right;
    double S_1 = (upper_bound - lower_bound) / 6 * (F_lower_bound + 4 * F_mid_point + F_upper_bound);

    double error;

    if (S_2 >= S_1) {

        error = S_2 - S_1;

    }
    else {

        error = S_1 - S_2;

    }

    double integral;

    if (error < 15 * tolerance) {

        integral = S_2 + (S_2 - S_1) / 15;


    }
    else {

        double L_value = this->solve_Flux_integral(lower_bound, mid_point, tolerance / 2, Spacetimes);
        double R_value = this->solve_Flux_integral(mid_point, upper_bound, tolerance / 2, Spacetimes);

        integral = L_value + R_value;

    }

    return integral;
}

double Novikov_Thorne_Model::get_flux(double r, Spacetime_Base_Class* Spacetimes[]) {

    double State_Vector[2] = { r, M_PI_2 };

    if (e_metric == Wormhole) {

        r = sqrt(r * r + WH_R_THROAT * WH_R_THROAT);

    }

    Metric_type s_Metric = Spacetimes[e_metric]->get_metric(State_Vector);

    double metric_det = get_metric_det(s_Metric.Metric);
    double E_disk = disk_Energy(r, Spacetimes);
    double L_disk = disk_Angular_Momentum(r, Spacetimes);

    double Kepler = Keplerian_angular_velocity(r, Spacetimes);
    double dr_Kepler = dr_Keplerian_angular_velocity(r, Spacetimes);

    double Flux_coeff = -dr_Kepler / ((E_disk - Kepler * L_disk) * (E_disk - Kepler * L_disk)) / (4 * M_PI * sqrt(-metric_det));

    double Flux_integral = solve_Flux_integral(r_in, r, INTEGRAL_ACCURACY, Spacetimes);


    return Flux_coeff * Flux_integral;

}

double Novikov_Thorne_Model::get_r_in() { return r_in; };
double Novikov_Thorne_Model::get_r_out() { return r_out; };

/************************************************************
|                                                           |
| Optically Thin Toroidal Model Class Functions Definitions |
|                                                           |
************************************************************/

Optically_Thin_Toroidal_Model::Optically_Thin_Toroidal_Model(Disk_model_parameters* p_Disk_Parameters, Emission_law_parameters* p_Emission_Parameters) {

    if (NULL != p_Disk_Parameters) {

        this->s_Disk_Parameters = *p_Disk_Parameters;

    }

    if (NULL != p_Emission_Parameters) {

        this->s_Emission_Parameters = *p_Emission_Parameters;

    }

}

int Optically_Thin_Toroidal_Model::load_parameters(Disk_model_parameters* p_Disk_Parameters, Emission_law_parameters* p_Emission_Parameters) {

    if (NULL != p_Disk_Parameters) {

        this->s_Disk_Parameters = *p_Disk_Parameters;

    }
    else {

        std::cout << "p_Disk_Parameters is a NULL pointer! \n";

        exit(ERROR);
    }

    if (NULL != p_Emission_Parameters) {

        this->s_Emission_Parameters = *p_Emission_Parameters;

    }
    else {

        std::cout << "p_Emission_Parameters is a NULL pointer! \n";

        exit(ERROR);
    }

    return OK;

}

void Optically_Thin_Toroidal_Model::update_hotspot_position(int simulation_frame) {

    this->s_Disk_Parameters.Hotspot_position[e_phi] += 2 * M_PI / HOTSPOT_ANIMATION_NUMBER;

}

Disk_model_parameters Optically_Thin_Toroidal_Model::get_disk_params() { return this->s_Disk_Parameters; };

/* Disk Temperature Functions */

int Optically_Thin_Toroidal_Model::update_disk_temperature(double State_vector[]) {

    double& r = State_vector[e_r];
    double Radial_Cutoff{};

    double T = T_ELECTRON_EXACT_CGS * this->s_Disk_Parameters.Power_law_radial_scale / r;

    if (r < this->s_Disk_Parameters.Disk_r_cutoff) {

        Radial_Cutoff = (r - this->s_Disk_Parameters.Disk_r_cutoff) / this->s_Disk_Parameters.Disk_cutoff_scale;

        T *= exp(-Radial_Cutoff * Radial_Cutoff);

    }

    this->Disk_Temperature = T;

    if (T < 0 || isnan(T) || isinf(T)) {

        return ERROR;

    }

    return OK;

}

double Optically_Thin_Toroidal_Model::get_disk_temperature(double State_vector[]) {

    int result = this->update_disk_temperature(State_vector);

    if (ERROR == result) {

        std::cout << "Invalid Disk Temperature: " << this->Disk_Temperature << "\n";

        exit(ERROR);

    }

    return this->Disk_Temperature;

}

/* Disk Velocity Functions */

int Optically_Thin_Toroidal_Model::update_disk_velocity(double State_Vector[], Initial_conditions_type* s_Initial_Conditions) {

    double& r_source     = State_Vector[e_r];
    double& theta_source = State_Vector[e_theta];

    Metric_type s_Metric = s_Initial_Conditions->Spacetimes[e_metric]->get_metric(State_Vector);

    double rho = r_source * sin(theta_source);

    double sqrt_rho = sqrt(rho);
    double ell = sqrt_rho * sqrt_rho * sqrt_rho / (1 + rho);

    if (e_metric == Naked_Singularity) {

        double gamma = s_Initial_Conditions->Spacetimes[e_metric]->get_parameters().JNW_Gamma_Parameter;
        double r_singularity = 2 / gamma;

        ell *= pow(1 - r_singularity / r_source, gamma);

    }
    else if (e_metric == Wormhole) {


        ell *= (1 - WH_R_THROAT / r_source);

    }

    double u_t{}, u_phi{};

    double inv_metric[4][4]{};

    invert_metric(inv_metric, s_Metric.Metric);

    u_t = -1.0 / sqrt(-(inv_metric[0][0] - 2 * inv_metric[0][3] * ell + inv_metric[3][3] * ell * ell));
    u_phi = -u_t * ell;

    /*
    
    Convert U_source to contravariant components

    */

    this->Disk_velocity[e_t_coord]     = inv_metric[0][0] * u_t + inv_metric[0][3] * u_phi;
    this->Disk_velocity[e_r_coord]     = 0.0;
    this->Disk_velocity[e_theta_coord] = 0.0;
    this->Disk_velocity[e_phi_coord]   = inv_metric[3][3] * u_phi + inv_metric[3][0] * u_t;

    if (isnan(this->Disk_velocity[e_t_coord])   || 
        isinf(this->Disk_velocity[e_t_coord])   || 
        isnan(this->Disk_velocity[e_phi_coord]) ||
        isinf(this->Disk_velocity[e_phi_coord])) {

        return ERROR;

    }

    return OK;

}

double* Optically_Thin_Toroidal_Model::get_disk_velocity(double State_vector[], Initial_conditions_type* s_Initial_Conditions) {

    int result = this->update_disk_velocity(State_vector, s_Initial_Conditions);

    if (ERROR == result) {

        std::cout << "Invalid disk 4-velocity: "
            << "["
            << this->Disk_velocity[e_t_coord]
            << ", "
            << this->Disk_velocity[e_r_coord]
            << ", "
            << this->Disk_velocity[e_theta_coord]
            << ", "
            << this->Disk_velocity[e_phi_coord]
            << "]\n";

        exit(ERROR);

    }

    return this->Disk_velocity;

}

/* Disk Density Functions */

int Optically_Thin_Toroidal_Model::update_disk_hotspot(double State_Vector[]) {

    double& Hotspot_r     = this->s_Disk_Parameters.Hotspot_position[e_r];
    double& Hotspot_theta = this->s_Disk_Parameters.Hotspot_position[e_theta];
    double& Hotspot_phi   = this->s_Disk_Parameters.Hotspot_position[e_phi];

    double x_center = Hotspot_r * sin(Hotspot_theta) * cos(Hotspot_phi);
    double y_center = Hotspot_r * sin(Hotspot_theta) * sin(Hotspot_phi);
    double z_center = Hotspot_r * cos(Hotspot_theta);

    double& r = State_Vector[e_r];

    double x_photon = r * sin(State_Vector[e_theta]) * cos(State_Vector[e_phi]);
    double y_photon = r * sin(State_Vector[e_theta]) * sin(State_Vector[e_phi]);
    double z_photon = r * cos(State_Vector[e_theta]);

    double hotspot_density = exp(-(x_center - x_photon) * (x_center - x_photon) / HOTSPOT_SCALE / HOTSPOT_SCALE);
    hotspot_density *= exp(-(y_center - y_photon) * (y_center - y_photon) / HOTSPOT_SCALE / HOTSPOT_SCALE);
    hotspot_density *= exp(-(z_center - z_photon) * (z_center - z_photon) / HOTSPOT_SCALE / HOTSPOT_SCALE);
    hotspot_density *= HOTSPOT_REL_SCALE;

    this->Disk_hotspot_density = hotspot_density;

    if (isnan(Disk_hotspot_density) || 
        isinf(Disk_hotspot_density) || 
        Disk_hotspot_density < 0) {

        return ERROR;

    }

    return OK;

}

double Optically_Thin_Toroidal_Model::get_disk_hotspot(double State_Vector[]) {

    int result = this->update_disk_hotspot(State_Vector);

    if (ERROR != result) {

        return this->Disk_hotspot_density;

    }
    else {

        std::cout << "Invalid Hotspot density: " << this->Disk_hotspot_density << "\n";

        exit(ERROR);

    }

}

int Optically_Thin_Toroidal_Model::update_disk_density_profile(double State_Vector[]) {

    double& r = State_Vector[e_r];
    double rho = sin(State_Vector[e_theta]);
    double h = cos(State_Vector[e_theta]);

    double Height_Cutoff{};
    double Radial_Cutoff{};

    double electron_density_profile{};

    switch (e_disk_model) {

    case Power_law:

        Height_Cutoff = h / (this->s_Disk_Parameters.Disk_opening_angle * rho);
        electron_density_profile = exp(-Height_Cutoff * Height_Cutoff / 2) / (r / R_0) / (r / R_0);

        if (r < this->s_Disk_Parameters.Disk_r_cutoff) {

            Radial_Cutoff = (r - this->s_Disk_Parameters.Disk_r_cutoff) / this->s_Disk_Parameters.Disk_cutoff_scale;
            electron_density_profile *= exp(-Radial_Cutoff * Radial_Cutoff);

        }

        break;

    case Exponential_law:

        Height_Cutoff = h / this->s_Disk_Parameters.Exp_law_height_scale;
        Radial_Cutoff = r / this->s_Disk_Parameters.Exp_law_radial_scale;

        electron_density_profile = exp(-Radial_Cutoff * Radial_Cutoff / 2 - Height_Cutoff * Height_Cutoff / 2);

        break;

    default:

        std::cout << "Wrong disk density model!" << '\n';

    }

    if (this->s_Disk_Parameters.Hotspot_scale != 0.0f) {

        this->Disk_density_profile = electron_density_profile + get_disk_hotspot(State_Vector);

    }
    else {

        this->Disk_density_profile = electron_density_profile;

    }

    if (this->Disk_density_profile < 0) {

        return ERROR;

    }

    return OK;

}

double Optically_Thin_Toroidal_Model::get_disk_density_profile(double State_Vector[]) {

    int result = this->update_disk_density_profile(State_Vector);

    if (ERROR == result) {

        exit(ERROR);

    }

    return this->Disk_density_profile;

}

/* Disk Magnetic Field Functions */

double Optically_Thin_Toroidal_Model::get_magnetic_field(double B_field[4], 
                                                         double State_vector[], 
                                                         Initial_conditions_type* s_Initial_Conditions) {

    /*
    
    Computes the magnetic field 4-vector, measured by a comoving obverver (with 4-velocity Plasma_velocity), called B_field,
    in the basis of a static observer (with 4-velocity n^mu = {1, 0, 0, 0} ). 
    Returns the magnitude of the magnetic field in the plasma frame.

    NOTE: The magnitude of the magnetic field in these frames is different, because its not concerved under Lorentz boosts.
    In the static frame I only set the geometry of the field, then finally scale B_field (the plasma frame one) by B_Plasma_norm_CGS.
    
    NOTE: B_field is measured by a static observer, so:
    1) B_field is effectively in the coordinate basis (not any ZAMO!) - this means one cannot directly apply Lorentz boosts to it.

    */

    double electron_density = this->s_Disk_Parameters.Density_scale * this->get_disk_density_profile(State_vector);

    double B_Plasma_norm_CGS = sqrt(this->s_Disk_Parameters.Magnetization * C_LIGHT_CGS * C_LIGHT_CGS * electron_density * M_PROTON_CGS * 4 * M_PI);
    double B_Static[4] = {         0.0, 
                           MAG_FIELD_GEOMETRY[0],
                           MAG_FIELD_GEOMETRY[1],
                           MAG_FIELD_GEOMETRY[2] };

    double* Plasma_velocity = this->get_disk_velocity(State_vector, s_Initial_Conditions);

    /* 

    The only reference I could find for this is https://iopscience.iop.org/article/10.3847/1538-4357/ab718e/pdf 8e) and 8f).
    One arrives at the expression by projecting F^mu^nu (written in terms of B_Plasma and Plasma_veclocity),
    onto the 4-velocity of the static observer - this gives B_Static, then inverting the expression 
    to obtain B_Plasma = f(B_Static, Plasma_Velocity)
    
    */

    Metric_type s_Metric = s_Initial_Conditions->Spacetimes[e_metric]->get_metric(State_vector);

    for (int left_idx = 0; left_idx <= 3; left_idx++) {

        for (int right_idx = 0; right_idx <= 3; right_idx++) {

            B_field[e_t_coord] += s_Metric.Metric[left_idx][right_idx] * Plasma_velocity[left_idx] * B_Static[right_idx];

        }

    }

    for (int index = 1; index <= 3; index++) {

        B_field[index] = (B_Static[index] + B_field[e_t_coord] * Plasma_velocity[index]) / Plasma_velocity[e_t_coord];
            
    }
    
    return B_Plasma_norm_CGS;

}

double Optically_Thin_Toroidal_Model::get_electron_pitch_angle(double State_Vector[], double B_field_local[4], Initial_conditions_type* s_Initial_Conditions) {

    /*
    
    Computes the angle between the photon wave-3-vector and the magnetic field, measured by a comoving (with the plasma) observer.
    There is a neat invariant way to compute this directly from the coordinate basis 4-vectors.
    
    */

    double* Plasma_velocity = this->get_disk_velocity(State_Vector, s_Initial_Conditions);
    double Wave_vec_dot_Plasma_vec =            -1           * Plasma_velocity[e_t_coord] + 
                                     State_Vector[e_p_r]     * Plasma_velocity[e_r_coord] +
                                     State_Vector[e_p_theta] * Plasma_velocity[e_theta_coord] +
                                     State_Vector[e_p_phi]   * Plasma_velocity[e_phi_coord];

    Metric_type s_Metric  = s_Initial_Conditions->Spacetimes[e_metric]->get_metric(State_Vector);
    double B_field_norm_squared{};
    double B_field_dot_Plasma_vel{};

    /*
    
    TODO: Maybe make functions that do this, or functions that raise and lower indicies
    
    */

    for (int left_idx = 0; left_idx <= 3; left_idx++) {

        for (int right_idx = 0; right_idx <= 3; right_idx++) {

            B_field_norm_squared   += s_Metric.Metric[left_idx][right_idx] * B_field_local[left_idx] * B_field_local[right_idx];
            B_field_dot_Plasma_vel += s_Metric.Metric[left_idx][right_idx] * B_field_local[left_idx] * Plasma_velocity[right_idx];

        }

    }

    double Wave_vec_dot_B_field =           -1           * B_field_local[e_t_coord] +
                                 State_Vector[e_p_r]     * B_field_local[e_r_coord] +
                                 State_Vector[e_p_theta] * B_field_local[e_theta_coord] +
                                 State_Vector[e_p_phi]   * B_field_local[e_phi_coord];

    double cos_angle = 1;

    if (!isinf(1.0 / Wave_vec_dot_Plasma_vec) && !isinf(1.0 / B_field_norm_squared)) {

        cos_angle = (Wave_vec_dot_B_field + B_field_dot_Plasma_vel * Wave_vec_dot_Plasma_vec) / Wave_vec_dot_Plasma_vec / sqrt(B_field_norm_squared + B_field_dot_Plasma_vel * B_field_dot_Plasma_vel);

    }

    if (fabs(cos_angle) < 1.0) {

        return cos_angle;

    }
    else {

        return cos_angle / fabs(cos_angle);

    }

}

/* =============================================== Disk Transfer Functions Functions =============================================== */

void Optically_Thin_Toroidal_Model::get_faradey_fit_functions(double X,
                                                              double X_to_1_point_2,
                                                              double X_frac,
                                                              double faradey_fucntions[STOKES_PARAM_NUM]) {

    /* The reference for this implementation is from Appendix B2 of https://arxiv.org/pdf/1602.03184.pdf */

    faradey_fucntions[I] = 0.0f; // This is zero by definition
    faradey_fucntions[U] = 0.0f; // This is zero by definition
    faradey_fucntions[Q] = 0.0f;
    faradey_fucntions[V] = 0.0f;

    constexpr double TWO_TO_MINUS_ONE_THIRD = 0.7937005259840997373758528196361;
    constexpr double THREE_TO_23_OVER_6     = 67.44733739;

    if (X > std::numeric_limits<double>::min() && X < std::numeric_limits<double>::infinity()) {

        faradey_fucntions[Q] = 2.011 * exp(-X_frac / 4.7) - cos(X / 2) * exp(-X_to_1_point_2 / 2.73) - 0.011 * exp(-X / 47.2) +
            (0.011 * exp(-X / 47.2) - TWO_TO_MINUS_ONE_THIRD / THREE_TO_23_OVER_6 * 1e4 * M_PI * pow(X, -8.0 / 3)) / 2 * (1 + tanh(10 * log(X / 120)));
        faradey_fucntions[V] = 1.0 - 0.11 * log(1.0 + 0.035 * X);

    }

}

void Optically_Thin_Toroidal_Model::get_synchotron_emission_fit_function(Sync_emission_fit_functions e_fit_functions,
    double Emission_functions[4],
    double X,
    double sqrt_X,
    double cbrt_X) {

    Emission_functions[I] = 0;
    Emission_functions[U] = 0;
    Emission_functions[Q] = 0;
    Emission_functions[V] = 0;

    switch (e_fit_functions) {

    case Leung_2011:

        /*

        The below expressions were originally designed to work with f_s as defined in the paper = 4 / 27 * f_crit. This is annoying for my implementation (I define X = f / f_crit),
        so I convert the expressions to work with f_crit explicitly.

        Ontop of that the whole function scales linearly with frequency. This is replaced by the variable "X" in the second line (in order to take the redshift into account without passing it in as an argument, because I don't like that).
        This is then corrected in the caller function, by multiplying with the critical frequency.

        Ref: https://www.aanda.org/articles/aa/pdf/2022/11/aa44339-22.pdf

        */

        Emission_functions[I] = 0.0;

        if (!isnan(cbrt_X) && !isinf(cbrt_X) && !isnan(X) && !isinf(X)) {

            Emission_functions[I] = M_SQRT2 * M_PI * Q_ELECTRON_CGS * Q_ELECTRON_CGS / 3 / C_LIGHT_CGS;
            Emission_functions[I] *= (4. / 27 * X);
            Emission_functions[I] *= ((2.5980762113) + (1.3747296369986) * 1.887749 / cbrt_X); // The constants inside parentasies here are sqrt(27. / 4) and sqrt(cbrt(27. / 4)), respectfully
            Emission_functions[I] *= ((2.5980762113) + (1.3747296369986) * 1.887749 / cbrt_X); // The other constant is pow (2, 11. / 12)
            Emission_functions[I] *= exp(-1.889881574 * cbrt_X);                               // This constant is cbrt(27. / 4)
        }

        break;

    case Dexter_2016:

        /*

        Below for the I and Q functions, the frequency is replaced by the variable "X" (in order to take the redshift into account without passing it in as an argument, because I don't like that).
        This is then corrected in the caller function, by multiplying with the critical frequency.

        Ref: https://arxiv.org/pdf/1602.03184.pdf

        */

        if (!isnan(cbrt_X) && !isinf(cbrt_X) && !isnan(X) && !isinf(X) && !isnan(sqrt_X) && !isinf(sqrt_X)) {

            double exponenet = exp(-1.8899 * cbrt_X);

            Emission_functions[I] = Q_ELECTRON_CGS * Q_ELECTRON_CGS / 1.73205080757 / C_LIGHT_CGS * X * 2.5651 * (1 + 1.92 / cbrt_X + 0.9977 / cbrt_X / cbrt_X) * exponenet;
            Emission_functions[Q] = Q_ELECTRON_CGS * Q_ELECTRON_CGS / 1.73205080757 / C_LIGHT_CGS * X * 2.5651 * (1 + 0.932 / cbrt_X + 0.4998 / cbrt_X / cbrt_X) * exponenet;
            Emission_functions[V] = 4 * Q_ELECTRON_CGS * Q_ELECTRON_CGS / 1.73205080757 / C_LIGHT_CGS / 3 * X * (1.8138 / X + 3.423 / cbrt_X / cbrt_X + 0.02955 / sqrt_X + 2.0377 / cbrt_X) * exponenet;
        }

        break;

    default:

        break;
    }

}

void Optically_Thin_Toroidal_Model::get_transfer_fit_functions(double Emission_fucntions[STOKES_PARAM_NUM],
                                                               double Faradey_functions[STOKES_PARAM_NUM],
                                                               Emission_functions_arguments* Arguments) {

    if (INCLUDE_POLARIZATION){

        this->get_synchotron_emission_fit_function(Dexter_2016, Emission_fucntions, Arguments->X_emission, Arguments->X_1_2_emission, Arguments->X_1_3_emission);

        this->get_faradey_fit_functions(Arguments->X_faradey, Arguments->X_to_1_point_2_faradey, Arguments->X_to_1_point_035_faradey, Faradey_functions);

    }
    else {

        this->get_synchotron_emission_fit_function(Leung_2011, Emission_fucntions, Arguments->X_emission, Arguments->X_1_2_emission, Arguments->X_1_3_emission);

    }

}

void Optically_Thin_Toroidal_Model::get_synchotron_transfer_functions(double State_vector[],
                                                                      Initial_conditions_type* s_Initial_conditions,
                                                                      double Emission_fucntions[STOKES_PARAM_NUM],
                                                                      double Faradey_functions[STOKES_PARAM_NUM],
                                                                      double Absorbtion_functions[STOKES_PARAM_NUM]) {

    /* The transfer functions arrays needs to be manually cleared, because this function only adds to it. */
    for (int stokes_index = 0; stokes_index <= STOKES_PARAM_NUM - 1; stokes_index++) {

        Emission_fucntions[stokes_index] = 0.0;
        Faradey_functions[stokes_index] = 0.0;
        Absorbtion_functions[stokes_index] = 0.0;
    }

    /* Electron Density in CGS */
    double electron_density = this->s_Disk_Parameters.Density_scale * this->get_disk_density_profile(State_vector);

    /* Dimentionless Electron Temperature */
    double T_electron = this->get_disk_temperature(State_vector);
    double T_electron_dim = BOLTZMANN_CONST_CGS * T_electron / M_ELECTRON_CGS / C_LIGHT_CGS / C_LIGHT_CGS;
    double K0_Bessel = std::cyl_bessel_k(0.0, 1.0 / T_electron_dim);
    double K1_Bessel = std::cyl_bessel_k(1.0, 1.0 / T_electron_dim);
    double K2_Bessel = std::cyl_bessel_k(2.0, 1.0 / T_electron_dim);

    /* Magnetic Field */
    double B_field_local[4]{};
    double B_CGS = this->get_magnetic_field(B_field_local, State_vector, s_Initial_conditions);

    /* Disk Coordinate Velocity */
    double* U_source_coord = this->get_disk_velocity(State_vector, s_Initial_conditions);

    /* Redshit */
    double redshift = Redshift(State_vector, U_source_coord);

    /* Cyclotron Frequency */
    double f_cyclo = Q_ELECTRON_CGS * B_CGS / (2 * M_PI * M_ELECTRON_CGS * C_LIGHT_CGS);

    /* The "averaged" critical frequency (without the sin(theta) term - that gets added on later from a pre-computed table) */
    double f_crit_no_sin = 3. / 2 * f_cyclo * T_electron_dim * T_electron_dim;

    Emission_functions_arguments Arguments_angle_uncorrected{};

    /* Both the emission and faradey function expressions are in terms of an dimentionless variable X, but the definitions for X are different */
    Arguments_angle_uncorrected.X_emission = 1e100;
    Arguments_angle_uncorrected.X_faradey  = 1e100;

    if (f_crit_no_sin > std::numeric_limits<double>::min()) {

        Arguments_angle_uncorrected.X_emission = OBS_FREQUENCY_CGS / f_crit_no_sin / redshift;
        Arguments_angle_uncorrected.X_faradey = T_electron_dim * sqrt(M_SQRT2 * 1e3 * f_cyclo / (OBS_FREQUENCY_CGS / redshift));

    }

    /* All the functions are normalized by a Bessel function, so I check if I can divide by it */
    if (!isinf(1e10 / K2_Bessel)) {

       /* == Compute all the weird powers of X outside the averaging loop == */

        Arguments_angle_uncorrected.X_1_2_emission = sqrt(Arguments_angle_uncorrected.X_emission);
        Arguments_angle_uncorrected.X_1_3_emission = cbrt(Arguments_angle_uncorrected.X_emission);

        Arguments_angle_uncorrected.X_to_1_point_2_faradey = pow(Arguments_angle_uncorrected.X_faradey, 1.2);
        Arguments_angle_uncorrected.X_to_1_point_035_faradey = pow(Arguments_angle_uncorrected.X_faradey, 1.035f);

       /* ================================================================== */

        if (AVERAGE_EMISSION_PITCH_ANGLE) {

            for (int averaging_idx = 0; averaging_idx <= NUM_SAMPLES_TO_AVG - 1; averaging_idx++) {

                Emission_functions_arguments Arguments_angle_corrected{};

                double sin_pitch_angle = this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[averaging_idx];
                double cos_pitch_angle = this->s_Precomputed_e_pitch_angles.cos_electron_pitch_angles[averaging_idx];

                double temp_emission_functions[STOKES_PARAM_NUM]{};
                double temp_faradey_functions[STOKES_PARAM_NUM]{};

                Arguments_angle_corrected.X_emission     = Arguments_angle_uncorrected.X_emission / sin_pitch_angle;
                Arguments_angle_corrected.X_1_2_emission = Arguments_angle_uncorrected.X_1_2_emission * this->s_Precomputed_e_pitch_angles.one_over_sqrt_sin[averaging_idx];
                Arguments_angle_corrected.X_1_3_emission = Arguments_angle_uncorrected.X_1_3_emission * this->s_Precomputed_e_pitch_angles.one_over_cbrt_sin[averaging_idx];

                Arguments_angle_corrected.X_faradey                = Arguments_angle_corrected.X_faradey / this->s_Precomputed_e_pitch_angles.one_over_sqrt_sin[averaging_idx];
                Arguments_angle_corrected.X_to_1_point_035_faradey = Arguments_angle_corrected.X_to_1_point_035_faradey / this->s_Precomputed_e_pitch_angles.one_over_sin_to_1_point_035[averaging_idx];
                Arguments_angle_corrected.X_to_1_point_2_faradey   = Arguments_angle_corrected.X_to_1_point_2_faradey /  this->s_Precomputed_e_pitch_angles.one_over_sin_to_1_point_2_over_2[averaging_idx];

                this->get_transfer_fit_functions(temp_emission_functions, temp_faradey_functions, &Arguments_angle_corrected);

                /* ================================================ The emission functions ================================================ */

                temp_emission_functions[I] *= electron_density * (f_crit_no_sin * sin_pitch_angle) / K2_Bessel; // Scale the emission by the remaining position-dependant factors
                temp_emission_functions[I] *= sin_pitch_angle * M_PI / NUM_SAMPLES_TO_AVG / 2;                  // Scale by the factors coming form averaging over all emission orientations

                temp_emission_functions[Q] *= electron_density * (f_crit_no_sin * sin_pitch_angle) / K2_Bessel; // Scale the emission by the remaining position-dependant factors
                temp_emission_functions[Q] *= sin_pitch_angle * M_PI / NUM_SAMPLES_TO_AVG / 2;                  // Scale by the factors coming form averaging over all emission orientations
                 
                if (!isinf(1e2 / sin_pitch_angle)) {

                    temp_emission_functions[V] *= electron_density * (f_crit_no_sin * sin_pitch_angle) / K2_Bessel; // Scale the emission by the remaining position-dependant factors
                    temp_emission_functions[V] *= (1. / T_electron_dim) * cos_pitch_angle / sin_pitch_angle;        // The V component has some extra angle dependance
                    temp_emission_functions[V] *= sin_pitch_angle * M_PI / NUM_SAMPLES_TO_AVG / 2;                  // Scale by the factors coming form averaging over all emission orientations

                }

                Emission_fucntions[I] += temp_emission_functions[I];
                Emission_fucntions[Q] += temp_emission_functions[Q];
                Emission_fucntions[U]  = 0.0;
                Emission_fucntions[V] += temp_emission_functions[V];

                /* ================================================ The faradey functions ================================================ */
                /* Originally derived in https://iopscience.iop.org/article/10.1086/592326/pdf - expressions 25, 26 and 33 */

                double const omega_plasma_squared = 4 * M_PI * electron_density * Q_ELECTRON_CGS * Q_ELECTRON_CGS / M_ELECTRON_CGS;

                temp_faradey_functions[Q] *= omega_plasma_squared * (2 * M_PI * f_cyclo) * (2 * M_PI * f_cyclo) * sin_pitch_angle * sin_pitch_angle * (K1_Bessel / K2_Bessel + 6 * T_electron_dim);
                temp_faradey_functions[Q] /= 2 * C_LIGHT_CGS * (2 * M_PI * OBS_FREQUENCY_CGS / redshift) * (2 * M_PI * OBS_FREQUENCY_CGS / redshift) * (2 * M_PI * OBS_FREQUENCY_CGS / redshift);
                temp_faradey_functions[Q] *= sin_pitch_angle * M_PI / NUM_SAMPLES_TO_AVG / 2; // Scale by the factors coming form averaging over all emission orientations

                temp_faradey_functions[V] *= omega_plasma_squared * (2 * M_PI * f_cyclo) * cos_pitch_angle * (K0_Bessel) / K2_Bessel;
                temp_faradey_functions[V] /= C_LIGHT_CGS * (2 * M_PI * OBS_FREQUENCY_CGS / redshift) * (2 * M_PI * OBS_FREQUENCY_CGS / redshift);
                temp_faradey_functions[V] *= sin_pitch_angle * M_PI / NUM_SAMPLES_TO_AVG / 2; // Scale by the factors coming form averaging over all emission orientations

                // The I and U components are 0 by definition

                Faradey_functions[I]  = 0.0f;
                Faradey_functions[Q] += temp_faradey_functions[Q];
                Faradey_functions[U]  = 0.0f;
                Faradey_functions[V] += temp_faradey_functions[V];

            }

        }
        else {
            
            double cos_pitch_angle = get_electron_pitch_angle(State_vector, B_field_local, s_Initial_conditions);
            double sin_pitch_angle = sqrt(1.0 - cos_pitch_angle * cos_pitch_angle);

            double one_over_sqrt_sin = 1.0 / sqrt(sin_pitch_angle);
            double one_over_cbrt_sin = 1.0 / cbrt(sin_pitch_angle);

            Emission_functions_arguments Arguments_angle_corrected{};

            Arguments_angle_corrected.X_emission     = Arguments_angle_uncorrected.X_emission / sin_pitch_angle;
            Arguments_angle_corrected.X_1_2_emission = Arguments_angle_uncorrected.X_1_2_emission * one_over_sqrt_sin;
            Arguments_angle_corrected.X_1_3_emission = Arguments_angle_uncorrected.X_1_3_emission * one_over_cbrt_sin;

            Arguments_angle_corrected.X_faradey                = Arguments_angle_uncorrected.X_faradey / one_over_sqrt_sin;
            Arguments_angle_corrected.X_to_1_point_035_faradey = Arguments_angle_uncorrected.X_to_1_point_035_faradey * pow(sin_pitch_angle, 1.035);
            Arguments_angle_corrected.X_to_1_point_2_faradey   = Arguments_angle_uncorrected.X_to_1_point_2_faradey * pow(sin_pitch_angle, 1.2);

            this->get_transfer_fit_functions(Emission_fucntions, Faradey_functions, &Arguments_angle_corrected);

            /* ================================================ The emission functions ================================================ */

            Emission_fucntions[I] *= electron_density * (f_crit_no_sin * sin_pitch_angle) / K2_Bessel; // Scale the emission by the remaining position-dependant factors
            Emission_fucntions[Q] *= electron_density * (f_crit_no_sin * sin_pitch_angle) / K2_Bessel; // Scale the emission by the remaining position-dependant factors

            if (!isinf(1.0 / sin_pitch_angle) && !isnan(1.0 / sin_pitch_angle)) {

                Emission_fucntions[V] *= electron_density * (f_crit_no_sin * sin_pitch_angle) / K2_Bessel; // Scale the emission by the remaining position-dependant factors
                Emission_fucntions[V] *= (1. / T_electron_dim) * cos_pitch_angle / sin_pitch_angle;        // The V component has some extra angle dependance

            }

            /* ================================================ The faradey functions ================================================ */
            /* Originally derived in https://iopscience.iop.org/article/10.1086/592326/pdf - expressions 25, 26 and 33 */

             // The I and U components are 0 by definition

            double const omega_plasma_squared = 4 * M_PI * electron_density * Q_ELECTRON_CGS * Q_ELECTRON_CGS / M_ELECTRON_CGS;

            Faradey_functions[I] = 0.0f;

            Faradey_functions[Q] *= omega_plasma_squared * (2 * M_PI * f_cyclo) * (2 * M_PI * f_cyclo) * sin_pitch_angle * sin_pitch_angle * (K1_Bessel / K2_Bessel + 6 * T_electron_dim);
            Faradey_functions[Q] /= 2 * C_LIGHT_CGS * (2 * M_PI * OBS_FREQUENCY_CGS / redshift) * (2 * M_PI * OBS_FREQUENCY_CGS / redshift) * (2 * M_PI * OBS_FREQUENCY_CGS / redshift);

            Faradey_functions[U] = 0.0f;

            Faradey_functions[V] *= omega_plasma_squared * (2 * M_PI * f_cyclo) * cos_pitch_angle * (K0_Bessel) / K2_Bessel;
            Faradey_functions[V] /= C_LIGHT_CGS * (2 * M_PI * OBS_FREQUENCY_CGS / redshift) * (2 * M_PI * OBS_FREQUENCY_CGS / redshift);


        }

    }

    /* ================================================ The absorbtion functions ================================================ */

    double Planck_function_CGS = get_planck_function_CGS(OBS_FREQUENCY_CGS / redshift, this->get_disk_temperature(State_vector));

    if (Planck_function_CGS > std::numeric_limits<double>::min()) {

        for (int stokes_index = 0; stokes_index <= STOKES_PARAM_NUM - 1; stokes_index++) {

            Absorbtion_functions[stokes_index] = Emission_fucntions[stokes_index] / Planck_function_CGS;

        }

    }

}

void Optically_Thin_Toroidal_Model::get_emission_function_synchotron_phenomenological(double State_vector[], Initial_conditions_type* s_Initial_conditions, double Emission_functions[4]) {

    double& emission_power_law = this->s_Emission_Parameters.Emission_power_law;
    double& emission_scale = this->s_Emission_Parameters.Emission_scale;

    double* U_source_coord = this->get_disk_velocity(State_vector, s_Initial_conditions);
    double redshift = Redshift(State_vector, U_source_coord);

    Emission_functions[I] = emission_scale * get_disk_density_profile(State_vector) * pow(redshift, emission_power_law);
    Emission_functions[Q] = 0.0;
    Emission_functions[U] = 0.0;
    Emission_functions[V] = 0.0;

}

void Optically_Thin_Toroidal_Model::get_radiative_transfer_functions(double State_Vector[e_State_Number],
                                                                     Initial_conditions_type* s_Initial_conditions, 
                                                                     double Emission_functions[4], 
                                                                     double Faradey_functions[4],
                                                                     double Absorbtion_functions[4]) {
    switch (e_emission) {

    case(Synchotron_phenomenological):

        this->get_emission_function_synchotron_phenomenological(State_Vector, 
                                                                s_Initial_conditions, 
                                                                Emission_functions);
        
        this->get_absorbtion_function_phenomenological(State_Vector,
                                                       s_Initial_conditions,
                                                       Absorbtion_functions);

        for (int index = I; index <= STOKES_PARAM_NUM - 1; index++) {

            Faradey_functions[index] = 0.0;

        }

        break;

    default:

        this->get_synchotron_transfer_functions(State_Vector, s_Initial_conditions, Emission_functions, Faradey_functions, Absorbtion_functions);

        break;
    }

}

void Optically_Thin_Toroidal_Model::get_absorbtion_function_phenomenological(double State_vector[], 
                                                                             Initial_conditions_type* s_Initial_conditions,
                                                                             double Absorbtion_functions[STOKES_PARAM_NUM]) {

    double& abs_coeff = this->s_Emission_Parameters.Absorbtion_coeff;
    double& emission_scale = this->s_Emission_Parameters.Emission_scale;
    double& source_f_power_law = this->s_Emission_Parameters.Source_f_power_law;
    double& emission_power_law = this->s_Emission_Parameters.Emission_power_law;

    double* U_source = this->get_disk_velocity(State_vector, s_Initial_conditions);
    double redshift = Redshift(State_vector, U_source);

    Absorbtion_functions[I] = abs_coeff * emission_scale * this->get_disk_density_profile(State_vector) * pow(redshift, source_f_power_law + emission_power_law);
    Absorbtion_functions[Q] = 0;
    Absorbtion_functions[U] = 0;
    Absorbtion_functions[V] = 0;

}

void Optically_Thin_Toroidal_Model::precompute_electron_pitch_angles() {

    for (int index = 0; index <= NUM_SAMPLES_TO_AVG - 1; index++) {

        double pitch_angle = double(index) / NUM_SAMPLES_TO_AVG * M_PI;
        this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index] = sin(pitch_angle);
        this->s_Precomputed_e_pitch_angles.cos_electron_pitch_angles[index] = cos(pitch_angle);

        if (this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index] != 0) {

            this->s_Precomputed_e_pitch_angles.one_over_sqrt_sin[index] = 1. / sqrt(this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index]);
            this->s_Precomputed_e_pitch_angles.one_over_cbrt_sin[index] = 1. / cbrt(this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index]);
            this->s_Precomputed_e_pitch_angles.one_over_sin_to_1_point_035[index] = 1. / pow(this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index], 1.035f);
            this->s_Precomputed_e_pitch_angles.one_over_sin_to_1_point_2_over_2[index] = 1. / pow(this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index], 1.2f / 2);
        }

    }

}

double get_planck_function_CGS(double Frequency, double Temperature) {

    return 2 * PLANCK_CONSTANT_CGS * Frequency * Frequency * Frequency / C_LIGHT_CGS / C_LIGHT_CGS / (exp(PLANCK_CONSTANT_CGS * Frequency / BOLTZMANN_CONST_CGS / Temperature) - 1);

}

