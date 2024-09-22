#pragma once

#define _USE_MATH_DEFINES

#include "General_math_functions.h"
#include "General_GR_functions.h"
#include "Disk_Models.h"
#include "Spacetimes.h"
#include "Constants.h"

#include "gsl/gsl_sf_hyperg.h"

/***************************************************
|                                                  |
| Novikov-Thorne Model Class Functions Definitions |
|                                                  | 
***************************************************/

Novikov_Thorne_Model::Novikov_Thorne_Model(double x, double y, Spacetime_Base_Class* Spacetime) {

    r_in = x;

    if (x == NULL) {

        r_in = Spacetime->get_ISCO()[Outer];

    }

    r_out = y;

};

double Novikov_Thorne_Model::Keplerian_angular_velocity(double r, Spacetime_Base_Class* Spacetimes) {

    double State_Vector[2] = { r, M_PI_2 };

    Metric_type s_dr_Metric = Spacetimes->get_dr_metric(State_Vector);

    return (-s_dr_Metric.Metric[0][3] + sqrt(s_dr_Metric.Metric[0][3] * s_dr_Metric.Metric[0][3] - s_dr_Metric.Metric[0][0] * s_dr_Metric.Metric[3][3])) / s_dr_Metric.Metric[3][3];

}

double Novikov_Thorne_Model::dr_Keplerian_angular_velocity(double r, Spacetime_Base_Class* Spacetime) {

    double State_Vector[2] = { r, M_PI_2 };

    Metric_type s_dr_Metric = Spacetime->get_dr_metric(State_Vector);
    Metric_type s_d2r_Metric = Spacetime->get_d2r_metric(State_Vector);

    double root = sqrt(s_dr_Metric.Metric[0][3] * s_dr_Metric.Metric[0][3] - s_dr_Metric.Metric[0][0] * s_dr_Metric.Metric[3][3]);

    double Kepler = this->Keplerian_angular_velocity(r, Spacetime);

    return  - Kepler / s_dr_Metric.Metric[3][3] * s_d2r_Metric.Metric[3][3] + (-s_d2r_Metric.Metric[0][3]
            + 1.0 / root / 2 * (2 * s_dr_Metric.Metric[0][3] * s_d2r_Metric.Metric[0][3] - s_dr_Metric.Metric[0][0] * s_d2r_Metric.Metric[3][3]
            - s_d2r_Metric.Metric[0][0] * s_dr_Metric.Metric[3][3])) / s_dr_Metric.Metric[3][3];

}

double Novikov_Thorne_Model::Redshift(double J, double State_Vector[], double r_obs, double theta_obs, Spacetime_Base_Class* Spacetime) {

    double& r_source = State_Vector[e_r];
    double& theta_source = State_Vector[e_theta];

    /*
    Get the observer 4-velocity
    */

    double State_Vector_obs[2] = { r_obs, theta_obs };

    Metric_type s_Metric_obs = Spacetime->get_metric(State_Vector_obs);

    double U_obs[4] = { 1.0 / s_Metric_obs.Lapse_function, 0 ,0 , s_Metric_obs.Shift_function / s_Metric_obs.Lapse_function };

    /*
    Get the source 4-velocity
    */

    Metric_type s_Metric_source = Spacetime->get_metric(State_Vector);

    double Kepler = this->Keplerian_angular_velocity(r_source, Spacetime);

    double Gamma = 1 / sqrt(-s_Metric_source.Metric[0][0] - 2 * s_Metric_source.Metric[0][3] * Kepler - s_Metric_source.Metric[3][3] * Kepler * Kepler);

    if (isnan(Gamma) || isinf(Gamma) || isnan(Kepler) || isinf(Kepler)) {

        std::cout << "Invalid NT disk 4-velocity: "
                  << "Gamma = "
                  << Gamma
                  << "\n"
                  << "Kepler = "
                  << "\n"
                  << Kepler
                  << "\n";

        exit(ERROR);

    }

    double U_source[4] = { Gamma, 0, 0, Gamma * Kepler };

    return  (-U_obs[0] + U_obs[3] * J) / (-U_source[0] + U_source[3] * J);

}

double Novikov_Thorne_Model::disk_Energy(double r, Spacetime_Base_Class* Spacetime) {

    double State_Vector[2] = { r, M_PI_2 };

    Metric_type s_Metric_source = Spacetime->get_metric(State_Vector);

    double Kepler = this->Keplerian_angular_velocity(r, Spacetime);

    double root = sqrt(-s_Metric_source.Metric[0][0] - 2 * s_Metric_source.Metric[0][3] * Kepler - s_Metric_source.Metric[3][3] * Kepler * Kepler);

    return  -(s_Metric_source.Metric[0][0] + s_Metric_source.Metric[0][3] * Kepler) / root;

}

double Novikov_Thorne_Model::disk_Angular_Momentum(double r, Spacetime_Base_Class* Spacetime) {

    double State_Vector[2] = { r, M_PI_2 };

    Metric_type s_Metric_source = Spacetime->get_metric(State_Vector);

    double Kepler = this->Keplerian_angular_velocity(r, Spacetime);

    double root = sqrt(-s_Metric_source.Metric[0][0] - 2 * s_Metric_source.Metric[0][3] * Kepler - s_Metric_source.Metric[3][3] * Kepler * Kepler);

    return  (s_Metric_source.Metric[3][3] * Kepler + s_Metric_source.Metric[0][3]) / root;

}

double Novikov_Thorne_Model::Flux_integrand(double r, Spacetime_Base_Class* Spacetime) {

    double State_Vector[2] = { r, M_PI_2 };

    Metric_type s_Metric = Spacetime->get_metric(State_Vector);
    Metric_type s_dr_Metric = Spacetime->get_dr_metric(State_Vector);

    double Kepler = this->Keplerian_angular_velocity(r, Spacetime);
    double dr_Kepler = this->dr_Keplerian_angular_velocity(r, Spacetime);

    double root = sqrt(-s_Metric.Metric[0][0] - 2 * s_Metric.Metric[0][3] * Kepler - s_Metric.Metric[3][3] * Kepler * Kepler);
    double dr_root = (-s_dr_Metric.Metric[0][0] - 2 * (s_dr_Metric.Metric[0][3] * Kepler + s_Metric.Metric[0][3] * dr_Kepler)
        - s_dr_Metric.Metric[3][3] * Kepler * Kepler - 2 * s_Metric.Metric[3][3] * Kepler * dr_Kepler);

    double E = this->disk_Energy(r, Spacetime);
    double L = this->disk_Angular_Momentum(r, Spacetime);

    double dr_L = (s_dr_Metric.Metric[3][3] * Kepler + s_Metric.Metric[3][3] * dr_Kepler + s_dr_Metric.Metric[0][3]) / root - L / root / root / 2 * dr_root;

    return (E - Kepler * L) * dr_L;

}

double Novikov_Thorne_Model::solve_Flux_integral(double lower_bound, double upper_bound, double tolerance, Spacetime_Base_Class* Spacetime) {

    double mid_point = (lower_bound + upper_bound) / 2;
    double left_mid_point = (lower_bound + mid_point) / 2;
    double right_mid_point = (mid_point + upper_bound) / 2;

    double F_lower_bound = this->Flux_integrand(lower_bound, Spacetime);
    double F_mid_point = this->Flux_integrand(mid_point, Spacetime);
    double F_upper_bound = this->Flux_integrand(upper_bound, Spacetime);

    double F_left_mid = this->Flux_integrand(left_mid_point, Spacetime);
    double F_right_mid = this->Flux_integrand(right_mid_point, Spacetime);

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

        double L_value = this->solve_Flux_integral(lower_bound, mid_point, tolerance / 2, Spacetime);
        double R_value = this->solve_Flux_integral(mid_point, upper_bound, tolerance / 2, Spacetime);

        integral = L_value + R_value;

    }

    return integral;
}

double Novikov_Thorne_Model::get_flux(double r, Spacetime_Base_Class* Spacetime) {

    double State_Vector[2] = { r, M_PI_2 };

    if (e_metric == Wormhole) {

        r = sqrt(r * r + WH_R_THROAT * WH_R_THROAT);

    }

    Metric_type s_Metric = Spacetime->get_metric(State_Vector);

    double metric_det = get_eq_induced_metric_det(s_Metric.Metric);
    double E_disk = disk_Energy(r, Spacetime);
    double L_disk = disk_Angular_Momentum(r, Spacetime);

    double Kepler = Keplerian_angular_velocity(r, Spacetime);
    double dr_Kepler = dr_Keplerian_angular_velocity(r, Spacetime);

    double Flux_coeff = -dr_Kepler / ((E_disk - Kepler * L_disk) * (E_disk - Kepler * L_disk)) / (4 * M_PI * sqrt(-metric_det));

    double Flux_integral = solve_Flux_integral(r_in, r, INTEGRAL_ACCURACY, Spacetime);


    return Flux_coeff * Flux_integral;

}

double Novikov_Thorne_Model::get_r_in() { return r_in; };
double Novikov_Thorne_Model::get_r_out() { return r_out; };

/************************************************************
|                                                           |
| Optically Thin Toroidal Model Class Functions Definitions |
|                                                           |
************************************************************/

/* ==================================================== Disk Temperature Functions ===================================================== */

double Generic_Optically_Thin_Model::get_disk_temperature(double State_vector[]) {


    double& r = State_vector[e_r];
    double Radial_Cutoff{};

    double T = T_ELECTRON_EXACT_CGS * this->s_Disk_Parameters.Power_law_radial_scale / r;

    if (r < this->s_Disk_Parameters.Disk_r_cutoff) {

        Radial_Cutoff = (r - this->s_Disk_Parameters.Disk_r_cutoff) / this->s_Disk_Parameters.Disk_cutoff_scale;

        T *= exp(-Radial_Cutoff * Radial_Cutoff);

    }

    this->Disk_Temperature = T;

    if (T < 0 || isnan(T) || isinf(T)) {

        std::cout << "Invalid Disk Temperature: " << this->Disk_Temperature << "\n";

        exit(ERROR);

    }

    return this->Disk_Temperature;

}

/* ====================================================== Disk Velocity Functions ====================================================== */

double* Generic_Optically_Thin_Model::get_disk_velocity(double State_Vector[], Simulation_Context_type* p_Sim_Context) {

    double& r_source = State_Vector[e_r];
    double& theta_source = State_Vector[e_theta];

    Metric_type s_Metric = p_Sim_Context->p_Spacetime->get_metric(State_Vector);

    double rho = r_source * sin(theta_source);

    double sqrt_rho = sqrt(rho);
    double ell = sqrt_rho * sqrt_rho * sqrt_rho / (1 + rho);

    if (e_metric == Naked_Singularity) {

        double gamma = p_Sim_Context->p_Spacetime->get_parameters().JNW_Gamma_Parameter;
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

    this->Disk_velocity[e_t_coord] = inv_metric[0][0] * u_t + inv_metric[0][3] * u_phi;
    this->Disk_velocity[e_r_coord] = 0.0;
    this->Disk_velocity[e_theta_coord] = 0.0;
    this->Disk_velocity[e_phi_coord] = inv_metric[3][3] * u_phi + inv_metric[3][0] * u_t;

    if (isnan(this->Disk_velocity[e_t_coord]) ||
        isinf(this->Disk_velocity[e_t_coord]) ||
        isnan(this->Disk_velocity[e_phi_coord]) ||
        isinf(this->Disk_velocity[e_phi_coord])) {

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

/* ======================================================= Disk Density Functions ====================================================== */

double Generic_Optically_Thin_Model::get_disk_hotspot_density_profile(double State_Vector[]) {

    double& Hotspot_r = this->s_Disk_Parameters.Hotspot_position[e_r];
    double& Hotspot_theta = this->s_Disk_Parameters.Hotspot_position[e_theta];
    double& Hotspot_phi = this->s_Disk_Parameters.Hotspot_position[e_phi];

    double x_center = Hotspot_r * sin(Hotspot_theta) * cos(Hotspot_phi);
    double y_center = Hotspot_r * sin(Hotspot_theta) * sin(Hotspot_phi);
    double z_center = Hotspot_r * cos(Hotspot_theta);

    double& r = State_Vector[e_r];

    double x_photon = r * sin(State_Vector[e_theta]) * cos(State_Vector[e_phi]);
    double y_photon = r * sin(State_Vector[e_theta]) * sin(State_Vector[e_phi]);
    double z_photon = r * cos(State_Vector[e_theta]);

    double hotspot_density_profile = exp(-(x_center - x_photon) * (x_center - x_photon) / HOTSPOT_SCALE / HOTSPOT_SCALE);
    hotspot_density_profile *= exp(-(y_center - y_photon) * (y_center - y_photon) / HOTSPOT_SCALE / HOTSPOT_SCALE);
    hotspot_density_profile *= exp(-(z_center - z_photon) * (z_center - z_photon) / HOTSPOT_SCALE / HOTSPOT_SCALE);
    hotspot_density_profile *= HOTSPOT_REL_SCALE;

    this->Hotspot_density_profile = hotspot_density_profile;

    if (isnan(Hotspot_density_profile) ||
        isinf(Hotspot_density_profile) ||
              Hotspot_density_profile < 0) {

        std::cout << "Invalid hotspot density profile: " << this->Hotspot_density_profile << "\n";

        exit(ERROR);

    }

    return this->Hotspot_density_profile;


}

double Generic_Optically_Thin_Model::get_disk_density_profile(double State_Vector[]) {

    double& r  = State_Vector[e_r];
    double rho = sin(State_Vector[e_theta]);
    double h   = cos(State_Vector[e_theta]);

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

        std::cout << "Unsupported disk density model!" << '\n';

        exit(ERROR);

    }

    if (this->s_Disk_Parameters.Hotspot_scale != 0.0f) {

        this->Disk_density_profile = electron_density_profile + get_disk_hotspot_density_profile(State_Vector);

    }
    else {

        this->Disk_density_profile = electron_density_profile;

    }

    if (isnan(Disk_density_profile) ||
        isinf(Disk_density_profile) ||
              Disk_density_profile < 0) {

        std::cout << "Invalid disk density profile: " << this->Disk_density_profile << "\n";

        exit(ERROR);

    }

    return this->Disk_density_profile;

}

/* =================================================== Disk Magnetic Field Functions =================================================== */

double Generic_Optically_Thin_Model::get_magnetic_field(double B_field[4],
                                                        double State_Vector[],
                                                        Simulation_Context_type* p_Sim_Context) {

    /*
    
    Computes the magnetic field 4-vector, measured by a comoving obverver (with 4-velocity Plasma_velocity), called B_field,
    in the basis of a static observer (with 4-velocity n_mu = {1, 0, 0, 0} ). 
    Returns the magnitude of the magnetic field in the plasma frame.

    NOTE: The magnitude of the magnetic field in these frames is different, because its not concerved under Lorentz boosts.
    In the static frame I only set the geometry of the field, then finally scale B_field (the plasma frame one) by B_Plasma_norm_CGS.
    
    NOTE: B_field is measured by a static observer, so:
    1) B_field is effectively in the coordinate basis (not any ZAMO!) - this means one cannot directly apply Lorentz boosts to it.

    */

    double electron_density = this->s_Disk_Parameters.Density_scale * this->get_disk_density_profile(State_Vector);

    double B_Plasma_norm_CGS = sqrt(this->s_Disk_Parameters.Magnetization * C_LIGHT_CGS * C_LIGHT_CGS * electron_density * M_PROTON_CGS * 4 * M_PI);
    double B_Static[4] = {         0.0, 
                           MAG_FIELD_GEOMETRY[0],
                           MAG_FIELD_GEOMETRY[1],
                           MAG_FIELD_GEOMETRY[2] };

    double* Plasma_velocity = this->get_disk_velocity(State_Vector, p_Sim_Context);

    /* 

    The only reference I could find for this is https://iopscience.iop.org/article/10.3847/1538-4357/ab718e/pdf 8e) and 8f).
    One arrives at the expression by projecting F^mu^nu (written in terms of B_Plasma and Plasma_veclocity),
    onto the 4-velocity of the static observer - this gives B_Static, then inverting the expression 
    to obtain B_Plasma = f(B_Static, Plasma_Velocity)
    
    */

    Metric_type s_Metric = p_Sim_Context->p_Spacetime->get_metric(State_Vector);

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

double Generic_Optically_Thin_Model::get_electron_pitch_angle(double B_field_local[4], 
                                                              double State_Vector[], 
                                                              Simulation_Context_type* p_Sim_Context) {

    /*
    
    Computes the angle between the photon wave-3-vector and the magnetic field, measured by a comoving (with the plasma) observer.
    There is a neat invariant way to compute this directly from the coordinate basis 4-vectors.
    
    */

    double* Plasma_velocity = this->get_disk_velocity(State_Vector, p_Sim_Context);
    double Wave_vec_dot_Plasma_vec =           -1           * Plasma_velocity[e_t_coord] + 
                                    State_Vector[e_p_r]     * Plasma_velocity[e_r_coord] +
                                    State_Vector[e_p_theta] * Plasma_velocity[e_theta_coord] +
                                    State_Vector[e_p_phi]   * Plasma_velocity[e_phi_coord];

    Metric_type s_Metric  = p_Sim_Context->p_Spacetime->get_metric(State_Vector);
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

/* =============================================== Thermal Synchotron Transfer Functions =============================================== */

void Generic_Optically_Thin_Model::get_thermal_synchotron_faradey_fit_functions(double X,
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

void Generic_Optically_Thin_Model::get_thermal_synchotron_emission_fit_functions(Thermal_Syncotron_fit_selector e_fit_functions,
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

void Generic_Optically_Thin_Model::get_thermal_synchotron_fit_functions(double Emission_fucntions[STOKES_PARAM_NUM],
                                                                        double Faradey_functions[STOKES_PARAM_NUM],
                                                                        Thermal_emission_f_arguments* Emission_args,
                                                                        Thermal_faradey_f_arguments* Faradey_args) {

    if (INCLUDE_POLARIZATION){

        this->get_thermal_synchotron_emission_fit_functions(Dexter_2016, Emission_fucntions, Emission_args->X, Emission_args->sqrt_X, Emission_args->cbrt_X);
        this->get_thermal_synchotron_faradey_fit_functions(Faradey_args->X, Faradey_args->X_to_1_point_2, Faradey_args->X_to_1_point_035, Faradey_functions);

    }
    else {

        this->get_thermal_synchotron_emission_fit_functions(Leung_2011, Emission_fucntions, Emission_args->X, Emission_args->sqrt_X, Emission_args->cbrt_X);

    }
}


void Generic_Optically_Thin_Model::get_thermal_synchotron_transfer_functions(double State_Vector[],
                                                                             Simulation_Context_type* p_Sim_Context,
                                                                             double Emission_fucntions[STOKES_PARAM_NUM],
                                                                             double Faradey_functions[STOKES_PARAM_NUM],
                                                                             double Absorbtion_functions[STOKES_PARAM_NUM]) {

    /* === The transfer functions arrays needs to be manually cleared, because this function only adds to it. === */
    for (int stokes_index = 0; stokes_index <= STOKES_PARAM_NUM - 1; stokes_index++) {

        Emission_fucntions[stokes_index]   = 0.0;
        Faradey_functions[stokes_index]    = 0.0;
        Absorbtion_functions[stokes_index] = 0.0;

    }

    /* Electron Density in CGS */
    double electron_density = this->s_Disk_Parameters.Density_scale * this->get_disk_density_profile(State_Vector);

    /* Dimentionless Electron Temperature */
    double T_electron = this->get_disk_temperature(State_Vector);
    double T_electron_dim = BOLTZMANN_CONST_CGS * T_electron / M_ELECTRON_CGS / C_LIGHT_CGS / C_LIGHT_CGS;
    double K0_Bessel = std::cyl_bessel_k(0.0, 1.0 / T_electron_dim);
    double K1_Bessel = std::cyl_bessel_k(1.0, 1.0 / T_electron_dim);
    double K2_Bessel = std::cyl_bessel_k(2.0, 1.0 / T_electron_dim);

    /* Magnetic Field */
    double B_field_local[4]{};
    double B_CGS = this->get_magnetic_field(B_field_local, State_Vector, p_Sim_Context);

    /* Redshit */
    double redshift = Redshift(State_Vector, this->get_disk_velocity(State_Vector, p_Sim_Context), p_Sim_Context->p_Observer);

    /* Cyclotron Frequency */
    double f_cyclo = Q_ELECTRON_CGS * B_CGS / (2 * M_PI * M_ELECTRON_CGS * C_LIGHT_CGS);

    /* The "averaged" critical frequency (without the sin(theta) term - that gets added on later from a pre-computed table) */
    double f_crit_no_sin = 3. / 2 * f_cyclo * T_electron_dim * T_electron_dim;

    Thermal_emission_f_arguments Emission_args_ang_uncorrected{};
    Thermal_faradey_f_arguments Faradey_args_ang_uncorrected{};

    /* Both the emission and faradey function expressions are in terms of an dimentionless variable X, but the definitions for X are different */
    Emission_args_ang_uncorrected.X = 1e100;
    Faradey_args_ang_uncorrected.X  = 1e100;

    if (f_crit_no_sin > std::numeric_limits<double>::min()) {

        Emission_args_ang_uncorrected.X = OBS_FREQUENCY_CGS / f_crit_no_sin / redshift;
        Faradey_args_ang_uncorrected.X  = T_electron_dim * sqrt(M_SQRT2 * 1e3 * f_cyclo / (OBS_FREQUENCY_CGS / redshift));

    }

    /* All the functions are normalized by a Bessel function, so I check if I can divide by it */
    if (!isinf(1e10 / K2_Bessel)) {

       /* == Compute all the weird powers of X outside the averaging loop == */

        Emission_args_ang_uncorrected.sqrt_X = sqrt(Emission_args_ang_uncorrected.X);
        Emission_args_ang_uncorrected.cbrt_X = cbrt(Emission_args_ang_uncorrected.X);

        Faradey_args_ang_uncorrected.X_to_1_point_2   = pow(Faradey_args_ang_uncorrected.X, 1.2f);
        Faradey_args_ang_uncorrected.X_to_1_point_035 = pow(Faradey_args_ang_uncorrected.X, 1.035f);

       /* ================================================================== */

        if (AVERAGE_EMISSION_PITCH_ANGLE) {

            for (int averaging_idx = 0; averaging_idx <= NUM_SAMPLES_TO_AVG - 1; averaging_idx++) {

                Thermal_emission_f_arguments Emission_args_ang_corrected{};
                Thermal_faradey_f_arguments Faradey_args_ang_corrected{};

                double& sin_pitch_angle = this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[averaging_idx];
                double& cos_pitch_angle = this->s_Precomputed_e_pitch_angles.cos_electron_pitch_angles[averaging_idx];

                double temp_emission_functions[STOKES_PARAM_NUM]{};
                double temp_faradey_functions[STOKES_PARAM_NUM]{};

                Emission_args_ang_corrected.X      = Emission_args_ang_uncorrected.X / sin_pitch_angle;
                Emission_args_ang_corrected.sqrt_X = Emission_args_ang_uncorrected.sqrt_X * this->s_Precomputed_e_pitch_angles.one_over_sqrt_sin[averaging_idx];
                Emission_args_ang_corrected.cbrt_X = Emission_args_ang_uncorrected.cbrt_X * this->s_Precomputed_e_pitch_angles.one_over_cbrt_sin[averaging_idx];

                Faradey_args_ang_corrected.X                = Faradey_args_ang_uncorrected.X / this->s_Precomputed_e_pitch_angles.one_over_sqrt_sin[averaging_idx];
                Faradey_args_ang_corrected.X_to_1_point_035 = Faradey_args_ang_uncorrected.X_to_1_point_035 / this->s_Precomputed_e_pitch_angles.one_over_sin_to_1_point_035[averaging_idx];
                Faradey_args_ang_corrected.X_to_1_point_2   = Faradey_args_ang_uncorrected.X_to_1_point_2 / this->s_Precomputed_e_pitch_angles.one_over_sin_to_1_point_2_over_2[averaging_idx];

                this->get_thermal_synchotron_fit_functions(temp_emission_functions, temp_faradey_functions, &Emission_args_ang_corrected, &Faradey_args_ang_corrected);

                /* ================================================ The emission functions ================================================ */

                temp_emission_functions[I] *= electron_density * (f_crit_no_sin * sin_pitch_angle) / K2_Bessel; // Scale the emission by the remaining position-dependant factors
                temp_emission_functions[I] *= sin_pitch_angle * M_PI / NUM_SAMPLES_TO_AVG / 2;                  // Scale by the factors coming form averaging over all emission orientations

                temp_emission_functions[Q] *= electron_density * (f_crit_no_sin * sin_pitch_angle) / K2_Bessel; // Scale the emission by the remaining position-dependant factors
                temp_emission_functions[Q] *= sin_pitch_angle * M_PI / NUM_SAMPLES_TO_AVG / 2;                  // Scale by the factors coming form averaging over all emission orientations

                if (!isinf(1e2 / sin_pitch_angle) && !isnan(1e2 / sin_pitch_angle)) {

                    temp_emission_functions[V] *= electron_density * (f_crit_no_sin * sin_pitch_angle) / K2_Bessel; // Scale the emission by the remaining position-dependant factors
                    temp_emission_functions[V] *= (1. / T_electron_dim) * cos_pitch_angle / sin_pitch_angle;        // The V component has some extra angle dependance
                    temp_emission_functions[V] *= sin_pitch_angle * M_PI / NUM_SAMPLES_TO_AVG / 2;                  // Scale by the factors coming form averaging over all emission orientations

                }

                Emission_fucntions[I] += temp_emission_functions[I];
                Emission_fucntions[Q] += temp_emission_functions[Q];
                Emission_fucntions[U] = 0.0;
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

                Faradey_functions[I] = 0.0f;
                Faradey_functions[Q] += temp_faradey_functions[Q];
                Faradey_functions[U] = 0.0f;
                Faradey_functions[V] += temp_faradey_functions[V];

            }
        }
        else {
            
            double cos_pitch_angle = get_electron_pitch_angle(B_field_local, State_Vector, p_Sim_Context);
            double sin_pitch_angle = sqrt(1.0 - cos_pitch_angle * cos_pitch_angle);

            double one_over_sqrt_sin = 1.0 / sqrt(sin_pitch_angle);
            double one_over_cbrt_sin = 1.0 / cbrt(sin_pitch_angle);

            Thermal_emission_f_arguments Emission_args_ang_corrected{};
            Thermal_faradey_f_arguments Faradey_args_ang_corrected{};

            Emission_args_ang_corrected.X      = Emission_args_ang_uncorrected.X / sin_pitch_angle;
            Emission_args_ang_corrected.sqrt_X = Emission_args_ang_uncorrected.sqrt_X * one_over_sqrt_sin;
            Emission_args_ang_corrected.cbrt_X = Emission_args_ang_uncorrected.cbrt_X * one_over_cbrt_sin;

            Faradey_args_ang_corrected.X                = Faradey_args_ang_uncorrected.X / one_over_sqrt_sin;
            Faradey_args_ang_corrected.X_to_1_point_035 = Faradey_args_ang_uncorrected.X_to_1_point_035 * pow(sin_pitch_angle, 1.035);
            Faradey_args_ang_corrected.X_to_1_point_2   = Faradey_args_ang_uncorrected.X_to_1_point_2 * pow(sin_pitch_angle, 1.2);

            this->get_thermal_synchotron_fit_functions(Emission_fucntions, Faradey_functions, &Emission_args_ang_corrected, &Faradey_args_ang_corrected);

            /* ================================================ The emission functions ================================================ */

            Emission_fucntions[I] *= electron_density * (f_crit_no_sin * sin_pitch_angle) / K2_Bessel; // Scale the emission by the remaining position-dependant factors
            Emission_fucntions[Q] *= electron_density * (f_crit_no_sin * sin_pitch_angle) / K2_Bessel; // Scale the emission by the remaining position-dependant factors

            if (!isinf(1e2 / sin_pitch_angle) && !isnan(1e2 / sin_pitch_angle)) {

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

    double Planck_function_CGS = get_planck_function_CGS(OBS_FREQUENCY_CGS / redshift, this->get_disk_temperature(State_Vector));

    if (Planck_function_CGS > std::numeric_limits<double>::min()) {

        for (int stokes_index = 0; stokes_index <= STOKES_PARAM_NUM - 1; stokes_index++) {

            Absorbtion_functions[stokes_index] = Emission_fucntions[stokes_index] / Planck_function_CGS;

        }
    }
}

/* ========================================== Kappa Synchotron Transfer Functions ========================================== */

void Generic_Optically_Thin_Model::get_kappa_synchotron_emission_fit_functions(double Emission_functions[STOKES_PARAM_NUM],
                                                                               double X,
                                                                               double sqrt_X,
                                                                               double cbrt_X,
                                                                               double X_to_7_over_20,
                                                                               double kappa,
                                                                               double sin_emission_angle,
                                                                               double T_electron_dim){

    // The reference for these expressions is https://arxiv.org/pdf/1602.08749, equations (35), (36), (37) and (38).

    // ------------------------------------------------------------------------ Low frequency fit ------------------------------------------------------------------------ //

    constexpr double THREE_TO_7_OVER_3 = 12.980246132766677;

    double Emission_functions_low[STOKES_PARAM_NUM]{};
    double Common_factor_low = cbrt_X * sin_emission_angle * (4 * M_PI / THREE_TO_7_OVER_3) * std::tgamma(kappa - 4.0 / 3) / std::tgamma(kappa - 2);

    Emission_functions_low[I] =  Common_factor_low;
    Emission_functions_low[Q] = -Common_factor_low / 2;
    Emission_functions_low[U] = 0;
    Emission_functions_low[V] = -Common_factor_low * (9.0 / 16 * pow(pow(sin_emission_angle, -12.0 / 5) - 1, 12.0 / 25)) * pow(kappa, -66.0/125) / T_electron_dim / X_to_7_over_20;

    // ----------------------------------------------------------------------- High frequency fit ------------------------------------------------------------------------ //

    double Emission_functions_high[STOKES_PARAM_NUM]{};
    double Common_factor_high = pow(X, -(kappa - 2) / 2) * sin_emission_angle * pow(3, (kappa - 1) / 2) * (kappa - 2) * (kappa - 1) / 4 * std::tgamma(kappa / 4 - 1.0 / 3) * std::tgamma(kappa / 4 + 4.0 / 3);

    Emission_functions_high[I] =  Common_factor_high;
    Emission_functions_high[Q] = -Common_factor_high * (16.0 / 25 + kappa / 50);
    Emission_functions_high[U] = 0;
    Emission_functions_high[V] = -Common_factor_high * (49.0 / 64 * pow(pow(sin_emission_angle, -5.0 / 2) - 1, 11.0 / 25)) * pow(kappa, -11.0 / 25) / T_electron_dim / sqrt_X;

    // ------------------------------------------------------------------------ Bridging function ------------------------------------------------------------------------ //

    double power_I = 3 * pow(kappa, -3.0 / 2);
    Emission_functions[I] = pow(pow(Emission_functions_low[I], -power_I) + pow(Emission_functions_high[I], -power_I), -1.0 / power_I);

    double power_Q = 3.7 * pow(kappa, -8.0 / 5);
    Emission_functions[Q] = pow(pow(Emission_functions_low[Q], -power_Q) + pow(Emission_functions_high[Q], -power_Q), -1.0 / power_Q);

    Emission_functions[U] = 0;

    double power_V = 13.0 / 5 * pow(kappa, -36.0 / 25);
    Emission_functions[V] = pow(pow(Emission_functions_low[V], -power_V) + pow(Emission_functions_high[V], -power_V), -1.0 / power_V);

}

void Generic_Optically_Thin_Model::get_kappa_synchotron_absorbtion_fit_functions(double Absorbtion_functions[STOKES_PARAM_NUM],
                                                                                 double X,
                                                                                 double sqrt_X,
                                                                                 double cbrt_X,
                                                                                 double X_to_7_over_20,
                                                                                 double kappa,
                                                                                 double sin_emission_angle,
                                                                                 double T_electron_dim){

    // The reference for these expressions is https://arxiv.org/pdf/1602.08749, equations (39), (40), (41) and (42).

    // ------------------------------------------------------------------------ Low frequency fit ------------------------------------------------------------------------ //

    constexpr double THREE_TO_1_OVER_6 = 1.2009369551760027;
    constexpr double GAMMA_OF_5_OVER_3 = 0.90274529295;

    // Below are the coefficients, present in the 2F1 hypergoemetric function from expression (39). Bcause in general |(-kappa * T_electron_dim)| > 1, I will use the algebraic relation 
    // from Abramowitz and Stegun 15.3.8. Note that our equivalent to the argument z in 15.3.8 is strictly real and negative, so taking fractional powers of it returns a real number.

    double a =  kappa - 1.0 / 3;
    double b =  kappa + 1.0;
    double c =  kappa + 2.0 / 3;
    double z = -kappa * T_electron_dim;

    double _2F1 = pow(1 - z, -a) * std::tgamma(c) / std::tgamma(b) * std::tgamma(b - a) / std::tgamma(c - a) * gsl_sf_hyperg_2F1(a, c - b, a - b + 1, 1.0 / (1 - z))
                + pow(1 - z, -b) * std::tgamma(c) / std::tgamma(a) * std::tgamma(a - b) / std::tgamma(c - b) * gsl_sf_hyperg_2F1(b, c - a, b - a + 1, 1.0 / (1 - z));

    double Absorbtion_functions_low[STOKES_PARAM_NUM]{};
    double Common_factor_low =  1.0 / cbrt_X / cbrt_X * THREE_TO_1_OVER_6 * 10.0 / 41 * 2 * M_PI / pow(T_electron_dim * kappa, 10.0 / 3 - kappa) * (kappa - 2) * (kappa - 1) * kappa / (3 * kappa - 1) 
                             * GAMMA_OF_5_OVER_3 * _2F1;

    Absorbtion_functions_low[I] =  Common_factor_low;
    Absorbtion_functions_low[Q] = -Common_factor_low * 25.0 / 48;
    Absorbtion_functions_low[U] = 0;
    Absorbtion_functions_low[V] = -Common_factor_low * pow((pow(sin_emission_angle, -114.0 / 50) - 1), 223.0 / 500) / X_to_7_over_20 * pow(kappa, -7.0 / 10);

    // ----------------------------------------------------------------------- High frequency fit ------------------------------------------------------------------------ //

    double Absorbtion_functions_high[STOKES_PARAM_NUM]{};
    double Common_factor_high = pow(X, -(1 + kappa) / 2) * M_PI * (2 / M_2_SQRTPI) / 3 * (kappa - 2) * (kappa - 1) * kappa / (kappa * T_electron_dim) / (kappa * T_electron_dim) / (kappa * T_electron_dim)
                              * (2 * std::tgamma(2 + kappa / 2) / (2 + kappa) - 1);

    Absorbtion_functions_high[I] =  Common_factor_low * (pow(3.0 / kappa, 19.0 / 4) + 3.0 / 5);
    Absorbtion_functions_high[Q] = -Common_factor_low * (441 * pow(kappa, - 144.0/ 25) + 11.0 / 20);
    Absorbtion_functions_high[U] = 0;
    Absorbtion_functions_high[V] = -Common_factor_low * 143.0 / 10 * pow(T_electron_dim, -116.0 / 125) * sqrt(pow(sin_emission_angle, -41.0 / 20) - 1) * (169 * pow(kappa, -8) + 13.0 / 2500 * kappa - 1.0 / 200 + 47.0 / 200 / kappa) / sqrt_X;

    // ------------------------------------------------------------------------ Bridging function ------------------------------------------------------------------------ //

    double power_I = pow(-7.0 / 4 + 8.0 / 5 * kappa, -43.0 / 50);
    Absorbtion_functions[I] = pow(pow(Absorbtion_functions_low[I], -power_I) + pow(Absorbtion_functions_high[I], -power_I), -1.0 / power_I);

    double power_Q = 7.0 / 5 * pow(kappa, -23.0 / 20);
    Absorbtion_functions[Q] = pow(pow(Absorbtion_functions_low[Q], -power_Q) + pow(Absorbtion_functions_high[Q], -power_Q), -1.0 / power_Q);

    Absorbtion_functions[U] = 0;

    double power_V = 61.0 / 50 * pow(kappa, -142.0 / 125) + 7.0 / 1000;
    Absorbtion_functions[V] = pow(pow(Absorbtion_functions_low[V], -power_V) + pow(Absorbtion_functions_high[V], -power_V), -1.0 / power_V);


}

/* ========================================== Phenomenological Synchotron Transfer Functions ========================================== */

void Generic_Optically_Thin_Model::get_phenomenological_synchotron_functions(double State_Vector[],
                                                                             Simulation_Context_type* p_Sim_Context, 
                                                                             double Emission_functions[STOKES_PARAM_NUM],
                                                                             double Faradey_functions[STOKES_PARAM_NUM],
                                                                             double Absorbtion_functions[STOKES_PARAM_NUM]) {

    double& emission_power_law = this->s_Emission_Parameters.Emission_power_law;
    double& source_f_power_law = this->s_Emission_Parameters.Source_f_power_law;
    double& emission_scale     = this->s_Emission_Parameters.Emission_scale;
    double& abs_coeff          = this->s_Emission_Parameters.Absorbtion_coeff;

    double* U_source_coord = this->get_disk_velocity(State_Vector, p_Sim_Context);
    double redshift = Redshift(State_Vector, U_source_coord, p_Sim_Context->p_Observer);

    Emission_functions[I] = emission_scale * get_disk_density_profile(State_Vector) * pow(redshift, emission_power_law);
    Emission_functions[Q] = 0.0;
    Emission_functions[U] = 0.0;
    Emission_functions[V] = 0.0;

    Absorbtion_functions[I] = abs_coeff * emission_scale * this->get_disk_density_profile(State_Vector) * pow(redshift, source_f_power_law + emission_power_law);
    Absorbtion_functions[Q] = 0.0;
    Absorbtion_functions[U] = 0.0;
    Absorbtion_functions[V] = 0.0;

    Faradey_functions[I] = 0.0;
    Faradey_functions[Q] = 0.0;
    Faradey_functions[U] = 0.0;
    Faradey_functions[V] = 0.0;

}

/* ============================================ Main "Selector" For The Transfer Functions ============================================ */

void Generic_Optically_Thin_Model::get_radiative_transfer_functions(double State_Vector[e_State_Number],
                                                                    Simulation_Context_type* p_Sim_Context, 
                                                                    double Emission_functions[STOKES_PARAM_NUM],
                                                                    double Faradey_functions[STOKES_PARAM_NUM],
                                                                    double Absorbtion_functions[STOKES_PARAM_NUM]) {
    switch (e_emission) {

    case(Phenomenological_synchotron):

        this->get_phenomenological_synchotron_functions(State_Vector, p_Sim_Context, Emission_functions, Faradey_functions, Absorbtion_functions);
        break;

    default:

        this->get_thermal_synchotron_transfer_functions(State_Vector, p_Sim_Context, Emission_functions, Faradey_functions, Absorbtion_functions);
        break;
    }
}

/* ========================================================== Misc Functions ========================================================== */

void Generic_Optically_Thin_Model::precompute_electron_pitch_angles() {

    for (int index = 0; index <= NUM_SAMPLES_TO_AVG - 1; index++) {

        double pitch_angle = double(index) / NUM_SAMPLES_TO_AVG * M_PI;
        this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index] = sin(pitch_angle);
        this->s_Precomputed_e_pitch_angles.cos_electron_pitch_angles[index] = cos(pitch_angle);

        if (this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index] != 0) {

            // Used in the thermal synchotron emission functions

            this->s_Precomputed_e_pitch_angles.one_over_sqrt_sin[index] = 1. / sqrt(this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index]);
            this->s_Precomputed_e_pitch_angles.one_over_cbrt_sin[index] = 1. / cbrt(this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index]);

            // Used in the thermal synchotron Faradey functions

            this->s_Precomputed_e_pitch_angles.one_over_sin_to_1_point_035[index] = 1. / pow(this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index], 1.035f);
            this->s_Precomputed_e_pitch_angles.one_over_sin_to_1_point_2_over_2[index] = 1. / pow(this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index], 1.2f / 2);

        }
    }
}

int Generic_Optically_Thin_Model::load_parameters(Disk_model_parameters* p_Disk_Parameters, Emission_law_parameters* p_Emission_Parameters) {

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

double get_planck_function_CGS(double Frequency, double Temperature) {

    return 2 * PLANCK_CONSTANT_CGS * Frequency * Frequency * Frequency / C_LIGHT_CGS / C_LIGHT_CGS / (exp(PLANCK_CONSTANT_CGS * Frequency / BOLTZMANN_CONST_CGS / Temperature) - 1);

}
