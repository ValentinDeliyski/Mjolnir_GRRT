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

Novikov_Thorne_Model::Novikov_Thorne_Model(NT_parameters_type NT_params, Spacetime_Base_Class* Spacetime) {

    r_in = NT_params.r_in;

    if (NT_params.r_in < 0) {

        r_in = Spacetime->get_ISCO()[Outer];

    }

    r_out = NT_params.r_out;

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

/* ==================================================== Temperature Functions ===================================================== */

double Generic_Optically_Thin_Model::get_disk_temperature(double State_vector[]) {


    double& r            = State_vector[e_r];
    double& T_scale      = this->s_Disk_params.Electron_temperature_scale;
    double& R_0          = this->s_Disk_params.Power_law_temperature_R_0;
    double& r_cutoff     = this->s_Disk_params.Power_law_temperature_R_cutoff;
    double& Cutoff_scale = this->s_Disk_params.Power_law_temperature_cutoff_scale;
    double& Power_law    = this->s_Disk_params.Power_law_temperature_radial_power_law;

    double Radial_Cutoff{};

    double Disk_temperature_profile = 1.0 / pow(r / R_0, Power_law);

    if (r < r_cutoff) {

        Radial_Cutoff = (r - r_cutoff) / Cutoff_scale;

        Disk_temperature_profile *= exp(-Radial_Cutoff * Radial_Cutoff);

    }

    if (Disk_temperature_profile < 0 || isnan(Disk_temperature_profile) || isinf(Disk_temperature_profile)) {

        std::cout << "Invalid Disk Temperature profile: " << Disk_temperature_profile << "\n";

        exit(ERROR);

    }

    this->Disk_Temperature = T_scale * Disk_temperature_profile;

    return this->Disk_Temperature;

}

double Generic_Optically_Thin_Model::get_hotspot_temperature(double State_Vector[]) {

    double& Hotspot_r     = this->s_Hotspot_params.Position[e_r];
    double& Hotspot_theta = this->s_Hotspot_params.Position[e_theta];
    double& Hotspot_phi   = this->s_Hotspot_params.Position[e_phi];

    double& Hotspot_spread    = this->s_Hotspot_params.Temperature_spread;

    double& photon_r = State_Vector[e_r];

    // I these this more than once, so I precompute them
    double sin_theta = sin(State_Vector[e_theta]);
    double sin_hotspot_theta = sin(Hotspot_theta);

    double x_center = Hotspot_r * sin_hotspot_theta * cos(Hotspot_phi);
    double y_center = Hotspot_r * sin_hotspot_theta * sin(Hotspot_phi);
    double z_center = Hotspot_r * cos(Hotspot_theta);

    double x_photon = photon_r * sin_theta * cos(State_Vector[e_phi]);
    double y_photon = photon_r * sin_theta * sin(State_Vector[e_phi]);
    double z_photon = photon_r * cos(State_Vector[e_theta]);

    double exponent_argument = - (x_center - x_photon) * (x_center - x_photon) / Hotspot_spread / Hotspot_spread / 2
                               - (y_center - y_photon) * (y_center - y_photon) / Hotspot_spread / Hotspot_spread / 2
                               - (z_center - z_photon) * (z_center - z_photon) / Hotspot_spread / Hotspot_spread / 2;

    double Hotspot_temperature_profile = exp(exponent_argument);

    if (isnan(Hotspot_temperature_profile) || isinf(Hotspot_temperature_profile) || Hotspot_temperature_profile < 0) {

        std::cout << "Invalid hotspot temperature profile: " << Hotspot_temperature_profile << "\n";

        exit(ERROR);

    }

    this->Hotspot_Temperature = this->s_Hotspot_params.Electron_temperature_scale * Hotspot_temperature_profile;

    return this->Hotspot_Temperature;

}

/* ====================================================== Velocity Functions ====================================================== */

double* Generic_Optically_Thin_Model::get_disk_velocity(double State_Vector[], Simulation_Context_type* p_Sim_Context) {

    double& r_source = State_Vector[e_r];
    double& theta_source = State_Vector[e_theta];

    Metric_type s_Metric = p_Sim_Context->p_Spacetime->get_metric(State_Vector);

    double rho = r_source * sin(theta_source);

    double sqrt_rho = sqrt(rho);
    double ell = sqrt_rho * sqrt_rho * sqrt_rho / (1 + rho);

    if (Naked_Singularity == p_Sim_Context->e_Spacetime) {

        double gamma = p_Sim_Context->p_Init_Conditions->Metric_params.JNW_Gamma_Parameter;
        double r_singularity = 2. / gamma;

        ell *= pow(1. - r_singularity / r_source, gamma);

    }
    else if (Wormhole == p_Sim_Context->e_Spacetime) {


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
    //this->Disk_velocity[e_t_coord] = sqrt(r_source / (r_source - 3));
    this->Disk_velocity[e_r_coord] = 0.0;
    this->Disk_velocity[e_theta_coord] = 0.0;
    //this->Disk_velocity[e_phi_coord] = sqrt(r_source / (r_source - 3)) / pow(r_source, 3. / 2);
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

/* ======================================================= Density Functions ====================================================== */

double Generic_Optically_Thin_Model::get_hotspot_density(double State_Vector[]) {

    double& Hotspot_r     = this->s_Hotspot_params.Position[e_r];
    double& Hotspot_theta = this->s_Hotspot_params.Position[e_theta];
    double& Hotspot_phi   = this->s_Hotspot_params.Position[e_phi];

    double& Hotspot_spread = this->s_Hotspot_params.Density_spread;

    double& photon_r = State_Vector[e_r];

    // I need these more than once, so I precompute them here
    double sin_theta = sin(State_Vector[e_theta]);
    double sin_hotspot_theta = sin(Hotspot_theta);

    double x_center = Hotspot_r * sin_hotspot_theta * cos(Hotspot_phi);
    double y_center = Hotspot_r * sin_hotspot_theta * sin(Hotspot_phi);
    double z_center = Hotspot_r * cos(Hotspot_theta);

    double x_photon = photon_r * sin_theta * cos(State_Vector[e_phi]);
    double y_photon = photon_r * sin_theta * sin(State_Vector[e_phi]);
    double z_photon = photon_r * cos(State_Vector[e_theta]);

    double exponent_argument = -(x_center - x_photon) * (x_center - x_photon) / Hotspot_spread / Hotspot_spread / 2.
                               -(y_center - y_photon) * (y_center - y_photon) / Hotspot_spread / Hotspot_spread / 2. 
                               -(z_center - z_photon) * (z_center - z_photon) / Hotspot_spread / Hotspot_spread / 2.;

    double Hotspot_density_profile = exp(exponent_argument);
  
    if (isnan(Hotspot_density_profile) ||
        isinf(Hotspot_density_profile) ||
              Hotspot_density_profile < 0) {

        std::cout << "Invalid hotspot density profile: " << Hotspot_density_profile << "\n";

        exit(ERROR);

    }

    this->Hotspot_density = this->s_Hotspot_params.Electron_density_scale * Hotspot_density_profile;

    return this->Hotspot_density;

}

double Generic_Optically_Thin_Model::get_disk_density(double State_Vector[]) {

    double& r  = State_Vector[e_r];
    double rho = sin(State_Vector[e_theta]);
    double h   = cos(State_Vector[e_theta]);

    double Height_Cutoff{};
    double Radial_Cutoff{};
    double Disk_density_profile{};

    double& Opening_angle  = this->s_Disk_params.Power_law_disk_opening_angle;
    double& R_0            = this->s_Disk_params.Power_law_density_R_0;
    double& r_cutoff       = this->s_Disk_params.Power_law_density_R_cutoff;
    double& Cutoff_scale   = this->s_Disk_params.Power_law_density_cutoff_scale;
    double& Disk_power_law = this->s_Disk_params.Power_law_density_radial_power_law;

    double& Height_scale = this->s_Disk_params.Exp_law_density_height_scale;
    double& Radial_scale = this->s_Disk_params.Exp_law_density_radial_scale;

    switch (this->s_Disk_params.Density_profile_type) {

    case e_Power_law_profile:

        Height_Cutoff = h / (Opening_angle * rho);
        Disk_density_profile = exp(-Height_Cutoff * Height_Cutoff / 2) / pow(r / R_0, Disk_power_law);

        if (r < r_cutoff) {

            Radial_Cutoff = (r - r_cutoff) / Cutoff_scale;
            Disk_density_profile *= exp(-Radial_Cutoff * Radial_Cutoff);

        }

        break;

    case e_Exponential_law_profile:

        Height_Cutoff = h / Height_scale;
        Radial_Cutoff = r / Radial_scale;

        Disk_density_profile = exp(-Radial_Cutoff * Radial_Cutoff / 2 - Height_Cutoff * Height_Cutoff / 2);

        break;

    default:

        std::cout << "Unsupported disk density model!" << '\n';

        exit(ERROR);

    }

    if (isnan(Disk_density) ||
        isinf(Disk_density) ||
              Disk_density < 0) {

        std::cout << "Invalid disk density profile: " << Disk_density << "\n";

        exit(ERROR);

    }

    this->Disk_density = this->s_Disk_params.Electron_density_scale * Disk_density_profile;

    return this->Disk_density;

}

/* =================================================== Disk Magnetic Field Functions =================================================== */

double Generic_Optically_Thin_Model::get_total_magnetic_field(double B_coord_frame[4],
                                                              double State_Vector[],
                                                              Simulation_Context_type* p_Sim_Context) {

    /*
    
    Computes the magnetic field 4-vector, measured by a comoving obverver (with 4-velocity Plasma_velocity), called B_coord_frame,
    in the basis of a static observer (with 4-velocity n_mu = {1, 0, 0, 0} ). 
    Returns the magnitude of the magnetic field in the plasma frame.

    NOTE: The magnitude of the magnetic field in these frames is different, because its not concerved under Lorentz boosts.
    In the plasma frame I set the geometry of the field, then scale it by X_B_Plasma_norm_CGS.
    
    NOTE: B_coord_frame is measured by a static observer, so:
    1) B_coord_frame is effectively in the coordinate basis (not any ZAMO!) - this means one cannot directly apply Lorentz boosts to it.

    2) The magnetic field norm in here is considered to be in CGS

    */

    double Disk_density    = this->get_disk_density(State_Vector);
    double Hotspot_density = this->get_hotspot_density(State_Vector);

    double Disk_B_plasma_frame_norm_CGS    = sqrt(this->s_Disk_params.Magnetization * C_LIGHT_CGS * C_LIGHT_CGS * Disk_density * M_PROTON_CGS * 4 * M_PI);
    double Hotspot_B_plasma_frame_norm_CGS = sqrt(this->s_Hotspot_params.Magnetization * C_LIGHT_CGS * C_LIGHT_CGS * Hotspot_density * M_PROTON_CGS * 4 * M_PI); 

    double Disk_B_plasma_frame[4] = {                               0.0, 
                                     Disk_B_plasma_frame_norm_CGS * this->s_Disk_params.Mag_field_geometry[0],
                                     Disk_B_plasma_frame_norm_CGS * this->s_Disk_params.Mag_field_geometry[1],
                                     Disk_B_plasma_frame_norm_CGS * this->s_Disk_params.Mag_field_geometry[2] };

    double Hotspot_B_plasma_frame[4] = {                                  0.0,
                                        Hotspot_B_plasma_frame_norm_CGS * this->s_Hotspot_params.Mag_field_geometry[0],
                                        Hotspot_B_plasma_frame_norm_CGS * this->s_Hotspot_params.Mag_field_geometry[1],
                                        Hotspot_B_plasma_frame_norm_CGS * this->s_Hotspot_params.Mag_field_geometry[2] };

    /*
    
    Compute the the resultant magnetic field (the "background" disk field + the hotspot one) in the plasma frame.
    
    */

    double Interpolated_B_plasma_frame[4]{};
    double Interpolated_B_plasma_frame_norm_CGS{};

    for (int index = 1; index <= 3; index++) {

        /*
        
        Here I interpolate between the "background" disk field and the hotspot one, based on the electron density.
        
        */

        if (!isnan(Disk_density + Hotspot_density)) {

            Interpolated_B_plasma_frame[index] = (Disk_density * Disk_B_plasma_frame[index] + Hotspot_density * Hotspot_B_plasma_frame[index]) / (Disk_density + Hotspot_density);

        }
        else {

            Interpolated_B_plasma_frame[index] = Disk_B_plasma_frame[index];

        }

        Interpolated_B_plasma_frame_norm_CGS += Interpolated_B_plasma_frame[index] * Interpolated_B_plasma_frame[index];

    }

    Interpolated_B_plasma_frame_norm_CGS = sqrt(Interpolated_B_plasma_frame_norm_CGS);

    // Get the plasma velocity, so we can map the magnetic field back to the coordinate frame

    double* Plasma_velocity = this->get_disk_velocity(State_Vector, p_Sim_Context);

    /* 

    The only reference I could find for this is https://iopscience.iop.org/article/10.3847/1538-4357/ab718e/pdf 8e) and 8f).
    One arrives at the expression by projecting F^mu^nu (written in terms of B_coord_frame and Plasma_velocity),
    onto the 4-velocity of the static observer - this gives B_plasma_frame, then inverting the expression 
    to obtain B_coord_frame = f(B_plasma_frame, Plasma_Velocity)
    
    */

    Metric_type s_Metric = p_Sim_Context->p_Spacetime->get_metric(State_Vector);

    for (int left_idx = 0; left_idx <= 3; left_idx++) {

        for (int right_idx = 0; right_idx <= 3; right_idx++) {

            B_coord_frame[e_t_coord] += s_Metric.Metric[left_idx][right_idx] * Plasma_velocity[left_idx] * Interpolated_B_plasma_frame[right_idx];

        }

    }

    for (int index = 1; index <= 3; index++) {

        B_coord_frame[index] = ((Disk_B_plasma_frame[index] + Hotspot_B_plasma_frame[index]) + B_coord_frame[e_t_coord] * Plasma_velocity[index]) / Plasma_velocity[e_t_coord];
            
    }

    double Omega = Plasma_velocity[e_phi_coord] / Plasma_velocity[e_t_coord];
    double A = 1. / sqrt(-(s_Metric.Metric[e_t_coord][e_t_coord] + Omega * Omega * s_Metric.Metric[e_phi_coord][e_phi_coord]));

    B_coord_frame[e_t_coord] = A * sqrt(-s_Metric.Metric[e_phi_coord][e_phi_coord] / s_Metric.Metric[e_t_coord][e_t_coord]) * Omega;
    B_coord_frame[e_r_coord] = 0;
    B_coord_frame[e_theta_coord] = 0;
    B_coord_frame[e_phi_coord] = A / sqrt(-s_Metric.Metric[e_phi_coord][e_phi_coord] / s_Metric.Metric[e_t_coord][e_t_coord]);

    return Interpolated_B_plasma_frame_norm_CGS;

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

    double cos_angle = 1.0; 

    if (!isinf(1.0 / Wave_vec_dot_Plasma_vec) && !isinf(1.0 / B_field_norm_squared)) {

        cos_angle = Wave_vec_dot_B_field / Wave_vec_dot_Plasma_vec / sqrt(B_field_norm_squared);

    }

    if (fabs(cos_angle) < 1.0) {

        return cos_angle;

    }
    else {

        return cos_angle / fabs(cos_angle);

    }

}

/* =============================================== Thermal Synchotron Transfer Functions =============================================== */

void Generic_Optically_Thin_Model::evaluate_thermal_synchotron_transfer_functions(double Density,
                                                                                  double T_electron_dim,
                                                                                  double f_cyclo,
                                                                                  double sin_pitch_angle,
                                                                                  double cos_pitch_angle,
                                                                                  double Emission_functions[STOKES_PARAM_NUM],
                                                                                  double Faradey_functions[STOKES_PARAM_NUM],
                                                                                  Thermal_emission_f_arguments Emission_args,
                                                                                  Thermal_faradey_f_arguments Faradey_args) {

    for (int index = I; index <= STOKES_PARAM_NUM - 1; index++) {

        Emission_functions[index] = 0.0;
        Faradey_functions[index]  = 0.0;

    }

    /* ============ Extract the observational frequency (corrected with the redshift) from the emission arguments ============ */

    double& frequency = Emission_args.frequency;

    /* =================== These Bessel functions pop up as normalization factors in the expressions below =================== */

    double const K0_Bessel = std::cyl_bessel_k(0.0, 1.0 / T_electron_dim);
    double const K1_Bessel = std::cyl_bessel_k(1.0, 1.0 / T_electron_dim);
    double const K2_Bessel = std::cyl_bessel_k(2.0, 1.0 / T_electron_dim);
    double const f_crit = 3. / 2 * f_cyclo * T_electron_dim * T_electron_dim * sin_pitch_angle;

    double const omega_plasma_squared = 4 * M_PI * Density * Q_ELECTRON_CGS * Q_ELECTRON_CGS / M_ELECTRON_CGS;

    /* All the functions are normalized by a Bessel function, so I check if I can divide by it */
    if (!isinf(1e10 / K2_Bessel)) {

        this->get_thermal_synchotron_fit_functions(Emission_functions, Faradey_functions, &Emission_args, &Faradey_args);

        /* ================================================ The emission functions ================================================ */

        Emission_functions[I] *= Density * f_crit / K2_Bessel; // Scale the emission by the remaining position-dependant factors
        Emission_functions[Q] *= Density * f_crit / K2_Bessel; // Scale the emission by the remaining position-dependant factors


        if (!isinf(1e2 / sin_pitch_angle) && !isnan(1e2 / sin_pitch_angle)) {

            Emission_functions[V] *= Density * f_crit / K2_Bessel;                               // Scale the emission by the remaining position-dependant factors
            Emission_functions[V] *= (1. / T_electron_dim) * cos_pitch_angle / sin_pitch_angle;  // The V component has some extra angle dependance

        }

        /* ================================================ The faradey functions ================================================ */
        /* Originally derived in https://iopscience.iop.org/article/10.1086/592326/pdf - expressions 25, 26 and 33 */

        Faradey_functions[Q] *= omega_plasma_squared * (2 * M_PI * f_cyclo) * (2 * M_PI * f_cyclo) * sin_pitch_angle * sin_pitch_angle * (K1_Bessel / K2_Bessel + 6 * T_electron_dim);
        Faradey_functions[Q] /= 2 * C_LIGHT_CGS * (2 * M_PI * frequency) * (2 * M_PI * frequency) * (2 * M_PI * frequency);
        Faradey_functions[V] *= omega_plasma_squared * (2 * M_PI * f_cyclo) * cos_pitch_angle * (K0_Bessel) / K2_Bessel;
        Faradey_functions[V] /= C_LIGHT_CGS * (2 * M_PI * frequency) * (2 * M_PI * frequency);

    }

}

void Generic_Optically_Thin_Model::get_thermal_synchotron_transfer_functions(double State_vector[],
                                                                             Simulation_Context_type* p_Sim_Context,
                                                                             double Emission_functions[STOKES_PARAM_NUM],
                                                                             double Faradey_functions[STOKES_PARAM_NUM],
                                                                             double Absorbtion_functions[STOKES_PARAM_NUM],
                                                                             double Density,
                                                                             double Temperature,
                                                                             double* B_field,
                                                                             double B_field_norm) {

    /* Observation Frequency */
    double& const obs_frequency = p_Sim_Context->p_Init_Conditions->Observer_params.obs_frequency;

    /* Dimentionless Electron Temperature */
    double const T_electron_dim = BOLTZMANN_CONST_CGS * Temperature / M_ELECTRON_CGS / C_LIGHT_CGS / C_LIGHT_CGS;

    /* Redshit */
    double const redshift = Redshift(State_vector, this->get_disk_velocity(State_vector, p_Sim_Context), p_Sim_Context->p_Observer);

    /* Cyclotron Frequency */
    double const f_cyclo = Q_ELECTRON_CGS * B_field_norm / (2 * M_PI * M_ELECTRON_CGS * C_LIGHT_CGS);

    /* The "averaged" critical frequency (without the sin(theta) term - that gets added on later from a pre-computed table) */
    double const f_crit_no_sin = 3. / 2 * f_cyclo * T_electron_dim * T_electron_dim;

    Thermal_emission_f_arguments Emission_args_ang_uncorrected{};
    Thermal_faradey_f_arguments Faradey_args_ang_uncorrected{};

    /* Both the emission and faradey function expressions are in terms of an dimentionless variable X, but the definitions for X are different */
    Emission_args_ang_uncorrected = {1e100,  // X
                                     1e100,  // sqrt_X
                                     1e100,  // cbrt_X
                                     1e100}; // frequency

    Faradey_args_ang_uncorrected = { 1e100,   // X
                                     1e100,   // X_to_1_point_035
                                     1e100,   // X_to_1_point_2
                                     1e100 }; // frequency

    /* Compute all the wierd powers of X outside the pitch angle averaging loop */

    if (f_crit_no_sin > std::numeric_limits<double>::min()) {

        Emission_args_ang_uncorrected.X      = obs_frequency / f_crit_no_sin / redshift;
        Emission_args_ang_uncorrected.sqrt_X = sqrt(Emission_args_ang_uncorrected.X);
        Emission_args_ang_uncorrected.cbrt_X = cbrt(Emission_args_ang_uncorrected.X);

        Faradey_args_ang_uncorrected.X                = T_electron_dim * sqrt(M_SQRT2 * 1e3 * f_cyclo / (obs_frequency / redshift));
        Faradey_args_ang_uncorrected.X_to_1_point_2   = pow(Faradey_args_ang_uncorrected.X, 1.2f);
        Faradey_args_ang_uncorrected.X_to_1_point_035 = pow(Faradey_args_ang_uncorrected.X, 1.035f);

    }

    if (p_Sim_Context->p_Init_Conditions->Average_electron_pitch_angle) {

        /* ============ This loop averages over the emission pitch angle, which it gets from a pre-computed table ============ */

        for (int averaging_idx = 0; averaging_idx <= NUM_SAMPLES_TO_AVG - 1; averaging_idx++) {

            Thermal_emission_f_arguments Emission_args_ang_corrected{};
            Thermal_faradey_f_arguments Faradey_args_ang_corrected{};

            double& sin_pitch_angle = this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[averaging_idx];
            double& cos_pitch_angle = this->s_Precomputed_e_pitch_angles.cos_electron_pitch_angles[averaging_idx];

            Emission_args_ang_corrected.X         = Emission_args_ang_uncorrected.X / sin_pitch_angle;
            Emission_args_ang_corrected.sqrt_X    = Emission_args_ang_uncorrected.sqrt_X * this->s_Precomputed_e_pitch_angles.one_over_sqrt_sin[averaging_idx];
            Emission_args_ang_corrected.cbrt_X    = Emission_args_ang_uncorrected.cbrt_X * this->s_Precomputed_e_pitch_angles.one_over_cbrt_sin[averaging_idx];
            Emission_args_ang_corrected.frequency = obs_frequency / redshift;

            Faradey_args_ang_corrected.X                = Faradey_args_ang_uncorrected.X / this->s_Precomputed_e_pitch_angles.one_over_sqrt_sin[averaging_idx];
            Faradey_args_ang_corrected.X_to_1_point_035 = Faradey_args_ang_uncorrected.X_to_1_point_035 / this->s_Precomputed_e_pitch_angles.one_over_sin_to_1_point_035[averaging_idx];
            Faradey_args_ang_corrected.X_to_1_point_2   = Faradey_args_ang_uncorrected.X_to_1_point_2 / this->s_Precomputed_e_pitch_angles.one_over_sin_to_1_point_2_over_2[averaging_idx];
            Faradey_args_ang_corrected.frequency        = obs_frequency / redshift;

            double temp_emission_functions[STOKES_PARAM_NUM]{};
            double temp_faradey_functions[STOKES_PARAM_NUM]{};

            this->evaluate_thermal_synchotron_transfer_functions(Density, T_electron_dim, f_cyclo, sin_pitch_angle, cos_pitch_angle, temp_emission_functions, temp_faradey_functions, Emission_args_ang_corrected, Faradey_args_ang_corrected);

            // The U component is 0 by definition
            Emission_functions[I] += temp_emission_functions[I] * sin_pitch_angle * M_PI / NUM_SAMPLES_TO_AVG / 2;
            Emission_functions[Q] += temp_emission_functions[Q] * sin_pitch_angle * M_PI / NUM_SAMPLES_TO_AVG / 2;
            Emission_functions[U] = 0.0;
            Emission_functions[V] += temp_emission_functions[V] * sin_pitch_angle * M_PI / NUM_SAMPLES_TO_AVG / 2;

            // The I and U components are 0 by definition
            Faradey_functions[I] = 0.0f;
            Faradey_functions[Q] += temp_faradey_functions[Q] * sin_pitch_angle * M_PI / NUM_SAMPLES_TO_AVG / 2;
            Faradey_functions[U] = 0.0f;
            Faradey_functions[V] += temp_faradey_functions[V] * sin_pitch_angle * M_PI / NUM_SAMPLES_TO_AVG / 2;

        }
    }
    else {

        /* The magnetic field is the one measured by a comoving with the plasma observer, but expressed in the cooridante frame */
        
        double cos_pitch_angle = get_electron_pitch_angle(B_field, State_vector, p_Sim_Context);
        double sin_pitch_angle = sqrt(1.0 - cos_pitch_angle * cos_pitch_angle);

        double one_over_sqrt_sin = 1.0 / sqrt(sin_pitch_angle);
        double one_over_cbrt_sin = 1.0 / cbrt(sin_pitch_angle);

        Thermal_emission_f_arguments Emission_args_ang_corrected{};
        Thermal_faradey_f_arguments Faradey_args_ang_corrected{};

        Emission_args_ang_corrected.X         = Emission_args_ang_uncorrected.X / sin_pitch_angle;
        Emission_args_ang_corrected.sqrt_X    = Emission_args_ang_uncorrected.sqrt_X * one_over_sqrt_sin;
        Emission_args_ang_corrected.cbrt_X    = Emission_args_ang_uncorrected.cbrt_X * one_over_cbrt_sin;
        Emission_args_ang_corrected.frequency = obs_frequency / redshift;

        Faradey_args_ang_corrected.X                = Faradey_args_ang_uncorrected.X / one_over_sqrt_sin;
        Faradey_args_ang_corrected.X_to_1_point_035 = Faradey_args_ang_uncorrected.X_to_1_point_035 * pow(sin_pitch_angle, 1.035);
        Faradey_args_ang_corrected.X_to_1_point_2   = Faradey_args_ang_uncorrected.X_to_1_point_2 * pow(sin_pitch_angle, 1.2);
        Faradey_args_ang_corrected.frequency        = obs_frequency / redshift;

        this->evaluate_thermal_synchotron_transfer_functions(Density, T_electron_dim, f_cyclo, sin_pitch_angle, cos_pitch_angle, Emission_functions, Faradey_functions, Emission_args_ang_corrected, Faradey_args_ang_corrected);

    }

    /* ================================================ The absorbtion functions ================================================ */

    double Planck_function_CGS = get_planck_function_CGS(obs_frequency / redshift, this->get_disk_temperature(State_vector));

    if (Planck_function_CGS > std::numeric_limits<double>::min()) {

        for (int stokes_index = 0; stokes_index <= STOKES_PARAM_NUM - 1; stokes_index++) {

            Absorbtion_functions[stokes_index] = Emission_functions[stokes_index] / Planck_function_CGS;

        }
    }
    else {

        for (int stokes_index = 0; stokes_index <= STOKES_PARAM_NUM - 1; stokes_index++) {

            Absorbtion_functions[stokes_index] = 0.0;

        }

    }
}

/* ========================================== Kappa Synchotron Transfer Functions ========================================== */

void Generic_Optically_Thin_Model::evaluate_kappa_synchotron_transfer_functions(double Density,
                                                                                double f_cyclo,
                                                                                double Emission_functions[STOKES_PARAM_NUM],
                                                                                double Faradey_functions[STOKES_PARAM_NUM],
                                                                                double Absorbtion_functions[STOKES_PARAM_NUM],
                                                                                Kappa_transfer_f_arguments Transfer_args) {


    /* ===================== Extract the observation frequency (corrected with the redshift) from the transfer args ===================== */

    double frequency = Transfer_args.X * f_cyclo * (Transfer_args.T_electron_dim * Transfer_args.kappa) * 
                                                   (Transfer_args.T_electron_dim * Transfer_args.kappa) * 
                                                    Transfer_args.sin_emission_angle;

    for (int index = 0; index <= STOKES_PARAM_NUM - 1; index++) {

        Emission_functions[index] = 0.0;
        Faradey_functions[index] = 0.0;
        Absorbtion_functions[index] = 0.0;

    }

    this->get_kappa_synchotron_fit_functions(Emission_functions, Faradey_functions, Absorbtion_functions, &Transfer_args);

    Emission_functions[I] *= Density * Q_ELECTRON_CGS * Q_ELECTRON_CGS / C_LIGHT_CGS * f_cyclo; 
    Emission_functions[Q] *= Density * Q_ELECTRON_CGS * Q_ELECTRON_CGS / C_LIGHT_CGS * f_cyclo; 
    Emission_functions[V] *= Density * Q_ELECTRON_CGS * Q_ELECTRON_CGS / C_LIGHT_CGS * f_cyclo; 

    Absorbtion_functions[I] *= Density * Q_ELECTRON_CGS * Q_ELECTRON_CGS / M_ELECTRON_CGS / C_LIGHT_CGS / frequency;
    Absorbtion_functions[Q] *= Density * Q_ELECTRON_CGS * Q_ELECTRON_CGS / M_ELECTRON_CGS / C_LIGHT_CGS / frequency;
    Absorbtion_functions[V] *= Density * Q_ELECTRON_CGS * Q_ELECTRON_CGS / M_ELECTRON_CGS / C_LIGHT_CGS / frequency;

}


void Generic_Optically_Thin_Model::get_kappa_synchotron_transfer_functions(double State_Vector[],
                                                                           Simulation_Context_type* p_Sim_Context,
                                                                           double Emission_functions[STOKES_PARAM_NUM],
                                                                           double Faradey_functions[STOKES_PARAM_NUM],
                                                                           double Absorbtion_functions[STOKES_PARAM_NUM],
                                                                           double Density,
                                                                           double Temperature,
                                                                           double* B_field,
                                                                           double B_field_norm){

    /* === The transfer functions arrays needs to be manually cleared, because this function only adds to it. === */
    for (int stokes_index = 0; stokes_index <= STOKES_PARAM_NUM - 1; stokes_index++) { 

        Emission_functions[stokes_index] = 0.0; 
        Faradey_functions[stokes_index] = 0.0; 
        Absorbtion_functions[stokes_index] = 0.0;  

    }

    /* Observation frequency */

    double& obs_frequency = p_Sim_Context->p_Init_Conditions->Observer_params.obs_frequency;

    /* Dimentionless Electron Temperature */
    double T_electron_dim = BOLTZMANN_CONST_CGS * Temperature / M_ELECTRON_CGS / C_LIGHT_CGS / C_LIGHT_CGS;

    /* Redshit */
    double redshift = Redshift(State_Vector, this->get_disk_velocity(State_Vector, p_Sim_Context), p_Sim_Context->p_Observer);

    /* Cyclotron Frequency */
    double f_cyclo = Q_ELECTRON_CGS * B_field_norm / (2 * M_PI * M_ELECTRON_CGS * C_LIGHT_CGS);

    /* The "averaged" critical frequency (without the sin(theta) term - that gets added on later from a pre-computed table) */
    double f_k_no_sin = f_cyclo * (this->s_Emission_params.Kappa * T_electron_dim) * (this->s_Emission_params.Kappa * T_electron_dim);

    Kappa_transfer_f_arguments Transfer_args_uncorrected{};

    /* Init the argument X to something really large, which would correspond to no emission */
    Transfer_args_uncorrected = { 1e100, // X
                                  1e100, // sqrt_X
                                  1e100, // cbrt_X
                                  1e100, // X_to_7_over_20
                                  this->s_Emission_params.Kappa,
                                  0.0,
                                  T_electron_dim };

    if (f_k_no_sin > std::numeric_limits<double>::min()) {

        Transfer_args_uncorrected.X              = obs_frequency / f_k_no_sin / redshift;
        Transfer_args_uncorrected.sqrt_X         = sqrt(Transfer_args_uncorrected.X);
        Transfer_args_uncorrected.cbrt_X         = cbrt(Transfer_args_uncorrected.X);
        Transfer_args_uncorrected.X_to_7_over_20 = pow(Transfer_args_uncorrected.X, 7. / 20);

    }

    if (p_Sim_Context->p_Init_Conditions->Average_electron_pitch_angle) {

        /* ============ This loop averages over the emission pitch angle, which it gets from a pre-computed table ============ */

        for (int averaging_idx = 0; averaging_idx <= NUM_SAMPLES_TO_AVG - 1; averaging_idx++) {

            double& sin_pitch_angle = this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[averaging_idx];
            double& cos_pitch_angle = this->s_Precomputed_e_pitch_angles.cos_electron_pitch_angles[averaging_idx];

            Kappa_transfer_f_arguments Transfer_args_ang_corrected = { 1e100, // X
                                                                       1e100, // sqrt_X
                                                                       1e100, // cbrt_X
                                                                       1e100, // X_to_7_over_20
                                                                       this->s_Emission_params.Kappa,
                                                                       sin_pitch_angle,
                                                                       T_electron_dim };

            Transfer_args_ang_corrected.X = Transfer_args_uncorrected.X / sin_pitch_angle;
            Transfer_args_ang_corrected.sqrt_X = Transfer_args_uncorrected.sqrt_X * this->s_Precomputed_e_pitch_angles.one_over_sqrt_sin[averaging_idx];
            Transfer_args_ang_corrected.cbrt_X = Transfer_args_uncorrected.cbrt_X * this->s_Precomputed_e_pitch_angles.one_over_cbrt_sin[averaging_idx];
            Transfer_args_ang_corrected.X_to_7_over_20 = Transfer_args_uncorrected.X_to_7_over_20 * this->s_Precomputed_e_pitch_angles.one_over_sin_to_7_over_20[averaging_idx];

            double temp_emission_functions[STOKES_PARAM_NUM]{};
            double temp_faradey_functions[STOKES_PARAM_NUM]{};
            double temp_absorbtion_functions[STOKES_PARAM_NUM]{};

            this->evaluate_kappa_synchotron_transfer_functions(Density, f_cyclo, temp_emission_functions, temp_faradey_functions, temp_absorbtion_functions, Transfer_args_ang_corrected);

            // The U component is 0 by definition
            Emission_functions[I] += temp_emission_functions[I] * sin_pitch_angle * M_PI / NUM_SAMPLES_TO_AVG / 2;
            Emission_functions[Q] += temp_emission_functions[Q] * sin_pitch_angle * M_PI / NUM_SAMPLES_TO_AVG / 2;
            Emission_functions[U] = 0.0;
            Emission_functions[V] += temp_emission_functions[V] * sin_pitch_angle * M_PI / NUM_SAMPLES_TO_AVG / 2;

            // The I and U components are 0 by definition
            Faradey_functions[I] = 0.0f;
            Faradey_functions[Q] += temp_faradey_functions[Q] * sin_pitch_angle * M_PI / NUM_SAMPLES_TO_AVG / 2;
            Faradey_functions[U] = 0.0f;
            Faradey_functions[V] += temp_faradey_functions[V] * sin_pitch_angle * M_PI / NUM_SAMPLES_TO_AVG / 2;

        }
    }
    else {

        /* The magnetic field is the one measured by a comoving with the plasma observer, but expressed in the cooridante frame */

        double cos_pitch_angle = get_electron_pitch_angle(B_field, State_Vector, p_Sim_Context);
        double sin_pitch_angle = sqrt(1. - cos_pitch_angle * cos_pitch_angle);

        double one_over_sqrt_sin    = 1. / sqrt(sin_pitch_angle);
        double one_over_cbrt_sin    = 1. / cbrt(sin_pitch_angle);
        double one_over_7_to_20_sin = 1. / pow(sin_pitch_angle, 7. / 20);

        Kappa_transfer_f_arguments Transfer_args_ang_corrected = { 1e100, // X
                                                                   1e100, // sqrt_X
                                                                   1e100, // cbrt_X
                                                                   1e100, // X_to_7_over_20
                                                                   this->s_Emission_params.Kappa,
                                                                   sin_pitch_angle,
                                                                   T_electron_dim };

        Transfer_args_ang_corrected.X = Transfer_args_uncorrected.X / sin_pitch_angle;
        Transfer_args_ang_corrected.sqrt_X = Transfer_args_uncorrected.sqrt_X * one_over_sqrt_sin;
        Transfer_args_ang_corrected.cbrt_X = Transfer_args_uncorrected.cbrt_X * one_over_cbrt_sin;
        Transfer_args_ang_corrected.X_to_7_over_20 = Transfer_args_uncorrected.X_to_7_over_20 * one_over_7_to_20_sin;

        this->evaluate_kappa_synchotron_transfer_functions(Density, f_cyclo, Emission_functions, Faradey_functions, Absorbtion_functions, Transfer_args_ang_corrected);

    }

}

/* ========================================== Phenomenological Synchotron Transfer Functions ========================================== */

void Generic_Optically_Thin_Model::get_phenomenological_synchotron_functions(double State_Vector[],
                                                                             Simulation_Context_type* p_Sim_Context, 
                                                                             double Emission_functions[STOKES_PARAM_NUM],
                                                                             double Faradey_functions[STOKES_PARAM_NUM],
                                                                             double Absorbtion_functions[STOKES_PARAM_NUM],
                                                                             double Density) {

    double& emission_power_law = this->s_Emission_params.Phenomenological_emission_power_law;
    double& source_f_power_law = this->s_Emission_params.Phenomenological_source_f_power_law;
    double& emission_coeff     = this->s_Emission_params.Phenomenological_emission_coeff;
    double& abs_coeff          = this->s_Emission_params.Phenomenological_absorbtion_coeff;

    double* U_source_coord = this->get_disk_velocity(State_Vector, p_Sim_Context);
    double redshift        = Redshift(State_Vector, U_source_coord, p_Sim_Context->p_Observer);

    Emission_functions[I] = emission_coeff * Density / this->s_Disk_params.Electron_density_scale * pow(redshift, emission_power_law);
    Emission_functions[Q] = 0.0;
    Emission_functions[U] = 0.0;
    Emission_functions[V] = 0.0;

    Absorbtion_functions[I] = abs_coeff * emission_coeff * Density / this->s_Disk_params.Electron_density_scale * pow(redshift, source_f_power_law + emission_power_law);
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

    for (int index = I; index <= STOKES_PARAM_NUM - 1; index++) {

        Emission_functions[index] = 0.0;
        Faradey_functions[index] = 0.0;
        Absorbtion_functions[index] = 0.0;

    }

    /* Disk Electron Density in CGS */
    double Disk_density = this->get_disk_density(State_Vector);

    /* Disk Electron Temperature */
    double Disk_temp = this->get_disk_temperature(State_Vector);

    /* Hotspot Electron Density in CGS */
    double Hotspot_density = this->get_hotspot_density(State_Vector);

    /* Hotspot Electron Temperature */
    double Hotspot_temp = this->get_hotspot_temperature(State_Vector);

    double Interpolation_param = Disk_density / (Disk_density + Hotspot_density);


    /* 
    
    Total Magnetic Field (hotspot + "background" disk), measured by a comoving with the plasma observer, projected back onto a static observer (a.e. in the coordinate basis)
    
    */

    double B_field_local[4]{};
    double B_norm_CGS = this->get_total_magnetic_field(B_field_local, State_Vector, p_Sim_Context);

    double Disk_emission_functions[STOKES_PARAM_NUM]{};
    double Disk_faradey_functions[STOKES_PARAM_NUM]{};
    double Disk_absorbtion_functions[STOKES_PARAM_NUM]{};

    if (0.0 != p_Sim_Context->p_Init_Conditions->Disk_params.Electron_density_scale) {

        switch (p_Sim_Context->p_Init_Conditions->Disk_params.Ensamble_type) {

        case(e_Phenomenological_ensamble):

            this->get_phenomenological_synchotron_functions(State_Vector, p_Sim_Context, Disk_emission_functions, Disk_faradey_functions, Disk_absorbtion_functions, Interpolation_param * Disk_density);
            break;

        case(e_Kappa_ensamble):

            this->get_kappa_synchotron_transfer_functions(State_Vector, p_Sim_Context, Disk_emission_functions, Disk_faradey_functions, Disk_absorbtion_functions,
                Interpolation_param * Disk_density, Interpolation_param * Disk_temp, B_field_local, B_norm_CGS);
            break;

        default:

            this->get_thermal_synchotron_transfer_functions(State_Vector, p_Sim_Context, Disk_emission_functions, Disk_faradey_functions, Disk_absorbtion_functions,
                Interpolation_param * Disk_density, Interpolation_param * Disk_temp, B_field_local, B_norm_CGS);
            break;
        }
    }

    double Hotspot_emission_functions[STOKES_PARAM_NUM]{};
    double Hotspot_faradey_functions[STOKES_PARAM_NUM]{};
    double Hotspot_absorbtion_functions[STOKES_PARAM_NUM]{};

    if (0.0 != p_Sim_Context->p_GOT_Model->s_Hotspot_params.Electron_density_scale) {

        switch (p_Sim_Context->p_Init_Conditions->Hotspot_params.Ensamble_type) {

        case(e_Phenomenological_ensamble):

            this->get_phenomenological_synchotron_functions(State_Vector, p_Sim_Context, Hotspot_emission_functions, Hotspot_faradey_functions, Hotspot_absorbtion_functions, (1. - Interpolation_param) * Hotspot_density);
            break;

        case(e_Kappa_ensamble):

            this->get_kappa_synchotron_transfer_functions(State_Vector, p_Sim_Context, Hotspot_emission_functions, Hotspot_faradey_functions, Hotspot_absorbtion_functions,
                                                          (1. - Interpolation_param) * Hotspot_density, (1. - Interpolation_param) * Hotspot_temp, B_field_local, B_norm_CGS);
            break;

        default:

            this->get_thermal_synchotron_transfer_functions(State_Vector, p_Sim_Context, Hotspot_emission_functions, Hotspot_faradey_functions, Hotspot_absorbtion_functions,
                                                            (1. - Interpolation_param) * Hotspot_density, (1. - Interpolation_param) * Hotspot_temp, B_field_local, B_norm_CGS);
            break;
        }

    }

    for (int index = I; index <= STOKES_PARAM_NUM - 1; index++) {

        Emission_functions[index]   = Disk_emission_functions[index] + Hotspot_emission_functions[index];
        Faradey_functions[index]    = Disk_faradey_functions[index] + Hotspot_faradey_functions[index];
        Absorbtion_functions[index] = Disk_absorbtion_functions[index] + Hotspot_absorbtion_functions[index];

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

            this->s_Precomputed_e_pitch_angles.one_over_sin_to_1_point_035[index] = 1. / pow(this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index], 1.035);
            this->s_Precomputed_e_pitch_angles.one_over_sin_to_1_point_2_over_2[index] = 1. / pow(this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index], 1. / 2);

            // Used in the kappa synchotron emission functions

            this->s_Precomputed_e_pitch_angles.one_over_sin_to_7_over_20[index] = 1. / pow(this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index], 7. / 20);

        }
    }
}

int Generic_Optically_Thin_Model::load_parameters(Disk_model_parameters_type* p_Disk_params, Hotspot_model_parameters_type* p_Hotspot_params, Emission_model_parameters_type* p_Emission_params) {

    if (NULL != p_Disk_params) {

        this->s_Disk_params = *p_Disk_params;

    }
    else {

        std::cout << "p_Disk_params is a NULL pointer! \n";

        exit(ERROR);
    }

    if (NULL != p_Hotspot_params) {

        this->s_Hotspot_params = *p_Hotspot_params;

    }
    else {

        std::cout << "p_Hotspot_params is a NULL pointer! \n";

        exit(ERROR);
    }

    if (NULL != p_Emission_params) {

        this->s_Emission_params = *p_Emission_params;

    }
    else {

        std::cout << "p_Emission_params is a NULL pointer! \n";

        exit(ERROR);
    }

    return OK;

}

double get_planck_function_CGS(double Frequency, double Temperature) {

    return 2 * PLANCK_CONSTANT_CGS * Frequency * Frequency * Frequency / C_LIGHT_CGS / C_LIGHT_CGS / (exp(PLANCK_CONSTANT_CGS * Frequency / BOLTZMANN_CONST_CGS / Temperature) - 1.);

}
