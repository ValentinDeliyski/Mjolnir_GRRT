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

//! Copies over initial conditions from the Simulation Context struct to internal class variables for the sake of convenicence.
/*! Copies over initial conditions from the Simulation Context struct to internal class variables for the sake of convenicence.
 * 
 *   \param [in] p_Sim_Context - Pointer to the Simulation Context struct.
 *   \return Nothing.
 */
Novikov_Thorne_Model::Novikov_Thorne_Model(Simulation_Context_type* p_Sim_Context) {

    this->r_in  = p_Sim_Context->p_Init_Conditions->NT_params.r_in;
    this->r_out = p_Sim_Context->p_Init_Conditions->NT_params.r_out;
    this->flux_integral_accuracy = p_Sim_Context->p_Init_Conditions->Integrator_params.Simpson_accuracy;
    this->p_Spacetime = p_Sim_Context->p_Spacetime;
    this->e_Spacetime = p_Sim_Context->p_Init_Conditions->Metric_params.e_Spacetime;

};

//! Evaluates the Keplarian angular velocity of the Novikov-Thorne disk model.
/*! Evaluates the Keplarian angular velocity of the Novikov-Thorne disk model.
 *
 *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
 *   \return The Keplarian angular velocity.
 */
double Novikov_Thorne_Model::Keplerian_angular_velocity(const double* const State_Vector) {

    Metric_type s_dr_Metric = this->p_Spacetime->get_dr_metric(State_Vector);

    return (-s_dr_Metric.Metric[0][3] + sqrt(s_dr_Metric.Metric[0][3] * s_dr_Metric.Metric[0][3] - s_dr_Metric.Metric[0][0] * s_dr_Metric.Metric[3][3])) / s_dr_Metric.Metric[3][3];

}

//! Evaluates the radial derivative of the Keplarian angular velocity of the Novikov-Thorne disk model.
/*! Evaluates the radial derivative of the Keplarian angular velocity of the Novikov-Thorne disk model.
 *
 *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
 *   \return The radial derivative of the Keplarian angular velocity.
 */
double Novikov_Thorne_Model::dr_Keplerian_angular_velocity(const double* const State_Vector) {

    Metric_type s_dr_Metric = this->p_Spacetime->get_dr_metric(State_Vector);
    Metric_type s_d2r_Metric = this->p_Spacetime->get_d2r_metric(State_Vector);

    double root = sqrt(s_dr_Metric.Metric[0][3] * s_dr_Metric.Metric[0][3] - s_dr_Metric.Metric[0][0] * s_dr_Metric.Metric[3][3]);

    double Kepler = this->Keplerian_angular_velocity(State_Vector);

    return  - Kepler / s_dr_Metric.Metric[3][3] * s_d2r_Metric.Metric[3][3] + (-s_d2r_Metric.Metric[0][3]
            + 1.0 / root / 2 * (2 * s_dr_Metric.Metric[0][3] * s_d2r_Metric.Metric[0][3] - s_dr_Metric.Metric[0][0] * s_d2r_Metric.Metric[3][3]
            - s_d2r_Metric.Metric[0][0] * s_dr_Metric.Metric[3][3])) / s_dr_Metric.Metric[3][3];

}

//! Evaluates the redshift of the Novikov-Thorne disk model.
/*! Evaluates the redshift of the Novikov-Thorne disk model.
 *
 *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
 *   \return The redshift.
 */
double Novikov_Thorne_Model::Redshift(const double* const State_Vector, double r_obs, double theta_obs) {

    const double& r_source = State_Vector[e_r];
    const double& theta_source = State_Vector[e_theta];

    /*
    Get the observer 4-velocity
    */

    double State_Vector_obs[4] = {0, r_obs, theta_obs, 0 };

    Metric_type s_Metric_obs = this->p_Spacetime->get_metric(State_Vector_obs);

    double U_obs[4] = { 1.0 / s_Metric_obs.Lapse_function, 0 ,0 , s_Metric_obs.Shift_function / s_Metric_obs.Lapse_function };

    /*
    Get the source 4-velocity
    */

    Metric_type s_Metric_source = this->p_Spacetime->get_metric(State_Vector);

    double Kepler = this->Keplerian_angular_velocity(State_Vector);

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

    return  (-U_obs[0] + U_obs[3] * State_Vector[e_p_phi]) / (-U_source[0] + U_source[3] * State_Vector[e_p_phi]);

}

//! Evaluates the energy of the Novikov-Thorne disk model.
/*! Evaluates the energy of the Novikov-Thorne disk model.
 *
 *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
 *   \return The energy.
 */
double Novikov_Thorne_Model::disk_Energy(const double* const State_Vector) {

    Metric_type s_Metric_source = this->p_Spacetime->get_metric(State_Vector);

    double Kepler = this->Keplerian_angular_velocity(State_Vector);

    double root = sqrt(-s_Metric_source.Metric[0][0] - 2 * s_Metric_source.Metric[0][3] * Kepler - s_Metric_source.Metric[3][3] * Kepler * Kepler);

    return  -(s_Metric_source.Metric[0][0] + s_Metric_source.Metric[0][3] * Kepler) / root;

}

//! Evaluates the angular momentum magnitude of the Novikov-Thorne disk model.
/*! Evaluates the angular momentum magnitude of the Novikov-Thorne disk model.
 *
 *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
 *   \return The magnitude of the angular momentum.
 */
double Novikov_Thorne_Model::disk_Angular_Momentum(const double* const State_Vector) {

    Metric_type s_Metric_source = this->p_Spacetime->get_metric(State_Vector);

    double Kepler = this->Keplerian_angular_velocity(State_Vector);

    double root = sqrt(-s_Metric_source.Metric[0][0] - 2 * s_Metric_source.Metric[0][3] * Kepler - s_Metric_source.Metric[3][3] * Kepler * Kepler);

    return  (s_Metric_source.Metric[3][3] * Kepler + s_Metric_source.Metric[0][3]) / root;

}

//! Evaluates the integrand of the integral that appears in the flux expression of the Novikov-Thorne disk model.
/*! Evaluates the integrand of the integral that appears in the flux expression of the Novikov-Thorne disk model.
 *
 *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
 *   \return The value of the integrand.
 */
double Novikov_Thorne_Model::Flux_integrand(const double* const State_Vector) {

    Metric_type s_Metric = this->p_Spacetime->get_metric(State_Vector);
    Metric_type s_dr_Metric = this->p_Spacetime->get_dr_metric(State_Vector);

    double Kepler = this->Keplerian_angular_velocity(State_Vector);
    double dr_Kepler = this->dr_Keplerian_angular_velocity(State_Vector);

    double root = sqrt(-s_Metric.Metric[0][0] - 2 * s_Metric.Metric[0][3] * Kepler - s_Metric.Metric[3][3] * Kepler * Kepler);
    double dr_root = (-s_dr_Metric.Metric[0][0] - 2 * (s_dr_Metric.Metric[0][3] * Kepler + s_Metric.Metric[0][3] * dr_Kepler)
        - s_dr_Metric.Metric[3][3] * Kepler * Kepler - 2 * s_Metric.Metric[3][3] * Kepler * dr_Kepler);

    double E = this->disk_Energy(State_Vector);
    double L = this->disk_Angular_Momentum(State_Vector);

    double dr_L = (s_dr_Metric.Metric[3][3] * Kepler + s_Metric.Metric[3][3] * dr_Kepler + s_dr_Metric.Metric[0][3]) / root - L / root / root / 2 * dr_root;

    return (E - Kepler * L) * dr_L;

}

//! Evaluates the integral that appears in the flux expression of the Novikov-Thorne disk model.
/*! Evaluates the integral that appears in the flux expression of the Novikov-Thorne disk model, using the adaptive Simpson method.
 *
 *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
 *   \return The value of the integral term.
 */
double Novikov_Thorne_Model::solve_Flux_integral(double r_in, const double* const State_Vector, double tolerance) {

    const double& lower_bound = r_in;
    const double& upper_bound = State_Vector[e_r];

    double mid_point          = (lower_bound + upper_bound) / 2;
    double left_of_mid_point  = (lower_bound + mid_point) / 2;
    double right_of_mid_point = (mid_point + upper_bound) / 2;

    double lower_bound_state_vector[4]{};
    double mid_point_state_vector[4]{};
    double left_of_mid_point_state_vector[4]{};
    double right_of_mid_point_state_vector[4]{};

    memcpy(lower_bound_state_vector, State_Vector, 4 * sizeof(double));
    memcpy(mid_point_state_vector, State_Vector, 4 * sizeof(double));
    memcpy(left_of_mid_point_state_vector, State_Vector, 4 * sizeof(double));
    memcpy(right_of_mid_point_state_vector, State_Vector, 4 * sizeof(double));

    lower_bound_state_vector[e_r] = r_in;
    mid_point_state_vector[e_r] = mid_point;
    left_of_mid_point_state_vector[e_r] = left_of_mid_point;
    right_of_mid_point_state_vector[e_r] = right_of_mid_point;

    double F_lower_bound = this->Flux_integrand(lower_bound_state_vector);
    double F_mid_point   = this->Flux_integrand(mid_point_state_vector);
    double F_upper_bound = this->Flux_integrand(State_Vector);

    double F_left_mid = this->Flux_integrand(left_of_mid_point_state_vector);
    double F_right_mid = this->Flux_integrand(right_of_mid_point_state_vector);

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

        double L_value = this->solve_Flux_integral(lower_bound, mid_point_state_vector, tolerance / 2);
        double R_value = this->solve_Flux_integral(mid_point, State_Vector, tolerance / 2);

        integral = L_value + R_value;

    }

    return integral;
}

//! Evaluates the flux of the Novikov-Thorne disk model
/*! Evaluates the flux of the Novikov-Thorne disk model
 * 
 *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
 *   \return The Novikov-Thorne flux in [M_dot / M^2].
 */
double Novikov_Thorne_Model::get_flux(const double* const State_Vector) {

    Metric_type s_Metric = this->p_Spacetime->get_metric(State_Vector);

    double metric_det = get_eq_induced_metric_det(s_Metric.Metric);
    double E_disk = disk_Energy(State_Vector);
    double L_disk = disk_Angular_Momentum(State_Vector);

    double Kepler = Keplerian_angular_velocity(State_Vector);
    double dr_Kepler = dr_Keplerian_angular_velocity(State_Vector);

    double Flux_coeff = -dr_Kepler / ((E_disk - Kepler * L_disk) * (E_disk - Kepler * L_disk)) / (4 * M_PI * sqrt(-metric_det));

    double Flux_integral = solve_Flux_integral(this->r_in, State_Vector, this->flux_integral_accuracy);

    return Flux_coeff * Flux_integral;

}

/***********************************************************
|                                                          |
| General Optically Thin Model Class Functions Definitions |
|                                                          |
***********************************************************/

/* ==================================================== Temperature Functions ===================================================== */

//! Computes the background accretion disk temperature
/*! Computes the background accretion disk at temperature the current photon position.
 *
 *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
 *   \return The temperature in [K].
 */
double Generic_Optically_Thin_Model::get_disk_temperature(const double* const State_Vector) {


    const double& r            = State_Vector[e_r];
    const double& T_scale      = this->s_Disk_params.Electron_temperature_scale;
    const double& R_0          = this->s_Disk_params.Power_law_temperature_R_0;
    const double& r_cutoff     = this->s_Disk_params.Power_law_temperature_R_cutoff;
    const double& Cutoff_scale = this->s_Disk_params.Power_law_temperature_cutoff_scale;
    const double& Power_law    = this->s_Disk_params.Power_law_temperature_radial_power_law;

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

    double Disk_Temperature = T_scale * Disk_temperature_profile;

    return Disk_Temperature;

}

//! Computes the hotspot temperature
/*! Computes the hotspot at temperature the current photon position.
 *
 *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
 *   \return The temperature in [K].
 */
double Generic_Optically_Thin_Model::get_hotspot_temperature(const double* const State_Vector) {

    const double& Hotspot_r      = this->s_Hotspot_params.Position[e_r - 1];
    const double& Hotspot_theta  = this->s_Hotspot_params.Position[e_theta - 1];
    const double& Hotspot_phi    = this->s_Hotspot_params.Position[e_phi - 1];
    const double& Spatial_spread = this->s_Hotspot_params.Temperature_spread;

    const double& photon_r = State_Vector[e_r];

    // I these this more than once, so I precompute them
    double sin_photon_theta  = sin(State_Vector[e_theta]);
    double sin_hotspot_theta = sin(Hotspot_theta);

    double x_center = Hotspot_r * sin_hotspot_theta * cos(Hotspot_phi);
    double y_center = Hotspot_r * sin_hotspot_theta * sin(Hotspot_phi);
    double z_center = Hotspot_r * cos(Hotspot_theta);

    double x_photon = photon_r * sin_photon_theta * cos(State_Vector[e_phi]);
    double y_photon = photon_r * sin_photon_theta * sin(State_Vector[e_phi]);
    double z_photon = photon_r * cos(State_Vector[e_theta]);

    double squred_distance_to_hotspot_center = (x_center - x_photon) * (x_center - x_photon) 
                                             + (y_center - y_photon) * (y_center - y_photon) 
                                             + (z_center - z_photon) * (z_center - z_photon);

    double Spatial_profile{};

    switch (this->s_Hotspot_params.Temperature_profile_type) {

    case e_Gaussian_profile:

        Spatial_profile = exp(-squred_distance_to_hotspot_center / Spatial_spread / Spatial_spread / 2);
        break;

    default:

        if (squred_distance_to_hotspot_center < this->s_Hotspot_params.Radius * this->s_Hotspot_params.Radius) {

            Spatial_profile = 1;

        }
        break;
    }

    double Temporal_profile = 1.0, Temporal_argument{};

    double& t_ref = this->s_Hotspot_params.Coord_time_at_max;
    double& t_sigma = this->s_Hotspot_params.Temporal_spread;

    if (0 != this->s_Hotspot_params.Temporal_spread) {

        Temporal_profile = exp(-(-State_Vector[e_t] - t_ref) * (-State_Vector[e_t] - t_ref) / t_sigma / t_sigma / 2);

    }

    double Hotspot_Temperature = this->s_Hotspot_params.Electron_temperature_scale * Spatial_profile * Temporal_profile;
    
    if (isnan(Hotspot_Temperature) || isinf(Hotspot_Temperature) || Hotspot_Temperature < 0) {

        std::cout << "Invalid hotspot temperature profile: " << Hotspot_Temperature << "\n";

        exit(ERROR);

    }

    return Hotspot_Temperature;

}

/* ====================================================== Velocity Functions ====================================================== */

//! Computes the emission medium's plasma 4-velocity
/*! Computes the emission medium's plasma 4-velocity
 *
 *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
 *   \param [in] p_Sim_Context - Pointer to the Simulation Context struct - used to call the metric function
 *   \param [in] Velocity_profile - Enum for the type of velocity profile
 *   \return Pointer to the 4-velocity vector
 */
double* Generic_Optically_Thin_Model::get_plasma_velocity(const double* const State_Vector, const Simulation_Context_type* const p_Sim_Context, Velocity_enums const Velocity_profile) {

    /* === Initialize some variables === */
    double Omega{}, rho{}, ell{}, u_t{}, u_phi{}, Normalization{}, inv_metric[4][4]{};
    static double Plasma_velocity[4]{};

    const double& r_source     = State_Vector[e_r];
    const double& theta_source = State_Vector[e_theta];

    Metric_type s_Metric = p_Sim_Context->p_Spacetime->get_metric(State_Vector);

    switch (Velocity_profile) {

    case e_Keplarian:

        Omega = sqrt(1 / r_source / r_source / r_source);
        Normalization = 1 / (-s_Metric.Metric[e_t][e_t] - 2 * s_Metric.Metric[e_t][e_phi] * Omega - s_Metric.Metric[e_phi][e_phi] * Omega * Omega);

        if (Normalization < 0) {

            return NULL;

        }

        Plasma_velocity[e_t] = sqrt(Normalization);
        Plasma_velocity[e_r] = 0.0;
        Plasma_velocity[e_theta] = 0.0;
        Plasma_velocity[e_phi] = Plasma_velocity[e_t] * Omega;

        break;

    default:

        rho = r_source * sin(theta_source);
        ell = sqrt(rho * rho * rho) / (1 + rho);

        if (Janis_Newman_Winicour == p_Sim_Context->p_Init_Conditions->Metric_params.e_Spacetime) {

            double& gamma = p_Sim_Context->p_Init_Conditions->Metric_params.JNW_Gamma_Parameter;
            double r_singularity = 2. / gamma;

            ell *= pow(1. - r_singularity / r_source, gamma);

        }
        else if (Wormhole == p_Sim_Context->p_Init_Conditions->Metric_params.e_Spacetime) {


            ell *= (1 - p_Sim_Context->p_Init_Conditions->Metric_params.R_throat / r_source);

        }

        invert_metric(inv_metric, s_Metric.Metric);

        u_t = -1.0 / sqrt(-(inv_metric[0][0] - 2 * inv_metric[0][3] * ell + inv_metric[3][3] * ell * ell));
        u_phi = -u_t * ell;

        /* Convert U_source to contravariant components */

        Plasma_velocity[e_t] = inv_metric[0][0] * u_t + inv_metric[0][3] * u_phi;
        Plasma_velocity[e_r] = 0.0;
        Plasma_velocity[e_theta] = 0.0;
        Plasma_velocity[e_phi] = inv_metric[3][3] * u_phi + inv_metric[3][0] * u_t;

        break;

    }

    if (isnan(Plasma_velocity[e_t]) ||
        isinf(Plasma_velocity[e_t]) ||
        isnan(Plasma_velocity[e_phi]) ||
        isinf(Plasma_velocity[e_phi])) {

        std::cout << "Invalid disk 4-velocity: "
            << "["
            << Plasma_velocity[e_t]
            << ", "
            << Plasma_velocity[e_r]
            << ", "
            << Plasma_velocity[e_theta]
            << ", "
            << Plasma_velocity[e_phi]
            << "]\n";

        exit(ERROR);

    }

    return Plasma_velocity;

}

/* ======================================================= Density Functions ====================================================== */

//! Computes the hotspot density
/*! Computes the hotspot density at the current photon position.
 *
 *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
 *   \return The density in [g/cm^3].
 */
double Generic_Optically_Thin_Model::get_hotspot_density(const double* const State_Vector) {

    double& Hotspot_r     = this->s_Hotspot_params.Position[e_r - 1];
    double& Hotspot_theta = this->s_Hotspot_params.Position[e_theta - 1];
    double& Hotspot_phi   = this->s_Hotspot_params.Position[e_phi - 1];

    double& Spatial_spread = this->s_Hotspot_params.Density_spread;

    const double& photon_r = State_Vector[e_r];

    // I need these more than once, so I precompute them here
    double sin_theta = sin(State_Vector[e_theta]);
    double sin_hotspot_theta = sin(Hotspot_theta);

    double x_center = Hotspot_r * sin_hotspot_theta * cos(Hotspot_phi);
    double y_center = Hotspot_r * sin_hotspot_theta * sin(Hotspot_phi);
    double z_center = Hotspot_r * cos(Hotspot_theta);

    double x_photon = photon_r * sin_theta * cos(State_Vector[e_phi]);
    double y_photon = photon_r * sin_theta * sin(State_Vector[e_phi]);
    double z_photon = photon_r * cos(State_Vector[e_theta]);

    double squred_distance_to_hotspot_center = (x_center - x_photon) * (x_center - x_photon)
                                             + (y_center - y_photon) * (y_center - y_photon)
                                             + (z_center - z_photon) * (z_center - z_photon);

    double Spatial_profile{};

    switch (this->s_Hotspot_params.Density_profile_type) {

    case e_Gaussian_profile:

        Spatial_profile = exp(-squred_distance_to_hotspot_center / Spatial_spread / Spatial_spread / 2);
        break;

    default:

        if (squred_distance_to_hotspot_center < this->s_Hotspot_params.Radius * this->s_Hotspot_params.Radius) {

            Spatial_profile = 1;

        }
        break;
    }

    double Temporal_profile = 1.0, Temporal_argument{};

    double& t_ref   = this->s_Hotspot_params.Coord_time_at_max;
    double& t_sigma = this->s_Hotspot_params.Temporal_spread;

    if (0 != this->s_Hotspot_params.Temporal_spread) {

        Temporal_profile = exp(-(-State_Vector[e_t] - t_ref) * (-State_Vector[e_t] - t_ref) / t_sigma / t_sigma / 2);

    }

    double Hotspot_Density = this->s_Hotspot_params.Electron_density_scale * Spatial_profile * Temporal_profile;

    if (isnan(Hotspot_Density) || isinf(Hotspot_Density) || Hotspot_Density < 0) {

        std::cout << "Invalid hotspot density profile: " << Hotspot_Density << "\n";

        exit(ERROR);

    }

    return Hotspot_Density;

}

//! Computes the background accretion disk density
/*! Computes the background accretion disk density at the current photon position.
 * 
 *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
 *   \return The density in [g/cm^3].
 */
double Generic_Optically_Thin_Model::get_disk_density(const double* const State_Vector) {

    const double& r  = State_Vector[e_r];
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

    double Disk_Density = this->s_Disk_params.Electron_density_scale * Disk_density_profile;

    if (isnan(Disk_Density) || isinf(Disk_Density) || Disk_Density < 0) {

        std::cout << "Invalid disk density profile: " << Disk_Density << "\n";

        exit(ERROR);

    }

    return Disk_Density;

}

/* =================================================== Disk Magnetic Field Functions =================================================== */

//! Computes the magnetic field 4-vector in the coordinate and plasma frames.
/*! Computes the magnetic field 4-vector, measured by a comoving obverver (with 4-velocity Plasma_velocity) in the following frames:
 *      1) That of a static observer (with 4-velocity n_mu = {1, 0, 0, 0} ) - a.e. the coordinate frame.
 *      2) The plasma rest frame.
 *
 *    NOTE: The magnitude of the magnetic field in these frames is different, because its not concerved under Lorentz boosts.
 *          In the plasma frame I set the geometry of the field, then scale it by B_Plasma_norm_CGS.
 *
 *    NOTE: The magnitudes of the magnetic fields in these frames are given in Gauss.
 * 
 *   \param [out] Magnetic_fields - Struct that holds the magnetic field 4-vector in the two frames.
 *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
 *   \param [in] p_Sim_Context - Pointer to the Simulation Context struct - used to call the metric function for dot products.
 *   \param [in] Density - The current emission medium density - used to compute the field magnitude in the plasma frame.
 *   \param [in] Magnetization - The current emission medium magnetization - used to compute the field magnitude in the plasma frame.
 *   \return Nothing.
 */
void Generic_Optically_Thin_Model::get_magnetic_field(Magnetic_fields_type* const Magnetic_fields,
                                                      const double* const State_Vector,
                                                      const Simulation_Context_type* const p_Sim_Context,
                                                      const double* const Plasma_Velocity,
                                                      const double Density,
                                                      const double Magnetization)  {

    Magnetic_fields->B_field_plasma_frame_norm = sqrt(Magnetization * C_LIGHT_CGS * C_LIGHT_CGS * Density * M_PROTON_CGS * 4 * M_PI);

    double Disk_B_plasma_frame[4] = { 0.0,
                                     Magnetic_fields->B_field_plasma_frame_norm * this->s_Disk_params.Mag_field_geometry[0],
                                     Magnetic_fields->B_field_plasma_frame_norm * this->s_Disk_params.Mag_field_geometry[1],
                                     Magnetic_fields->B_field_plasma_frame_norm * this->s_Disk_params.Mag_field_geometry[2] };

    /* 

    The only reference I could find for this is https://iopscience.iop.org/article/10.3847/1538-4357/ab718e/pdf 8e) and 8f).
    One arrives at the expression by projecting F^mu^nu (written in terms of B_coord_frame and Plasma_velocity),
    onto the 4-velocity of the static observer - this gives B_plasma_frame, then inverting the expression 
    to obtain B_coord_frame = f(B_plasma_frame, Plasma_Velocity)
    
    */

    Metric_type s_Metric = p_Sim_Context->p_Spacetime->get_metric(State_Vector);

    for (int left_idx = 0; left_idx <= 3; left_idx++) {

        for (int right_idx = 0; right_idx <= 3; right_idx++) {

            Magnetic_fields->B_field_coord_frame[e_t] += s_Metric.Metric[left_idx][right_idx] * Plasma_Velocity[left_idx] * Disk_B_plasma_frame[right_idx];
        }

    }

    for (int index = 1; index <= 3; index++) {

        Magnetic_fields->B_field_coord_frame[index] = (Disk_B_plasma_frame[index] + Magnetic_fields->B_field_coord_frame[e_t] * Plasma_Velocity[index]) / Plasma_Velocity[e_t];
       
    }

    //double Omega = Plasma_velocity[e_phi] / Plasma_velocity[e_t];
    //double A = 1. / sqrt(-(s_Metric.Metric[e_t][e_t] + Omega * Omega * s_Metric.Metric[e_phi][e_phi]));

    //B_coord_frame[e_t_coord] = A * sqrt(-s_Metric.Metric[e_phi_coord][e_phi_coord] / s_Metric.Metric[e_t_coord][e_t_coord]) * Omega;
    //B_coord_frame[e_r_coord] = 0;
    //B_coord_frame[e_theta_coord] = 0;
    //B_coord_frame[e_phi_coord] = A / sqrt(-s_Metric.Metric[e_phi_coord][e_phi_coord] / s_Metric.Metric[e_t_coord][e_t_coord]);

}

//! Computes the angle between the magnetic field and photon momentum 3-vectors in the plasma frame.
/*! Computes the angle between the magnetic field and photon momentum 3-vectors in the plasma frame. There is a neat invariant way 
 *   to do this by just operating on coordinate basis 4-vector using the projection tensor for an observer with 4-velocity = Plasma_velocity  .
 *
 *   \param [in] B_field_coord_frame - The magnetic field in the coordinate frame.
 *   \param [in] Plasma_velocity - The plasma velocity 4-vector.
 *   \param [in] State_Vector - Current photon state vector - used to get the photon momentum 4-vector.
 *   \param [in] p_Sim_Context - Pointer to the Simulation Context struct - used to call the metric function for dot products.
 *   \return Cosine of the angle between the magnetic field and photon momentum 3-vectors in the plasma frame.
 */
double Generic_Optically_Thin_Model::get_electron_pitch_angle(const double* const B_field_coord_frame, 
                                                              const double* const Plasma_velocity,
                                                              const double* const State_Vector, 
                                                              const Simulation_Context_type* const p_Sim_Context) {

    double Wave_vec_dot_Plasma_vec = State_Vector[e_p_t]     * Plasma_velocity[e_t] +
                                     State_Vector[e_p_r]     * Plasma_velocity[e_r] +
                                     State_Vector[e_p_theta] * Plasma_velocity[e_theta] +
                                     State_Vector[e_p_phi]   * Plasma_velocity[e_phi];

    Metric_type s_Metric  = p_Sim_Context->p_Spacetime->get_metric(State_Vector);
    double B_field_norm_squared{};
    double B_field_dot_Plasma_vel{};

    /*
    
    TODO: Maybe make functions that do this, or functions that raise and lower indicies
    
    */

    double test{};

    for (int left_idx = 0; left_idx <= 3; left_idx++) {

        for (int right_idx = 0; right_idx <= 3; right_idx++) {

            B_field_norm_squared   += s_Metric.Metric[left_idx][right_idx] * B_field_coord_frame[left_idx] * B_field_coord_frame[right_idx];
            B_field_dot_Plasma_vel += s_Metric.Metric[left_idx][right_idx] * B_field_coord_frame[left_idx] * Plasma_velocity[right_idx];

        }

    }

    double Wave_vec_dot_B_field = State_Vector[e_p_t]     * B_field_coord_frame[e_t] +
                                  State_Vector[e_p_r]     * B_field_coord_frame[e_r] +
                                  State_Vector[e_p_theta] * B_field_coord_frame[e_theta] +
                                  State_Vector[e_p_phi]   * B_field_coord_frame[e_phi];

    double cos_angle = 1.0; 

    if (!isinf(1.0 / Wave_vec_dot_Plasma_vec) && !isinf(1.0 / B_field_norm_squared)) {

        cos_angle = Wave_vec_dot_B_field / (fabs(Wave_vec_dot_Plasma_vec) * sqrt(B_field_norm_squared));

    }

    if (fabs(cos_angle) <= 1.0) {

        return acos(cos_angle);

    }
    else {

        return acos(cos_angle / fabs(cos_angle));

    }

}

/* =============================================== Thermal synchrotron Transfer Functions =============================================== */

//! Evaluates the thermal ensamble polarized synchrotron emission and Faradey functions.
/*! Evaluates the thermal ensamble polarized synchrotron emission and Faradey functions.
 *
 *   \param [in] Density - The current emission medium density in [g/cm^3].
 *   \param [in] T_electron_dim - The current emission medium dimentionless temperature.
 *   \param [in] f_cyclo - The current cyclotron frequency in [Hz].
 *   \param [in] sin_pitch_angle - The sine of the angle between the magnetic field and the photon momentum 3-vector in the plasma frame.
 *   \param [in] cos_pitch_angle - The cosine of the angle between the magnetic field and the photon momentum 3-vector in the plasma frame.
 *   \param [out] Emission_functions - Vector to hold the emission functions.
 *   \param [out] Faradey_functions - Vector to hold the Faradey functions.
 *   \param [in] Emission_args - Sturct to hold the arguments for evaluating the emission fit functions @see get_thermal_synchrotron_fit_functions.
 *   \param [in] Faradey_args - Struct to hold the aruments for evaluating the Faradey fit function @see get_thermal_synchrotron_fit_functions.
 *   \return Nothing.
 */
void Generic_Optically_Thin_Model::evaluate_thermal_synchrotron_transfer_functions(double Density,
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

        this->get_thermal_synchrotron_fit_functions(Emission_functions, Faradey_functions, &Emission_args, &Faradey_args);

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


//! Computes the necessary variables for evaluating the thermal ensamble polarized synchrotron transfer functions.
/*! Computes the necessary variables (cyclotron frequency, emission angles and so on) for evaluating the thermal ensamble polarized
 *   synchrotron transfer functions, based on the current photon position.
 *
 *   \param [in] State_Vector - The current photon state vector.
 *   \param [in] Plasma_velocity - The current emission medium plasma velocity.
 *   \param [in] p_Sim_Context - Pointer to the Simulation Context struct - used to access the initial conditions.
 *   \param [out] Emission_functions - Vector to hold the emission functions.
 *   \param [out] Faradey_functions - Vector to hold the Faradey functions.
 *   \param [out] Absorbtion_functions - Vector to hold the absorbtion functions.
 *   \param [in] Density - The current density of the emission medium.
 *   \param [in] Temperature - The current temperrature of the emission medium.
 *   \param [in] B_field_coord_frame - The magnetic field 4-vector in the coordinate frame.
 *   \param [in] B_field_plasma_frame_norm - The norm of the magnetic field in the plasma frame.
 *   \return Nothing.
 */
void Generic_Optically_Thin_Model::get_thermal_synchrotron_transfer_functions(const double* const State_Vector,
                                                                             const double* const Plasma_velocity,
                                                                             const Simulation_Context_type* const p_Sim_Context,
                                                                             double* const Emission_functions,
                                                                             double* const Faradey_functions,
                                                                             double* const Absorbtion_functions,
                                                                             double  const Density,
                                                                             double  const Temperature,
                                                                             double* const B_field_coord_frame,
                                                                             double  const B_field_plasma_frame_norm) {

    /* === Zero out the transfer functions just in case === */
    for (int stokes_index = 0; stokes_index <= STOKES_PARAM_NUM - 1; stokes_index++) {

        Emission_functions[stokes_index] = 0.0;
        Faradey_functions[stokes_index] = 0.0;
        Absorbtion_functions[stokes_index] = 0.0;

    }

    if (NULL == Plasma_velocity) {

        return;

    }

    double redshift = Redshift(State_Vector, Plasma_velocity, p_Sim_Context->p_Observer);

    if (isinf(redshift) || isnan(redshift) || isinf(1.0 / redshift)) {

        return;

    }

    /* Observation Frequency */
    double const obs_frequency = p_Sim_Context->p_Init_Conditions->Observer_params.obs_frequency;

    /* Dimentionless Electron Temperature */
    double const T_electron_dim = BOLTZMANN_CONST_CGS * Temperature / M_ELECTRON_CGS / C_LIGHT_CGS / C_LIGHT_CGS;

    /* Cyclotron Frequency */
    double const f_cyclo = Q_ELECTRON_CGS * B_field_plasma_frame_norm / (2 * M_PI * M_ELECTRON_CGS * C_LIGHT_CGS);

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

    int& Num_Samples_to_avg = p_Sim_Context->p_Init_Conditions->Emission_pitch_angle_samples_to_average;

    if (p_Sim_Context->p_Init_Conditions->Average_electron_pitch_angle) {

        /* ============ This loop averages over the emission pitch angle, which it gets from a pre-computed table ============ */

        for (int averaging_idx = 1; averaging_idx <= Num_Samples_to_avg - 1; averaging_idx++) {

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

            this->evaluate_thermal_synchrotron_transfer_functions(Density, T_electron_dim, f_cyclo, sin_pitch_angle, cos_pitch_angle, temp_emission_functions, temp_faradey_functions, Emission_args_ang_corrected, Faradey_args_ang_corrected);

            // The U component is 0 by definition
            Emission_functions[I] += temp_emission_functions[I] * sin_pitch_angle * M_PI / Num_Samples_to_avg / 2;
            Emission_functions[Q] += temp_emission_functions[Q] * sin_pitch_angle * M_PI / Num_Samples_to_avg / 2;
            Emission_functions[U] = 0.0;
            Emission_functions[V] += temp_emission_functions[V] * sin_pitch_angle * M_PI / Num_Samples_to_avg / 2;

            // The I and U components are 0 by definition
            Faradey_functions[I] = 0.0f;
            Faradey_functions[Q] += temp_faradey_functions[Q] * sin_pitch_angle * M_PI / Num_Samples_to_avg / 2;
            Faradey_functions[U] = 0.0f;
            Faradey_functions[V] += temp_faradey_functions[V] * sin_pitch_angle * M_PI / Num_Samples_to_avg / 2;

        }
    }
    else {

        /* The magnetic field is the one measured by a comoving with the plasma observer, but expressed in the cooridante frame */
        
        double pitch_angle = get_electron_pitch_angle(B_field_coord_frame, Plasma_velocity, State_Vector, p_Sim_Context);
        double sin_pitch_angle = sin(pitch_angle);

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

        this->evaluate_thermal_synchrotron_transfer_functions(Density, T_electron_dim, f_cyclo, sin_pitch_angle, cos(pitch_angle), Emission_functions, Faradey_functions, Emission_args_ang_corrected, Faradey_args_ang_corrected);

    }

    /* ================================================ The absorbtion functions ================================================ */

    double Planck_function_CGS = get_planck_function_CGS(obs_frequency / redshift, Temperature);

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

    /* Account for the relativistic doppler effet via the redshift */

    for (int stokes_idx = 0; stokes_idx <= STOKES_PARAM_NUM - 1; stokes_idx++) {

        Emission_functions[stokes_idx] *= redshift * redshift;
        Faradey_functions[stokes_idx] /= redshift;
        Absorbtion_functions[stokes_idx] /= redshift;
    }

}

/* ========================================== Kappa synchrotron Transfer Functions ========================================== */

//! Evaluates the kappa ensamble polarized synchrotron emission and Faradey functions.
/*! Evaluates the kappa ensamble polarized synchrotron emission and Faradey functions.
 *
 *   \param [in] Density - The current emission medium density in [g/cm^3].
 *   \param [in] f_cyclo - The current cyclotron frequency in [Hz].
 *   \param [out] Emission_functions - Vector to hold the emission functions.
 *   \param [out] Faradey_functions - Vector to hold the Faradey functions.
 *   \param [out] Absorbtion_functions - Vector to hold the absorbtion functions.
 *   \param [in] Transfer_args - Sturct to hold the arguments for evaluating the emission fit functions @see get_kappa_synchrotron_fit_functions.
 *   \return Nothing.
 */
void Generic_Optically_Thin_Model::evaluate_kappa_synchrotron_transfer_functions(double Density,
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

    this->get_kappa_synchrotron_fit_functions(Emission_functions, Faradey_functions, Absorbtion_functions, &Transfer_args);

    Emission_functions[I] *= Density * Q_ELECTRON_CGS * Q_ELECTRON_CGS / C_LIGHT_CGS * f_cyclo; 
    Emission_functions[Q] *= Density * Q_ELECTRON_CGS * Q_ELECTRON_CGS / C_LIGHT_CGS * f_cyclo; 
    Emission_functions[V] *= Density * Q_ELECTRON_CGS * Q_ELECTRON_CGS / C_LIGHT_CGS * f_cyclo; 

    if (!isnan(frequency)) {

        Absorbtion_functions[I] *= Density * Q_ELECTRON_CGS * Q_ELECTRON_CGS / M_ELECTRON_CGS / C_LIGHT_CGS / frequency;
        Absorbtion_functions[Q] *= Density * Q_ELECTRON_CGS * Q_ELECTRON_CGS / M_ELECTRON_CGS / C_LIGHT_CGS / frequency;
        Absorbtion_functions[V] *= Density * Q_ELECTRON_CGS * Q_ELECTRON_CGS / M_ELECTRON_CGS / C_LIGHT_CGS / frequency;

    }

}

//! Computes the necessary variables for evaluating the kappa ensamble polarized synchrotron transfer functions.
/*! Computes the necessary variables (cyclotron frequency, emission angles and so on) for evaluating the kappa ensamble polarized
 *   synchrotron transfer functions, based on the current photon position.
 *
 *   \param [in] State_Vector - The current photon state vector.
 *   \param [in] Plasma_velocity - The current emission medium plasma velocity.
 *   \param [in] p_Sim_Context - Pointer to the Simulation Context struct - used to access the initial conditions.
 *   \param [out] Emission_functions - Vector to hold the emission functions.
 *   \param [out] Faradey_functions - Vector to hold the Faradey functions.
 *   \param [out] Absorbtion_functions - Vector to hold the absorbtion functions.
 *   \param [in] Density - The current density of the emission medium.
 *   \param [in] Temperature - The current temperrature of the emission medium.
 *   \param [in] B_field_coord_frame - The magnetic field 4-vector in the coordinate frame.
 *   \param [in] B_field_plasma_frame_norm - The norm of the magnetic field in the plasma frame.
 *   \return Nothing.
 */
void Generic_Optically_Thin_Model::get_kappa_synchrotron_transfer_functions(const double* const State_Vector,
                                                                           const double* const Plasma_velocity,
                                                                           const Simulation_Context_type* const p_Sim_Context,
                                                                           double* const Emission_functions,
                                                                           double* const Faradey_functions,
                                                                           double* const Absorbtion_functions,
                                                                           double  const Density,
                                                                           double  const Temperature,
                                                                           double* const B_field_coord_frame,
                                                                           double  const B_field_plasma_frame_norm){

    /* === Zero out the transfer functions just in case === */
    for (int stokes_index = 0; stokes_index <= STOKES_PARAM_NUM - 1; stokes_index++) { 

        Emission_functions[stokes_index] = 0.0; 
        Faradey_functions[stokes_index] = 0.0; 
        Absorbtion_functions[stokes_index] = 0.0;  

    }

    if (NULL == Plasma_velocity) {

        return;

    }

    const double redshift = Redshift(State_Vector, Plasma_velocity, p_Sim_Context->p_Observer);

    if (isinf(redshift) || isnan(redshift) || isinf(1.0 / redshift)) {

        return;

    }

    /* Observation frequency */
    double& obs_frequency = p_Sim_Context->p_Init_Conditions->Observer_params.obs_frequency;

    /* Dimentionless Electron Temperature */
    double T_electron_dim = BOLTZMANN_CONST_CGS * Temperature / M_ELECTRON_CGS / C_LIGHT_CGS / C_LIGHT_CGS;

    /* Cyclotron Frequency */
    double f_cyclo = Q_ELECTRON_CGS * B_field_plasma_frame_norm / (2 * M_PI * M_ELECTRON_CGS * C_LIGHT_CGS);

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

    int& Num_Samples_to_avg = p_Sim_Context->p_Init_Conditions->Emission_pitch_angle_samples_to_average;

    if (p_Sim_Context->p_Init_Conditions->Average_electron_pitch_angle) {

        /* ============ This loop averages over the emission pitch angle, which it gets from a pre-computed table ============ */

        for (int averaging_idx = 1; averaging_idx <= Num_Samples_to_avg - 1; averaging_idx++) {

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

            this->evaluate_kappa_synchrotron_transfer_functions(Density, f_cyclo, temp_emission_functions, temp_faradey_functions, temp_absorbtion_functions, Transfer_args_ang_corrected);

            // The U component is 0 by definition
            Emission_functions[I] += temp_emission_functions[I] * sin_pitch_angle * M_PI / Num_Samples_to_avg / 2;
            Emission_functions[Q] += temp_emission_functions[Q] * sin_pitch_angle * M_PI / Num_Samples_to_avg / 2;
            Emission_functions[U] = 0.0;
            Emission_functions[V] += temp_emission_functions[V] * sin_pitch_angle * M_PI / Num_Samples_to_avg / 2;

            // The I and U components are 0 by definition
            Faradey_functions[I] = 0.0f;
            Faradey_functions[Q] += temp_faradey_functions[Q] * sin_pitch_angle * M_PI / Num_Samples_to_avg / 2;
            Faradey_functions[U] = 0.0f;
            Faradey_functions[V] += temp_faradey_functions[V] * sin_pitch_angle * M_PI / Num_Samples_to_avg / 2;

        }
    }
    else {

        /* The magnetic field is the one measured by a comoving with the plasma observer, but expressed in the cooridante frame */

        double pitch_angle = get_electron_pitch_angle(B_field_coord_frame, Plasma_velocity, State_Vector, p_Sim_Context);
        double sin_pitch_angle = sin(pitch_angle);

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
        this->evaluate_kappa_synchrotron_transfer_functions(Density, f_cyclo, Emission_functions, Faradey_functions, Absorbtion_functions, Transfer_args_ang_corrected);

    }

    /* Account for the relativistic doppler effet via the redshift */

    for (int stokes_idx = 0; stokes_idx <= STOKES_PARAM_NUM - 1; stokes_idx++) {

        Emission_functions[stokes_idx] *= redshift * redshift;
        Faradey_functions[stokes_idx] /= redshift;
        Absorbtion_functions[stokes_idx] /= redshift;
    }

}

/* ========================================== Phenomenological synchrotron Transfer Functions ========================================== */

//! Evaluates the phonomenological synchrotron transfer functions.
/*! Evaluates the phonomenological synchrotron transfer functions.
 *
 *   \param [in] State_Vector - The current photon state vector
 *   \param [in] Plasma_velocity - The current emission medium plasma velocity.
 *   \param [in] p_Sim_Context - Pointer to the Simulation Context struct - used to access the initial conditions.
 *   \param [out] Emission_functions - Vector to hold the emission functions.
 *   \param [out] Faradey_functions - Vector to hold the Faradey functions.
 *   \param [out] Absorbtion_functions - Vector to hold the absorbtion functions.
 *   \param [in] Density - The current density of the emission medium.
 *   \return Nothing.
 */
void Generic_Optically_Thin_Model::get_phenomenological_synchrotron_functions(const double* const State_Vector,
                                                                             const double* const Plasma_Veclocity,
                                                                             const Simulation_Context_type* const p_Sim_Context, 
                                                                             double* const Emission_functions,
                                                                             double* const Faradey_functions,
                                                                             double* const Absorbtion_functions,
                                                                             const double Density) {

    if (NULL == Plasma_Veclocity) {

        return;

    }

    double redshift = Redshift(State_Vector, Plasma_Veclocity, p_Sim_Context->p_Observer);

    if (isinf(redshift) || isnan(redshift) || isinf(1.0 / redshift)) {

        return;

    }

    double& emission_power_law = this->s_Emission_params.Phenomenological_emission_power_law;
    double& source_f_power_law = this->s_Emission_params.Phenomenological_source_f_power_law;
    double& emission_coeff     = this->s_Emission_params.Phenomenological_emission_coeff;
    double& abs_coeff          = this->s_Emission_params.Phenomenological_absorbtion_coeff;

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

    /* Account for the relativistic doppler effet via the redshift */

    for (int stokes_idx = 0; stokes_idx <= STOKES_PARAM_NUM - 1; stokes_idx++) {

        Emission_functions[stokes_idx] *= redshift * redshift;
        Faradey_functions[stokes_idx] /= redshift;
        Absorbtion_functions[stokes_idx] /= redshift;
    }

}

/* ============================================ Main "Selector" For The Transfer Functions ============================================ */

//! Main "Selector" For The Transfer Functions.
/*! Calculates the density, temperature, magnetic field and 4-velocity of the chosen emission medium and calls the respective transfer functions evaluation.
 *
 *   \param [in] State_Vector - The current photon state vector.
 *   \param [in] p_Sim_Context - Pointer to the Simulation Context struct - used to access the initial conditions.
 *   \param [out] Emission_functions - Vector to hold the emission functions.
 *   \param [out] Faradey_functions - Vector to hold the Faradey functions.
 *   \param [out] Absorbtion_functions - Vector to hold the absorbtion functions.
 *   \param [in] Emission_medium - Enum that specifies which emission medium to evaluate.
 *   \return Nothing.
 */
void Generic_Optically_Thin_Model::get_radiative_transfer_functions(const double* const State_Vector,
                                                                    const Simulation_Context_type* const p_Sim_Context, 
                                                                    double* const Emission_functions,
                                                                    double* const Faradey_functions,
                                                                    double* const Absorbtion_functions,
                                                                    Emission_medium_enums Emission_medium) {

    for (int index = I; index <= STOKES_PARAM_NUM - 1; index++) {

        Emission_functions[index] = 0.0;
        Faradey_functions[index] = 0.0;
        Absorbtion_functions[index] = 0.0;

    }

    Magnetic_fields_type Magnetic_fields{};
    Ensamble_enums Ensamble_type{};

    double Density{}, Temperature{}, B_field_norm_plasma_frame{}, Magnetization{};
    double* Plasma_Velocity{};
    double* B_field_coord_frame{};

    switch (Emission_medium) {

    case Disk:

        Density         = this->get_disk_density(State_Vector);
        Temperature     = this->get_disk_temperature(State_Vector);
        Plasma_Velocity = this->get_plasma_velocity(State_Vector, p_Sim_Context, this->s_Disk_params.Velocity_profile_type);
        Ensamble_type   = p_Sim_Context->p_Init_Conditions->Disk_params.Ensamble_type;
        Magnetization   = p_Sim_Context->p_Init_Conditions->Disk_params.Magnetization;

        this->get_magnetic_field(&Magnetic_fields, State_Vector, p_Sim_Context, Plasma_Velocity, Density, Magnetization);

        break;

    case Hotspot:

        Density         = this->get_hotspot_density(State_Vector);
        Temperature     = this->get_hotspot_temperature(State_Vector);
        Plasma_Velocity = this->get_plasma_velocity(State_Vector, p_Sim_Context, this->s_Hotspot_params.Velocity_profile_type);
        Ensamble_type   = p_Sim_Context->p_Init_Conditions->Hotspot_params.Ensamble_type;
        Magnetization   = p_Sim_Context->p_Init_Conditions->Hotspot_params.Magnetization;

        this->get_magnetic_field(&Magnetic_fields, State_Vector, p_Sim_Context, Plasma_Velocity, Density, Magnetization);

        break;

    default:

        std::cout << "Unsupported emissison medium - something broke in the get_radiative_transfer_functions function!" << "\n";

        exit(ERROR);

        break;

    }

    B_field_coord_frame       = Magnetic_fields.B_field_coord_frame;
    B_field_norm_plasma_frame = Magnetic_fields.B_field_plasma_frame_norm;

    switch (Ensamble_type) {

    case(e_Phenomenological_ensamble):

        this->get_phenomenological_synchrotron_functions(State_Vector, Plasma_Velocity, p_Sim_Context, Emission_functions, Faradey_functions, Absorbtion_functions, Density);
        break;

    case(e_Kappa_ensamble):

        this->get_kappa_synchrotron_transfer_functions(State_Vector, Plasma_Velocity, p_Sim_Context, Emission_functions, Faradey_functions, Absorbtion_functions,
                                                      Density, Temperature, B_field_coord_frame, B_field_norm_plasma_frame);
        break;

    default:

        this->get_thermal_synchrotron_transfer_functions(State_Vector, Plasma_Velocity, p_Sim_Context, Emission_functions, Faradey_functions, Absorbtion_functions,
                                                       Density, Temperature, B_field_coord_frame, B_field_norm_plasma_frame);
        break;
    }
  
}

/* ========================================================== Misc Functions ========================================================== */

//! Precomputes the electron pitch angles and their weird powers to use in averaging.
/*! Precomputes the electron pitch angles and their weird powers to use in averaging.
 *
 *   \param [in] p_Init_Conditions - Pointer to the struct that holds the initial conditions - used to determine how much memory to allocate.
 *   \return Nothing.
 */
void Generic_Optically_Thin_Model::precompute_electron_pitch_angles(Initial_conditions_type* p_Init_Conditions) {

    // ====================================================== Allocate memory for the arrays ====================================================== //

    this->s_Precomputed_e_pitch_angles.cos_electron_pitch_angles = new double[p_Init_Conditions->Emission_pitch_angle_samples_to_average];
    this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles = new double[p_Init_Conditions->Emission_pitch_angle_samples_to_average];

    this->s_Precomputed_e_pitch_angles.one_over_sqrt_sin = new double[p_Init_Conditions->Emission_pitch_angle_samples_to_average];
    this->s_Precomputed_e_pitch_angles.one_over_cbrt_sin = new double[p_Init_Conditions->Emission_pitch_angle_samples_to_average];

    this->s_Precomputed_e_pitch_angles.one_over_sin_to_1_point_035      = new double[p_Init_Conditions->Emission_pitch_angle_samples_to_average];
    this->s_Precomputed_e_pitch_angles.one_over_sin_to_1_point_2_over_2 = new double[p_Init_Conditions->Emission_pitch_angle_samples_to_average];

    this->s_Precomputed_e_pitch_angles.one_over_sin_to_7_over_20 = new double[p_Init_Conditions->Emission_pitch_angle_samples_to_average];

    // =========================================================================================================================================== //

    for (int index = 0; index <= p_Init_Conditions->Emission_pitch_angle_samples_to_average - 1; index++) {

        double pitch_angle = double(index) / p_Init_Conditions->Emission_pitch_angle_samples_to_average * M_PI;
        this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index] = sin(pitch_angle);
        this->s_Precomputed_e_pitch_angles.cos_electron_pitch_angles[index] = cos(pitch_angle);

        if (this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index] != 0) {

            // Used in the thermal synchrotron emission functions

            this->s_Precomputed_e_pitch_angles.one_over_sqrt_sin[index] = 1. / sqrt(this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index]);
            this->s_Precomputed_e_pitch_angles.one_over_cbrt_sin[index] = 1. / cbrt(this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index]);

            // Used in the thermal synchrotron Faradey functions

            this->s_Precomputed_e_pitch_angles.one_over_sin_to_1_point_035[index] = 1. / pow(this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index], 1.035);
            this->s_Precomputed_e_pitch_angles.one_over_sin_to_1_point_2_over_2[index] = 1. / pow(this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index], 1.2 / 2);

            // Used in the kappa synchrotron emission functions

            this->s_Precomputed_e_pitch_angles.one_over_sin_to_7_over_20[index] = 1. / pow(this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index], 7. / 20);

        }
    }
}


//!  Copies over the initial data from the Simulation Context struct to internal class variables for the sake of convenience
/*!< Copies over the initial data from the Simulation Context struct to internal class variables for the sake of convenience
 *
 *   \param [in] p_Sim_Context - Pointer to the Simulation Context struct.
 *   \return Nothing.
 */
int Generic_Optically_Thin_Model::load_parameters(Simulation_Context_type* p_Sim_Context) {

    this->Num_samples_to_avg = p_Sim_Context->p_Init_Conditions->Emission_pitch_angle_samples_to_average;
    this->Include_polarization = p_Sim_Context->p_Init_Conditions->Observer_params.include_polarization;

    if (NULL != &p_Sim_Context->p_Init_Conditions->Disk_params) {

        this->s_Disk_params = p_Sim_Context->p_Init_Conditions->Disk_params;

    }
    else {

        std::cout << "p_Disk_params is a NULL pointer! \n";

        exit(ERROR);
    }

    if (NULL != &p_Sim_Context->p_Init_Conditions->Hotspot_params) {

        this->s_Hotspot_params = p_Sim_Context->p_Init_Conditions->Hotspot_params;

    }
    else {

        std::cout << "p_Hotspot_params is a NULL pointer! \n";

        exit(ERROR);
    }

    if (NULL != &p_Sim_Context->p_Init_Conditions->Emission_params) {

        this->s_Emission_params = p_Sim_Context->p_Init_Conditions->Emission_params;

    }
    else {

        std::cout << "p_Emission_params is a NULL pointer! \n";

        exit(ERROR);
    }

    return OK;

}


//! Evaluates the Planck function in the frequency domain in CGS units
/*! Evaluates the Planck function in the frequency domain in CGS units
 *
 *   \param [in] Frequency - Emission frequency in [Hz].
 *   \param [in] Temperature - Emission medium temperature in [K].
 *   \return The value of the Planck function in CGS.
 */
double get_planck_function_CGS(double Frequency, double Temperature) {

    return 2 * PLANCK_CONSTANT_CGS * Frequency * Frequency * Frequency / C_LIGHT_CGS / C_LIGHT_CGS / (exp(PLANCK_CONSTANT_CGS * Frequency / BOLTZMANN_CONST_CGS / Temperature) - 1.);

}
