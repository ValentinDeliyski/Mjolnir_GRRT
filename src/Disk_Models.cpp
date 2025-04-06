#include "Disk_Models.h"

/***************************************************
|                                                  |
| Novikov-Thorne Model Class Functions Definitions |
|                                                  | 
***************************************************/

Novikov_Thorne_Model::Novikov_Thorne_Model(Simulation_Context_type* p_Sim_Context) {

    this->r_in  = p_Sim_Context->p_Init_Conditions->NT_params.r_in;
    this->r_out = p_Sim_Context->p_Init_Conditions->NT_params.r_out;
    this->flux_integral_accuracy = p_Sim_Context->p_Init_Conditions->Integrator_params.Simpson_accuracy;
    this->p_Spacetime = p_Sim_Context->p_Spacetime;
    this->e_Spacetime = p_Sim_Context->p_Init_Conditions->Metric_parameters.e_Spacetime;

};

double Novikov_Thorne_Model::Keplerian_angular_velocity(const double* const State_Vector) {

    Metric_type s_dr_Metric = this->p_Spacetime->get_dr_metric(State_Vector);

    return (-s_dr_Metric.Metric[0][3] + sqrt(s_dr_Metric.Metric[0][3] * s_dr_Metric.Metric[0][3] - s_dr_Metric.Metric[0][0] * s_dr_Metric.Metric[3][3])) / s_dr_Metric.Metric[3][3];

}

double Novikov_Thorne_Model::dr_Keplerian_angular_velocity(const double* const State_Vector) {

    Metric_type s_dr_Metric = this->p_Spacetime->get_dr_metric(State_Vector);
    Metric_type s_d2r_Metric = this->p_Spacetime->get_d2r_metric(State_Vector);

    double root = sqrt(s_dr_Metric.Metric[0][3] * s_dr_Metric.Metric[0][3] - s_dr_Metric.Metric[0][0] * s_dr_Metric.Metric[3][3]);

    double Kepler = this->Keplerian_angular_velocity(State_Vector);

    return  - Kepler / s_dr_Metric.Metric[3][3] * s_d2r_Metric.Metric[3][3] + (-s_d2r_Metric.Metric[0][3]
            + 1.0 / root / 2 * (2 * s_dr_Metric.Metric[0][3] * s_d2r_Metric.Metric[0][3] - s_dr_Metric.Metric[0][0] * s_d2r_Metric.Metric[3][3]
            - s_d2r_Metric.Metric[0][0] * s_dr_Metric.Metric[3][3])) / s_dr_Metric.Metric[3][3];

}

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
                  << "Omega = "
                  << "\n"
                  << Kepler
                  << "\n";

        exit(ERROR);

    }

    double U_source[4] = { Gamma, 0, 0, Gamma * Kepler };

    return  (-U_obs[0] + U_obs[3] * State_Vector[e_p_phi]) / (-U_source[0] + U_source[3] * State_Vector[e_p_phi]);

}

double Novikov_Thorne_Model::disk_Energy(const double* const State_Vector) {

    Metric_type s_Metric_source = this->p_Spacetime->get_metric(State_Vector);

    double Kepler = this->Keplerian_angular_velocity(State_Vector);

    double root = sqrt(-s_Metric_source.Metric[0][0] - 2 * s_Metric_source.Metric[0][3] * Kepler - s_Metric_source.Metric[3][3] * Kepler * Kepler);

    return  -(s_Metric_source.Metric[0][0] + s_Metric_source.Metric[0][3] * Kepler) / root;

}

double Novikov_Thorne_Model::disk_Angular_Momentum(const double* const State_Vector) {

    Metric_type s_Metric_source = this->p_Spacetime->get_metric(State_Vector);

    double Kepler = this->Keplerian_angular_velocity(State_Vector);

    double root = sqrt(-s_Metric_source.Metric[0][0] - 2 * s_Metric_source.Metric[0][3] * Kepler - s_Metric_source.Metric[3][3] * Kepler * Kepler);

    return  (s_Metric_source.Metric[3][3] * Kepler + s_Metric_source.Metric[0][3]) / root;

}

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

Hotspot_position_type Generic_Optically_Thin_Model::get_hotspot_position(const double* const State_Vector,
                                                                         const Simulation_Context_type* const p_Sim_Context){

    // The velocity is needed to integrrate the position of the hospot. 
    double* Hotspot_velocity = this->get_plasma_velocity(this->s_Hotspot_params.Position - 1,
                                                         p_Sim_Context, 
                                                         p_Sim_Context->p_Init_Conditions->Hotspot_params.Velocity_profile_type, 
                                                         p_Sim_Context->p_Init_Conditions->Hotspot_params.Radial_velocity_fraction);
    Hotspot_position_type Hotspot_position{};
    double Hotspot_ang_velocity{};

    if (NULL != Hotspot_velocity) { Hotspot_ang_velocity = Hotspot_velocity[e_phi] / Hotspot_velocity[e_t]; }

    Hotspot_position.Distance    = this->s_Hotspot_params.Position[e_r - 1];
    Hotspot_position.Inclination = M_PI_2;
    Hotspot_position.Azimuth     = this->s_Hotspot_params.Position[e_phi - 1] + Hotspot_ang_velocity * (-State_Vector[e_t] - this->s_Hotspot_params.Coord_time_offset);

    double sin_hotspot_inclination = sin(Hotspot_position.Inclination);

    Hotspot_position.x = Hotspot_position.Distance * sin_hotspot_inclination * cos(Hotspot_position.Azimuth);
    Hotspot_position.y = Hotspot_position.Distance * sin_hotspot_inclination * sin(Hotspot_position.Azimuth);
    Hotspot_position.z = Hotspot_position.Distance * cos(Hotspot_position.Inclination);

    return Hotspot_position;

}

/* ==================================================== Temperature Functions ===================================================== */

double Generic_Optically_Thin_Model::get_disk_temperature(const double* const State_Vector) const {


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

    return T_scale * Disk_temperature_profile;

}

double Generic_Optically_Thin_Model::get_hotspot_temperature(const double* const State_Vector, const Simulation_Context_type* const p_Sim_Context) {

    Hotspot_position_type Hotspot_position = this->get_hotspot_position(State_Vector, p_Sim_Context);

    const double& Hotspot_r      = Hotspot_position.Distance;
    const double& Hotspot_theta  = Hotspot_position.Inclination;
    const double& Hotspot_phi    = Hotspot_position.Azimuth;
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

    double& t_ref = this->s_Hotspot_params.Coord_time_offset;
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

double* Generic_Optically_Thin_Model::get_plasma_velocity(const double* const State_Vector, 
                                                          const Simulation_Context_type* const p_Sim_Context, 
                                                          Velocity_enums const Velocity_profile,
                                                          double const Radial_velocity_fraction) {

    /* The reference for this implementation is https://arxiv.org/pdf/2206.12066. */

    /* === Initialize some variables === */
    double Omega{}, rho{}, ell{}, u_t{}, u_r{}, u_phi{}, Normalization{}, inv_metric[4][4]{};
    static double Plasma_velocity[4]{};

    const double& r_source     = State_Vector[e_r];
    const double& theta_source = State_Vector[e_theta];

    Metric_type s_Metric = p_Sim_Context->p_Spacetime->get_metric(State_Vector);
    invert_metric(inv_metric, s_Metric.Metric);

    Metric_type s_dr_Metric = p_Sim_Context->p_Spacetime->get_dr_metric(State_Vector);

    switch (Velocity_profile) {

    case e_Keplarian:

        /* This velocity profile is defined only for orbit radii > ISCO. When the ray passes below ISCO I return a NULL pointer, which tells the 
           rest of the code to ignore the emission from this region. */
        if (fabs(State_Vector[e_r]) < p_Sim_Context->p_Spacetime->get_ISCO()[Inner]){ return NULL; }

        /* Interpolated contravariant radial velocity component -> Corresponds to equation (10a) from the reference, but beta_r -> 1 - beta_r. */
        u_r = -Radial_velocity_fraction * sqrt((-1 - inv_metric[e_t][e_t]) * inv_metric[e_r][e_r]);

        /* Interpolated azimuthal angular velocity -> Corresponds to equation (10b) from the reference, but with beta_phi = 1 - beta_r. */
        Omega = -s_dr_Metric.Metric[e_t][e_phi] / s_dr_Metric.Metric[e_phi][e_phi];
        Omega += sqrt(s_dr_Metric.Metric[e_t][e_phi] * s_dr_Metric.Metric[e_t][e_phi] - s_dr_Metric.Metric[e_t][e_t] * s_dr_Metric.Metric[e_phi][e_phi]) / s_dr_Metric.Metric[e_phi][e_phi];
        Omega = Omega + Radial_velocity_fraction * (inv_metric[e_t][e_phi] / inv_metric[e_t][e_t] - Omega);

        break;

    default:

        rho = r_source * sin(theta_source);
        ell = sqrt(rho * rho * rho) / (1 + rho);

        /* I have noticed that this velocity profile becomes ill-defined in some places for the metric in the below "if" clause. 
           I correct this by modifying the angular momentum profile by something that seems reasonable. */
        if (Janis_Newman_Winicour == p_Sim_Context->p_Init_Conditions->Metric_parameters.e_Spacetime) {

            double& gamma = p_Sim_Context->p_Init_Conditions->Metric_parameters.JNW_Gamma_Parameter;
            ell *= pow(1. -  2. / r_source / gamma, gamma);

        }
        else if (Wormhole == p_Sim_Context->p_Init_Conditions->Metric_parameters.e_Spacetime) {

            ell *= (1. - p_Sim_Context->p_Init_Conditions->Metric_parameters.R_throat / r_source);

        }

        u_t   = -1.0 / sqrt(-(inv_metric[e_t][e_t] - 2 * inv_metric[e_t][e_phi] * ell + inv_metric[e_phi][e_phi] * ell * ell));
        u_phi = -u_t * ell;

        /* Convert U_source to contravariant components to compute the circular velocity profile*/
        Plasma_velocity[e_t] = inv_metric[e_t][e_t] * u_t + inv_metric[e_t][e_phi] * u_phi;
        Plasma_velocity[e_r] = 0.0;
        Plasma_velocity[e_theta] = 0.0;
        Plasma_velocity[e_phi] = inv_metric[e_phi][e_phi] * u_phi + inv_metric[e_phi][e_t] * u_t;

        /* Interpolated contravariant radial velocity component -> Corresponds to equation (10a) from the reference, but beta_r -> 1 - beta_r. */
        u_r = -Radial_velocity_fraction * sqrt((-1 - inv_metric[e_t][e_t]) * inv_metric[e_r][e_r]);

        /* Interpolated azimuthal angular velocity -> Corresponds to equation (10b) from the reference, but with beta_phi = 1 - beta_r. */
        Omega = Plasma_velocity[e_phi] / Plasma_velocity[e_t] + Radial_velocity_fraction * (inv_metric[e_t][e_phi] / inv_metric[e_t][e_t] - Plasma_velocity[e_phi] / Plasma_velocity[e_t]);

        break;

    }

    /* Interpolate between the circular and radial velocity profile */
    Normalization = -1 / (s_Metric.Metric[e_t][e_t] + 2 * s_Metric.Metric[e_t][e_phi] * Omega + s_Metric.Metric[e_phi][e_phi] * Omega * Omega);

    Plasma_velocity[e_t]     = sqrt((1 + s_Metric.Metric[e_r][e_r] * u_r * u_r) * Normalization);
    Plasma_velocity[e_r]     = u_r;
    Plasma_velocity[e_theta] = 0.0;
    Plasma_velocity[e_phi]   = Plasma_velocity[e_t] * Omega;

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

        return NULL;

    }

    return Plasma_velocity;

}

double Generic_Optically_Thin_Model::get_hotspot_density(const double* const State_Vector, const Simulation_Context_type* const p_Sim_Context) {

    Hotspot_position_type Hotspot_position = this->get_hotspot_position(State_Vector, p_Sim_Context);

    const double& Hotspot_r = Hotspot_position.Distance;
    const double& Hotspot_theta = Hotspot_position.Inclination;
    const double& Hotspot_phi = Hotspot_position.Azimuth;

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

            Spatial_profile = 1.0;

        }
        break;
    }

    double Temporal_profile = 1.0, Temporal_argument{};

    double& t_ref   = this->s_Hotspot_params.Coord_time_offset;
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

void Generic_Optically_Thin_Model::get_magnetic_field(const double* const State_Vector,
                                                      const Simulation_Context_type* const p_Sim_Context,
                                                      Emission_medium_state_type* const Emission_medium_state)  {

    Emission_medium_state->Magnetic_fields.B_field_plasma_frame_norm = sqrt(Emission_medium_state->Magnetization * C_LIGHT_CGS * C_LIGHT_CGS * Emission_medium_state->Density * M_PROTON_CGS * 4 * M_PI);

    Emission_medium_state->Magnetic_fields.B_field_plasma_frame[e_t]     = 0.0;
    Emission_medium_state->Magnetic_fields.B_field_plasma_frame[e_r]     = Emission_medium_state->Magnetic_fields.B_field_plasma_frame_norm * Emission_medium_state->Magnetic_fields.Magnetic_field_geometry[e_r - 1];
    Emission_medium_state->Magnetic_fields.B_field_plasma_frame[e_theta] = Emission_medium_state->Magnetic_fields.B_field_plasma_frame_norm * Emission_medium_state->Magnetic_fields.Magnetic_field_geometry[e_theta - 1];
    Emission_medium_state->Magnetic_fields.B_field_plasma_frame[e_phi]   = Emission_medium_state->Magnetic_fields.B_field_plasma_frame_norm * Emission_medium_state->Magnetic_fields.Magnetic_field_geometry[e_phi - 1];

    /* 

    The only reference I could find for this is https://iopscience.iop.org/article/10.3847/1538-4357/ab718e/pdf 8e) and 8f).
    One arrives at the expression by projecting F^mu^nu (written in terms of B_coord_frame and Plasma_velocity),
    onto the 4-velocity of the static observer - this gives B_plasma_frame, then inverting the expression 
    to obtain B_coord_frame = f(B_plasma_frame, Plasma_Velocity)
    
    */

    Metric_type s_Metric = p_Sim_Context->p_Spacetime->get_metric(State_Vector);

    for (int left_idx = 0; left_idx <= 3; left_idx++) {

        for (int right_idx = 0; right_idx <= 3; right_idx++) {

            Emission_medium_state->Magnetic_fields.B_field_coord_frame[e_t] += s_Metric.Metric[left_idx][right_idx] * Emission_medium_state->Plasma_Velocity[left_idx] * Emission_medium_state->Magnetic_fields.B_field_plasma_frame[right_idx];
        }

    }

    for (int index = 1; index <= 3; index++) {

        Emission_medium_state->Magnetic_fields.B_field_coord_frame[index] = (Emission_medium_state->Magnetic_fields.B_field_plasma_frame[index] + Emission_medium_state->Magnetic_fields.B_field_coord_frame[e_t] * Emission_medium_state->Plasma_Velocity[index]) / Emission_medium_state->Plasma_Velocity[e_t];
       
    }

    //double Omega = Plasma_velocity[e_phi] / Plasma_velocity[e_t];
    //double A = 1. / sqrt(-(s_Metric.Metric[e_t][e_t] + Omega * Omega * s_Metric.Metric[e_phi][e_phi]));

    //B_coord_frame[e_t_coord] = A * sqrt(-s_Metric.Metric[e_phi_coord][e_phi_coord] / s_Metric.Metric[e_t_coord][e_t_coord]) * Omega;
    //B_coord_frame[e_r_coord] = 0;
    //B_coord_frame[e_theta_coord] = 0;
    //B_coord_frame[e_phi_coord] = A / sqrt(-s_Metric.Metric[e_phi_coord][e_phi_coord] / s_Metric.Metric[e_t_coord][e_t_coord]);

}

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

void Generic_Optically_Thin_Model::get_thermal_synchrotron_transfer_functions(const double* const State_Vector,
                                                                              const Simulation_Context_type* const p_Sim_Context,
                                                                              const Emission_medium_state_type* const p_Emission_medium_state,
                                                                              Transfer_functions_type* const p_Transfer_functions) {

    /* === Zero out the transfer functions just in case === */
    memset(p_Transfer_functions, 0, sizeof(Transfer_functions_type));

    double redshift = get_redshift(State_Vector, p_Emission_medium_state->Plasma_Velocity, p_Sim_Context->p_Observer);

    /* Check weather redshift is numerically OK to use in the transfer functions. */
    if (isinf(redshift) || isnan(redshift) || isinf(1.0 / redshift)) { return; }

    /* The dimensionless electron temperature. */
    double const T_electron_dim = BOLTZMANN_CONST_CGS * p_Emission_medium_state->Temperature / M_ELECTRON_CGS / C_LIGHT_CGS / C_LIGHT_CGS;

    /* The cyclotron frequency. */
    double const f_cyclo = Q_ELECTRON_CGS * p_Emission_medium_state->Magnetic_fields.B_field_plasma_frame_norm / (2 * M_PI * M_ELECTRON_CGS * C_LIGHT_CGS);

    /* The "averaged" rescaled (by a factor of 27 / 4) critical frequency (without the sin(theta) term.
       That gets added on later from a pre-computed table). */
    double const f_s_no_sin = 2. / 9 * f_cyclo * T_electron_dim * T_electron_dim;

    /* Check weather the rescaled critical frequency f_s is numerically OK to use in the transfer functions. */
    if (isinf(f_s_no_sin) || isnan(f_s_no_sin) || isinf(1.0 / f_s_no_sin)) { return; }

    Thermal_transfer_f_arguments_type Transfer_args_uncorrected{};

    /* Observation Frequency */
    double const obs_frequency = p_Sim_Context->p_Init_Conditions->Observer_params.obs_frequency;

    /* Compute all the wierd powers of X outside the pitch angle averaging loop. */
    Transfer_args_uncorrected.X      = obs_frequency / f_s_no_sin / redshift;
    Transfer_args_uncorrected.sqrt_X = sqrt(Transfer_args_uncorrected.X);
    Transfer_args_uncorrected.cbrt_X = cbrt(Transfer_args_uncorrected.X);
    Transfer_args_uncorrected.X_to_0_p_5175 = pow(Transfer_args_uncorrected.X, 0.5175);
    Transfer_args_uncorrected.X_to_0_p_6    = pow(Transfer_args_uncorrected.X, 0.6);
    Transfer_args_uncorrected.X_to_0_p_7515 = pow(Transfer_args_uncorrected.X, 0.7515);

    /* Compute the werid power of the electron temperature outside the pitch angle averaging loop. */
    Transfer_args_uncorrected.T_electron_dim          = T_electron_dim;
    Transfer_args_uncorrected.T_electron_dim_to_24_25 = pow(T_electron_dim, 24. / 25);

    /* We also need the current photon frequency to evalauate the emission functions. */
    Transfer_args_uncorrected.frequency = obs_frequency / redshift;

    /* Reference for the sake of readabiity. */
    int& Num_Samples_to_avg = p_Sim_Context->p_Init_Conditions->Emission_pitch_angle_samples_to_average;

    /* This structcs holds the transfer function args, corrected for the electron pitch angle. I am setting it equal to the uncorrected one
        so they copy over the variables that don't depend on the pitch angle. The rest get corrected inside the averaging loop. */
    Thermal_transfer_f_arguments_type Transfer_args_corrected = Transfer_args_uncorrected;

    if (p_Sim_Context->p_Init_Conditions->Average_electron_pitch_angle) {

        /* ============ This loop averages over the emission pitch angle, which it gets from a pre-computed table ============ */

        for (int averaging_idx = 1; averaging_idx <= Num_Samples_to_avg - 1; averaging_idx++) {

            /* References for the sake of readabiity. */
            double& sin_pitch_angle = this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[averaging_idx];
            double& cos_pitch_angle = this->s_Precomputed_e_pitch_angles.cos_electron_pitch_angles[averaging_idx];

            Transfer_args_corrected.X      = Transfer_args_uncorrected.X / sin_pitch_angle;
            Transfer_args_corrected.sqrt_X = Transfer_args_uncorrected.sqrt_X * this->s_Precomputed_e_pitch_angles.one_over_sqrt_sin[averaging_idx];
            Transfer_args_corrected.cbrt_X = Transfer_args_uncorrected.cbrt_X * this->s_Precomputed_e_pitch_angles.one_over_cbrt_sin[averaging_idx];
            Transfer_args_corrected.X_to_0_p_5175 = Transfer_args_uncorrected.X_to_0_p_5175 * this->s_Precomputed_e_pitch_angles.one_over_sin_to_0_p_5175[averaging_idx];
            Transfer_args_corrected.X_to_0_p_6    = Transfer_args_uncorrected.X_to_0_p_6 * this->s_Precomputed_e_pitch_angles.one_over_sin_to_0_p_6[averaging_idx];
            Transfer_args_corrected.X_to_0_p_7515 = Transfer_args_uncorrected.X_to_0_p_7515 * this->s_Precomputed_e_pitch_angles.one_over_sin_to_0_p_7515[averaging_idx];

            Transfer_args_corrected.sin_pitch_angle = sin_pitch_angle;
            Transfer_args_corrected.cos_pitch_angle = cos_pitch_angle;

            Transfer_functions_type temp_Transfer_functions{};

            this->evaluate_synchrotron_transfer_functions(e_Thermal_ensamble, p_Emission_medium_state, &Transfer_args_corrected, p_Sim_Context, &temp_Transfer_functions);

            // The U component is 0 by definition
            p_Transfer_functions->Emission_functions[I] += temp_Transfer_functions.Emission_functions[I] * sin_pitch_angle * M_PI / Num_Samples_to_avg / 2;
            p_Transfer_functions->Emission_functions[Q] += temp_Transfer_functions.Emission_functions[Q] * sin_pitch_angle * M_PI / Num_Samples_to_avg / 2;
            p_Transfer_functions->Emission_functions[V] += temp_Transfer_functions.Emission_functions[V] * sin_pitch_angle * M_PI / Num_Samples_to_avg / 2;

            // The U component is 0 by definition
            p_Transfer_functions->Absorbtion_functions[I] += temp_Transfer_functions.Absorbtion_functions[I] * sin_pitch_angle * M_PI / Num_Samples_to_avg / 2;
            p_Transfer_functions->Absorbtion_functions[Q] += temp_Transfer_functions.Absorbtion_functions[Q] * sin_pitch_angle * M_PI / Num_Samples_to_avg / 2;
            p_Transfer_functions->Absorbtion_functions[V] += temp_Transfer_functions.Absorbtion_functions[V] * sin_pitch_angle * M_PI / Num_Samples_to_avg / 2;

            // The I and U components are 0 by definition
            p_Transfer_functions->Faradey_functions[Q] += temp_Transfer_functions.Faradey_functions[Q] * sin_pitch_angle * M_PI / Num_Samples_to_avg / 2;
            p_Transfer_functions->Faradey_functions[V] += temp_Transfer_functions.Faradey_functions[V] * sin_pitch_angle * M_PI / Num_Samples_to_avg / 2;

        }
    }
    else {

        /* The magnetic field is the one measured by a comoving with the plasma observer, but expressed in the cooridante frame */
        double pitch_angle = get_electron_pitch_angle(p_Emission_medium_state->Magnetic_fields.B_field_coord_frame, p_Emission_medium_state->Plasma_Velocity, State_Vector, p_Sim_Context);
        double sin_pitch_angle = sin(pitch_angle);

        double one_over_sqrt_sin = 1.0 / sqrt(sin_pitch_angle);
        double one_over_cbrt_sin = 1.0 / cbrt(sin_pitch_angle);
        double one_over_sin_to_0_p_5175 = 1.0 / pow(sin_pitch_angle, 0.5175);
        double one_over_sin_to_0_p_6    = 1.0 / pow(sin_pitch_angle, 0.6);
        double one_over_sin_to_0_p_7515 = 1.0 / pow(sin_pitch_angle, 0.7515);

        Transfer_args_corrected.X      = Transfer_args_uncorrected.X / sin_pitch_angle;
        Transfer_args_corrected.sqrt_X = Transfer_args_uncorrected.sqrt_X * one_over_sqrt_sin;
        Transfer_args_corrected.cbrt_X = Transfer_args_uncorrected.cbrt_X * one_over_cbrt_sin;
        Transfer_args_corrected.X_to_0_p_5175 = Transfer_args_uncorrected.X_to_0_p_5175 * one_over_sin_to_0_p_5175;
        Transfer_args_corrected.X_to_0_p_6    = Transfer_args_uncorrected.X_to_0_p_6 * one_over_sin_to_0_p_6;
        Transfer_args_corrected.X_to_0_p_7515 = Transfer_args_uncorrected.X_to_0_p_7515 * one_over_sin_to_0_p_7515;

        Transfer_args_corrected.sin_pitch_angle = sin_pitch_angle;
        Transfer_args_corrected.cos_pitch_angle = cos(pitch_angle);

        this->evaluate_synchrotron_transfer_functions(e_Thermal_ensamble, p_Emission_medium_state, &Transfer_args_corrected, p_Sim_Context, p_Transfer_functions);

    }

    /* Account for the relativistic doppler effet via the redshift. */
    for (int stokes_idx = 0; stokes_idx <= STOKES_PARAM_NUM - 1; stokes_idx++) {

        p_Transfer_functions->Emission_functions[stokes_idx] *= redshift * redshift;
        p_Transfer_functions->Faradey_functions[stokes_idx] /= redshift;
        p_Transfer_functions->Absorbtion_functions[stokes_idx] /= redshift;
    }

}
 
/* ========================================== Kappa synchrotron Transfer Functions ========================================== */

void Generic_Optically_Thin_Model::get_kappa_synchrotron_transfer_functions(const double* const State_Vector,
                                                                           const Simulation_Context_type* const p_Sim_Context,
                                                                           const Emission_medium_state_type* const p_Emission_medium_state,
                                                                           Transfer_functions_type* const p_Transfer_functions){

    /* === Zero out the transfer functions just in case === */
    memset(p_Transfer_functions, 0, sizeof(Transfer_functions_type));

    const double redshift = get_redshift(State_Vector, p_Emission_medium_state->Plasma_Velocity, p_Sim_Context->p_Observer);

    if (isinf(redshift) || isnan(redshift) || isinf(1.0 / redshift)) { return; }

    /* Dimensionless Electron Temperature */
    double T_electron_dim = BOLTZMANN_CONST_CGS * p_Emission_medium_state->Temperature / M_ELECTRON_CGS / C_LIGHT_CGS / C_LIGHT_CGS;

    /* Cyclotron Frequency */
    double f_cyclo = Q_ELECTRON_CGS * p_Emission_medium_state->Magnetic_fields.B_field_plasma_frame_norm / (2 * M_PI * M_ELECTRON_CGS * C_LIGHT_CGS);

    /* The "averaged" critical frequency (without the sin(theta) term - that gets added on later from a pre-computed table) */
    double f_k_no_sin = f_cyclo * (this->s_Emission_params.Kappa * T_electron_dim) * (this->s_Emission_params.Kappa * T_electron_dim);

    if (isinf(f_k_no_sin) || isnan(f_k_no_sin) || isinf(1.0 / f_k_no_sin)) { return; }

    /* Observation frequency */
    double& obs_frequency = p_Sim_Context->p_Init_Conditions->Observer_params.obs_frequency;

    Kappa_transfer_f_arguments_type Transfer_args_uncorrected{};

    /* Populate the angle-uncorrected transfer arguments struct. */
    Transfer_args_uncorrected.X              = obs_frequency / f_k_no_sin / redshift;
    Transfer_args_uncorrected.sqrt_X         = sqrt(Transfer_args_uncorrected.X);
    Transfer_args_uncorrected.cbrt_X         = cbrt(Transfer_args_uncorrected.X);
    Transfer_args_uncorrected.X_to_7_over_20 = pow(Transfer_args_uncorrected.X, 7. / 20);
    Transfer_args_uncorrected.kappa          = this->s_Emission_params.Kappa;
    Transfer_args_uncorrected.T_electron_dim = T_electron_dim;
    Transfer_args_uncorrected.frequency      = obs_frequency / redshift;

    int& Num_Samples_to_avg = p_Sim_Context->p_Init_Conditions->Emission_pitch_angle_samples_to_average;

    /* This structcs holds the transfer function args, corrected for the electron pitch angle. I am setting it equal to the uncorrected one
        so they copy over the variables that don't depend on the pitch angle. The rest get corrected inside the averaging loop. */
    Kappa_transfer_f_arguments_type Transfer_args_corrected = Transfer_args_uncorrected;

    if (p_Sim_Context->p_Init_Conditions->Average_electron_pitch_angle) {

        /* ============ This loop averages over the emission pitch angle, which it gets from a pre-computed table ============ */

        for (int averaging_idx = 1; averaging_idx <= Num_Samples_to_avg - 1; averaging_idx++) {

            double& sin_pitch_angle = this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[averaging_idx];

            Transfer_args_corrected.X                  = Transfer_args_uncorrected.X / sin_pitch_angle;
            Transfer_args_corrected.sqrt_X             = Transfer_args_uncorrected.sqrt_X * this->s_Precomputed_e_pitch_angles.one_over_sqrt_sin[averaging_idx];
            Transfer_args_corrected.cbrt_X             = Transfer_args_uncorrected.cbrt_X * this->s_Precomputed_e_pitch_angles.one_over_cbrt_sin[averaging_idx];
            Transfer_args_corrected.X_to_7_over_20     = Transfer_args_uncorrected.X_to_7_over_20 * this->s_Precomputed_e_pitch_angles.one_over_sin_to_7_over_20[averaging_idx];
            Transfer_args_corrected.sin_emission_angle = sin_pitch_angle;

            Transfer_functions_type temp_Transfer_functions{};

            this->evaluate_synchrotron_transfer_functions(e_Kappa_ensamble, p_Emission_medium_state, &Transfer_args_corrected, p_Sim_Context, &temp_Transfer_functions);

            // The U component is 0 by definition
            p_Transfer_functions->Emission_functions[I] += temp_Transfer_functions.Emission_functions[I] * sin_pitch_angle * M_PI / Num_Samples_to_avg / 2;
            p_Transfer_functions->Emission_functions[Q] += temp_Transfer_functions.Emission_functions[Q] * sin_pitch_angle * M_PI / Num_Samples_to_avg / 2;
            p_Transfer_functions->Emission_functions[V] += temp_Transfer_functions.Emission_functions[V] * sin_pitch_angle * M_PI / Num_Samples_to_avg / 2;

            // The U component is 0 by definition
            p_Transfer_functions->Absorbtion_functions[I] += temp_Transfer_functions.Absorbtion_functions[I] * sin_pitch_angle * M_PI / Num_Samples_to_avg / 2;
            p_Transfer_functions->Absorbtion_functions[Q] += temp_Transfer_functions.Absorbtion_functions[Q] * sin_pitch_angle * M_PI / Num_Samples_to_avg / 2;
            p_Transfer_functions->Absorbtion_functions[V] += temp_Transfer_functions.Absorbtion_functions[V] * sin_pitch_angle * M_PI / Num_Samples_to_avg / 2;

            // The I and U components are 0 by definition
            p_Transfer_functions->Faradey_functions[Q] += temp_Transfer_functions.Faradey_functions[Q] * sin_pitch_angle * M_PI / Num_Samples_to_avg / 2;
            p_Transfer_functions->Faradey_functions[V] += temp_Transfer_functions.Faradey_functions[V] * sin_pitch_angle * M_PI / Num_Samples_to_avg / 2;

        }
    }
    else {

        /* The magnetic field is the one measured by a comoving with the plasma observer, but expressed in the cooridante frame */
        double pitch_angle = get_electron_pitch_angle(p_Emission_medium_state->Magnetic_fields.B_field_coord_frame, p_Emission_medium_state->Plasma_Velocity, State_Vector, p_Sim_Context);
        double sin_pitch_angle = sin(pitch_angle);

        double one_over_sqrt_sin    = 1. / sqrt(sin_pitch_angle);
        double one_over_cbrt_sin    = 1. / cbrt(sin_pitch_angle);
        double one_over_7_to_20_sin = 1. / pow(sin_pitch_angle, 7. / 20);

        Transfer_args_corrected.X                  = Transfer_args_uncorrected.X / sin_pitch_angle;
        Transfer_args_corrected.sqrt_X             = Transfer_args_uncorrected.sqrt_X * one_over_sqrt_sin;
        Transfer_args_corrected.cbrt_X             = Transfer_args_uncorrected.cbrt_X * one_over_cbrt_sin;
        Transfer_args_corrected.X_to_7_over_20     = Transfer_args_uncorrected.X_to_7_over_20 * one_over_7_to_20_sin;
        Transfer_args_corrected.sin_emission_angle = sin_pitch_angle;

        this->evaluate_synchrotron_transfer_functions(e_Kappa_ensamble, p_Emission_medium_state, &Transfer_args_corrected, p_Sim_Context, p_Transfer_functions);

    }

    /* Account for the relativistic doppler effet via the redshift */
    for (int stokes_idx = 0; stokes_idx <= STOKES_PARAM_NUM - 1; stokes_idx++) {

        p_Transfer_functions->Emission_functions[stokes_idx] *= redshift * redshift;
        p_Transfer_functions->Faradey_functions[stokes_idx]  /= redshift;
        p_Transfer_functions->Absorbtion_functions[stokes_idx] /= redshift;
    }

}

/* ========================================== Phenomenological synchrotron Transfer Functions ========================================== */

void Generic_Optically_Thin_Model::get_phenomenological_synchrotron_functions(const double* const State_Vector,
                                                                              const Simulation_Context_type* const p_Sim_Context, 
                                                                              const Emission_medium_state_type* const p_Emission_medium_state,
                                                                              Transfer_functions_type* const p_Transfer_functions) {

    /* === Zero out the transfer functions just in case. === */
    memset(p_Transfer_functions, 0, sizeof(Transfer_functions_type));

    Phenomenological_transfer_f_arguments_type Transfer_args{};

    Transfer_args.redshift = get_redshift(State_Vector, p_Emission_medium_state->Plasma_Velocity, p_Sim_Context->p_Observer);

    if (isinf(Transfer_args.redshift) || isnan(Transfer_args.redshift) || isinf(1.0 / Transfer_args.redshift)) { return; }

    Transfer_args.frequency = p_Sim_Context->p_Init_Conditions->Observer_params.obs_frequency / Transfer_args.redshift;
    Transfer_args.f_cyclo = Q_ELECTRON_CGS * p_Emission_medium_state->Magnetic_fields.B_field_plasma_frame_norm / (2 * M_PI * M_ELECTRON_CGS * C_LIGHT_CGS);

    this->evaluate_synchrotron_transfer_functions(e_Phenomenological_ensamble, p_Emission_medium_state, &Transfer_args, p_Sim_Context, p_Transfer_functions);

    /* Account for the relativistic doppler effet via the redshift. */
    for (int stokes_idx = 0; stokes_idx <= STOKES_PARAM_NUM - 1; stokes_idx++) {

        p_Transfer_functions->Emission_functions[stokes_idx] *= Transfer_args.redshift * Transfer_args.redshift;
        p_Transfer_functions->Faradey_functions[stokes_idx] /= Transfer_args.redshift;
        p_Transfer_functions->Absorbtion_functions[stokes_idx] /= Transfer_args.redshift;
    }
}

/* ============================================ Main "Selector" For The Transfer Functions ============================================ */

void Generic_Optically_Thin_Model::get_radiative_transfer_functions(const double* const State_Vector,
                                                                    const Simulation_Context_type* const p_Sim_Context,
                                                                    const Emission_medium_enums Emission_medium,
                                                                    Transfer_functions_type* const p_Transfer_functions) {

    /* === Zero out the transfer functions just in case. === */
    memset(p_Transfer_functions, 0, sizeof(Transfer_functions_type));

    Emission_medium_state_type Emission_medium_state{};

    switch (Emission_medium) {

    case Disk:

        Emission_medium_state.Density         = this->get_disk_density(State_Vector);
        Emission_medium_state.Temperature     = this->get_disk_temperature(State_Vector);
        Emission_medium_state.Ensamble_type   = p_Sim_Context->p_Init_Conditions->Disk_params.Ensamble_type;
        Emission_medium_state.Plasma_Velocity = this->get_plasma_velocity(State_Vector, p_Sim_Context, this->s_Disk_params.Velocity_profile_type, this->s_Disk_params.Radial_velocity_fraction);
        Emission_medium_state.Magnetization   = p_Sim_Context->p_Init_Conditions->Disk_params.Magnetization;

        memcpy(Emission_medium_state.Magnetic_fields.Magnetic_field_geometry, p_Sim_Context->p_Init_Conditions->Disk_params.Mag_field_geometry, 3 * sizeof(double));

        break;

    case Hotspot:

        Emission_medium_state.Density         = this->get_hotspot_density(State_Vector, p_Sim_Context);
        Emission_medium_state.Temperature     = this->get_hotspot_temperature(State_Vector, p_Sim_Context);
        Emission_medium_state.Ensamble_type   = p_Sim_Context->p_Init_Conditions->Hotspot_params.Ensamble_type;
        Emission_medium_state.Plasma_Velocity = this->get_plasma_velocity(State_Vector, p_Sim_Context, this->s_Hotspot_params.Velocity_profile_type, this->s_Hotspot_params.Radial_velocity_fraction);
        Emission_medium_state.Magnetization   = p_Sim_Context->p_Init_Conditions->Hotspot_params.Magnetization;

        memcpy(Emission_medium_state.Magnetic_fields.Magnetic_field_geometry, p_Sim_Context->p_Init_Conditions->Hotspot_params.Mag_field_geometry, 3 * sizeof(double));

        break;

    default:

        std::cout << "Unsupported emissison medium - something broke in the get_radiative_transfer_functions function!" << "\n";

        exit(ERROR);

        break;

    }

    if (NULL == Emission_medium_state.Plasma_Velocity) { return; }

    this->get_magnetic_field(State_Vector, p_Sim_Context, &Emission_medium_state);

    switch (Emission_medium_state.Ensamble_type) {

    case(e_Phenomenological_ensamble):

        this->get_phenomenological_synchrotron_functions(State_Vector, p_Sim_Context, &Emission_medium_state, p_Transfer_functions);
        break;

    case(e_Kappa_ensamble):

        this->get_kappa_synchrotron_transfer_functions(State_Vector, p_Sim_Context, &Emission_medium_state,  p_Transfer_functions);
        break;

    default:

        this->get_thermal_synchrotron_transfer_functions(State_Vector, p_Sim_Context, &Emission_medium_state, p_Transfer_functions);
        break;
    }
}

void Generic_Optically_Thin_Model::evaluate_synchrotron_transfer_functions(const Ensamble_enums e_Ensamble_type,
                                                                           const Emission_medium_state_type* const p_Emission_medium_state,
                                                                           const void* const p_Transfer_args,
                                                                           const Simulation_Context_type* const p_Sim_Context,
                                                                           Transfer_functions_type* const p_Transfer_functions) {

    /* The dimensionless frequency needs to get extracted from the p_Transfer_args pointer, but it first needs to be recast to not-void.
       This happens in the scopes of the switch statemeent below, so I create a variable here to store it. */
    double  frequency_dim{};
    double& obs_frequency = p_Sim_Context->p_Init_Conditions->Observer_params.obs_frequency;

    switch (e_Ensamble_type) {

    case e_Thermal_ensamble:

        this->get_thermal_synchrotron_emission_fit_functions(static_cast<const Thermal_transfer_f_arguments_type*>(p_Transfer_args), p_Transfer_functions->Emission_functions);
        this->get_thermal_synchrotron_absorbtion_fit_functions(static_cast<const Thermal_transfer_f_arguments_type*>(p_Transfer_args), p_Emission_medium_state, p_Transfer_functions->Emission_functions, p_Transfer_functions->Absorbtion_functions);
        this->get_thermal_synchrotron_faradey_fit_functions(static_cast<const Thermal_transfer_f_arguments_type*>(p_Transfer_args), p_Transfer_functions->Faradey_functions);

        frequency_dim = static_cast<const Thermal_transfer_f_arguments_type*>(p_Transfer_args)->frequency / obs_frequency;

        break;

    case e_Kappa_ensamble:

        this->get_kappa_synchrotron_emission_fit_functions(static_cast<const Kappa_transfer_f_arguments_type*>(p_Transfer_args), p_Transfer_functions->Emission_functions);
        this->get_kappa_synchrotron_absorbtion_fit_functions(static_cast<const Kappa_transfer_f_arguments_type*>(p_Transfer_args), p_Transfer_functions->Absorbtion_functions);

        // TODO: Add the Faradey function fits.

        frequency_dim = static_cast<const Kappa_transfer_f_arguments_type*>(p_Transfer_args)->frequency / obs_frequency;

        break;

    case e_Phenomenological_ensamble:

        this->get_phenomenological_synchrotron_fit_functions(static_cast<const Phenomenological_transfer_f_arguments_type*>(p_Transfer_args), p_Transfer_functions);

        frequency_dim = static_cast<const Phenomenological_transfer_f_arguments_type*>(p_Transfer_args)->frequency / obs_frequency;

        break;

    default:

        std::cout << "\n" << "Error! Unsupported emission model - something broke in the evaluate_synchrotron_transfer_functions function!" << "\n";
        exit(ERROR);

    }

    /* The below coefficients pop up in the dimentionless radiative transfer equation. */

    const double f_cyclo_dim    = Q_ELECTRON_CGS * p_Emission_medium_state->Magnetic_fields.B_field_plasma_frame_norm / (2 * M_PI * M_ELECTRON_CGS * C_LIGHT_CGS) / obs_frequency;
    const double distance_scale = p_Sim_Context->p_Init_Conditions->central_object_mass * M_SUN_SI * G_NEWTON_SI / C_LIGHT_SI / C_LIGHT_SI * METER_TO_CM;
    const double transport_matrix_ratio = Global_density_scale * Q_ELECTRON_CGS * Q_ELECTRON_CGS * distance_scale / C_LIGHT_CGS / obs_frequency / M_ELECTRON_CGS;

    /* ================================================ The emission functions ================================================ */

    p_Transfer_functions->Emission_functions[I] *= (p_Emission_medium_state->Density / Global_density_scale) * f_cyclo_dim;
    p_Transfer_functions->Emission_functions[Q] *= (p_Emission_medium_state->Density / Global_density_scale) * f_cyclo_dim;
    p_Transfer_functions->Emission_functions[V] *= (p_Emission_medium_state->Density / Global_density_scale) * f_cyclo_dim;

    /* ================================================ The absorbtion functions ================================================ */

    p_Transfer_functions->Absorbtion_functions[I] *= (p_Emission_medium_state->Density / Global_density_scale) * transport_matrix_ratio / frequency_dim;
    p_Transfer_functions->Absorbtion_functions[Q] *= (p_Emission_medium_state->Density / Global_density_scale) * transport_matrix_ratio / frequency_dim;
    p_Transfer_functions->Absorbtion_functions[V] *= (p_Emission_medium_state->Density / Global_density_scale) * transport_matrix_ratio / frequency_dim;

    /* ================================================ The faradey functions ================================================ */
    /* Originally derived in https://iopscience.iop.org/article/10.1086/592326/pdf - expressions 25, 26 and 33. */

    p_Transfer_functions->Faradey_functions[Q] *= -(p_Emission_medium_state->Density / Global_density_scale) * transport_matrix_ratio * (f_cyclo_dim / frequency_dim) * (f_cyclo_dim / frequency_dim) / frequency_dim;
    p_Transfer_functions->Faradey_functions[V] *= 2 * (p_Emission_medium_state->Density / Global_density_scale) * transport_matrix_ratio * (f_cyclo_dim / frequency_dim) / frequency_dim;

}

/* ========================================================== Misc Functions ========================================================== */

void Generic_Optically_Thin_Model::precompute_electron_pitch_angles(Initial_conditions_type* p_Init_Conditions) {

    // ====================================================== Allocate memory for the arrays ====================================================== //

    this->s_Precomputed_e_pitch_angles.cos_electron_pitch_angles = new double[p_Init_Conditions->Emission_pitch_angle_samples_to_average];
    this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles = new double[p_Init_Conditions->Emission_pitch_angle_samples_to_average];

    this->s_Precomputed_e_pitch_angles.one_over_sqrt_sin = new double[p_Init_Conditions->Emission_pitch_angle_samples_to_average];
    this->s_Precomputed_e_pitch_angles.one_over_cbrt_sin = new double[p_Init_Conditions->Emission_pitch_angle_samples_to_average];

    this->s_Precomputed_e_pitch_angles.one_over_sin_to_0_p_5175 = new double[p_Init_Conditions->Emission_pitch_angle_samples_to_average];
    this->s_Precomputed_e_pitch_angles.one_over_sin_to_0_p_6    = new double[p_Init_Conditions->Emission_pitch_angle_samples_to_average];
    this->s_Precomputed_e_pitch_angles.one_over_sin_to_0_p_7515 = new double[p_Init_Conditions->Emission_pitch_angle_samples_to_average];

    this->s_Precomputed_e_pitch_angles.one_over_sin_to_7_over_20 = new double[p_Init_Conditions->Emission_pitch_angle_samples_to_average];

    // =========================================================================================================================================== //

    for (int index = 0; index <= p_Init_Conditions->Emission_pitch_angle_samples_to_average - 1; index++) {

        double pitch_angle = double(index) / p_Init_Conditions->Emission_pitch_angle_samples_to_average * M_PI;
        this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index] = sin(pitch_angle);
        this->s_Precomputed_e_pitch_angles.cos_electron_pitch_angles[index] = cos(pitch_angle);

        if (this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index] != 0) {

            // Used in the thermal and kappa synchrotron emission functions.

            this->s_Precomputed_e_pitch_angles.one_over_sqrt_sin[index] = 1. / sqrt(this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index]);
            this->s_Precomputed_e_pitch_angles.one_over_cbrt_sin[index] = 1. / cbrt(this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index]);

            // Used in the thermal synchrotron Faradey functions.

            this->s_Precomputed_e_pitch_angles.one_over_sin_to_0_p_5175[index] = 1. / pow(this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index], 0.5175);
            this->s_Precomputed_e_pitch_angles.one_over_sin_to_0_p_6[index]    = 1. / pow(this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index], 0.6);
            this->s_Precomputed_e_pitch_angles.one_over_sin_to_0_p_7515[index] = 1. / pow(this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index], 0.7515);

            // Used in the kappa synchrotron emission functions.

            this->s_Precomputed_e_pitch_angles.one_over_sin_to_7_over_20[index] = 1. / pow(this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index], 7. / 20);

        }
    }
}

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
