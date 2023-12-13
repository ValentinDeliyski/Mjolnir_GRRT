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

        r_in = Spacetimes[e_metric]->get_ISCO()[Inner];

    }
    
    r_out = y;

};

double Novikov_Thorne_Model::Keplerian_angular_velocity(double r, Spacetime_Base_Class* Spacetimes[]) {

    double State_Vector[2] = { r, M_PI_2 };

    Spacetimes[e_metric]->set_ignore_flag(true);

    Metric_type s_dr_Metric = Spacetimes[e_metric]->get_dr_metric(State_Vector);

    Spacetimes[e_metric]->set_ignore_flag(false);

    return (-s_dr_Metric.Metric[0][3] + sqrt(s_dr_Metric.Metric[0][3] * s_dr_Metric.Metric[0][3] - s_dr_Metric.Metric[0][0] * s_dr_Metric.Metric[3][3])) / s_dr_Metric.Metric[3][3];

}

double Novikov_Thorne_Model::dr_Keplerian_angular_velocity(double r, Spacetime_Base_Class* Spacetimes[]) {

    double State_Vector[2] = { r, M_PI_2 };

    Spacetimes[e_metric]->set_ignore_flag(true);

    Metric_type s_dr_Metric  = Spacetimes[e_metric]->get_dr_metric(State_Vector);
    Metric_type s_d2r_Metric = Spacetimes[e_metric]->get_d2r_metric(State_Vector);

    double root = sqrt(s_dr_Metric.Metric[0][3] * s_dr_Metric.Metric[0][3] - s_dr_Metric.Metric[0][0] * s_dr_Metric.Metric[3][3]);

    double Kepler = this->Keplerian_angular_velocity(r, Spacetimes);

    Spacetimes[e_metric]->set_ignore_flag(false);

    return  -Kepler / s_dr_Metric.Metric[3][3] * s_d2r_Metric.Metric[3][3] + (-s_d2r_Metric.Metric[0][3] 
            + 1.0 / root / 2 * (2 * s_dr_Metric.Metric[0][3] * s_d2r_Metric.Metric[0][3] - s_dr_Metric.Metric[0][0] * s_d2r_Metric.Metric[3][3]
                                - s_d2r_Metric.Metric[0][0] * s_dr_Metric.Metric[3][3])) / s_dr_Metric.Metric[3][3];

}

double Novikov_Thorne_Model::Redshift(double J, double State_Vector[], double r_obs, double theta_obs,
    Spacetime_Base_Class* Spacetimes[]) {

    double& r_source     = State_Vector[e_r];
    double& theta_source = State_Vector[e_theta];

    Spacetimes[e_metric]->set_ignore_flag(true);

    /*
    Get the observer 4-velocity
    */

    double State_Vector_obs[2] = {r_obs, theta_obs};

    Metric_type s_Metric_obs = Spacetimes[e_metric]->get_metric(State_Vector_obs);

    double U_obs[4] = { 1.0 / s_Metric_obs.Lapse_function, 0 ,0 , s_Metric_obs.Shift_function / s_Metric_obs.Lapse_function };

    /*
    Get the source 4-velocity
    */

    Metric_type s_Metric_source = Spacetimes[e_metric]->get_metric(State_Vector);

    double Kepler = this->Keplerian_angular_velocity(r_source, Spacetimes);

    double gamma = 1 / sqrt(-s_Metric_source.Metric[0][0] - 2 * s_Metric_source.Metric[0][3] * Kepler - s_Metric_source.Metric[3][3] * Kepler * Kepler);

    if (gamma != gamma) {

        return 0;
    }

    double U_source[4] = { gamma, 0, 0, gamma * Kepler };

    Spacetimes[e_metric]->set_ignore_flag(false);

    return  (-U_obs[0] + U_obs[3] * J) / (-U_source[0] + U_source[3] * J);

}

double Novikov_Thorne_Model::disk_Energy(double r, Spacetime_Base_Class* Spacetimes[]) {

    double State_Vector[2] = { r_obs, M_PI_2 };

    Spacetimes[e_metric]->set_ignore_flag(true);

    Metric_type s_Metric_source = Spacetimes[e_metric]->get_metric(State_Vector);

    double Kepler = this->Keplerian_angular_velocity(r, Spacetimes);

    double root = sqrt(-s_Metric_source.Metric[0][0] - 2 * s_Metric_source.Metric[0][3] * Kepler - s_Metric_source.Metric[3][3] * Kepler * Kepler);

    Spacetimes[e_metric]->set_ignore_flag(false);

    return  -(s_Metric_source.Metric[0][0] + s_Metric_source.Metric[0][3] * Kepler) / root;

}

double Novikov_Thorne_Model::disk_Angular_Momentum(double r, Spacetime_Base_Class* Spacetimes[]) {

    double State_Vector[2] = { r_obs, M_PI_2 };
    
    Spacetimes[e_metric]->set_ignore_flag(true);

    Metric_type s_Metric_source = Spacetimes[e_metric]->get_metric(State_Vector);

    double Kepler = this->Keplerian_angular_velocity(r, Spacetimes);

    double root = sqrt(-s_Metric_source.Metric[0][0] - 2 * s_Metric_source.Metric[0][3] * Kepler - s_Metric_source.Metric[3][3] * Kepler * Kepler);

    Spacetimes[e_metric]->set_ignore_flag(false);

    return  (s_Metric_source.Metric[3][3] * Kepler + s_Metric_source.Metric[0][3]) / root;

}

double Novikov_Thorne_Model::Flux_integrand(double r, Spacetime_Base_Class* Spacetimes[]) {

    double State_Vector[2] = { r, M_PI_2 };

    Spacetimes[e_metric]->set_ignore_flag(true);

    Metric_type s_Metric    = Spacetimes[e_metric]->get_metric(State_Vector);
    Metric_type s_dr_Metric = Spacetimes[e_metric]->get_dr_metric(State_Vector);

    double Kepler    = this->Keplerian_angular_velocity(r, Spacetimes);
    double dr_Kepler = this->dr_Keplerian_angular_velocity(r, Spacetimes);

    double metric_det = get_metric_det(s_Metric.Metric);

    double root = sqrt(-s_Metric.Metric[0][0] - 2 * s_Metric.Metric[0][3] * Kepler - s_Metric.Metric[3][3] * Kepler * Kepler);
    double dr_root = (-s_dr_Metric.Metric[0][0] - 2 * (s_dr_Metric.Metric[0][3] * Kepler + s_Metric.Metric[0][3] * dr_Kepler) 
                   - s_dr_Metric.Metric[3][3] * Kepler * Kepler - 2 * s_Metric.Metric[3][3] * Kepler * dr_Kepler);

    double E = this->disk_Energy(r, Spacetimes);
    double L = this->disk_Angular_Momentum(r, Spacetimes);

    double dr_L = (s_dr_Metric.Metric[3][3] * Kepler + s_Metric.Metric[3][3] * dr_Kepler + s_dr_Metric.Metric[0][3]) / root - L / root / root / 2 * dr_root;

    Spacetimes[e_metric]->set_ignore_flag(false);

    return (E - Kepler * L) * dr_L;

}

double Novikov_Thorne_Model::solve_Flux_integral(double lower_bound, double upper_bound, double tolerance, Spacetime_Base_Class* Spacetimes[]) {

    Spacetimes[e_metric]->set_ignore_flag(true);

    double mid_point       = (lower_bound + upper_bound) / 2;
    double left_mid_point  = (lower_bound + mid_point) / 2;
    double right_mid_point = (mid_point + upper_bound) / 2;

    double F_lower_bound = this->Flux_integrand(lower_bound, Spacetimes);
    double F_mid_point   = this->Flux_integrand(mid_point, Spacetimes);
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

    Spacetimes[e_metric]->set_ignore_flag(false);

    return integral;
}

double Novikov_Thorne_Model::get_flux(double r, Spacetime_Base_Class* Spacetimes[]) {

    double State_Vector[2] = { r, M_PI_2 };

    Spacetimes[e_metric]->set_ignore_flag(true);

    Metric_type s_Metric = Spacetimes[e_metric]->get_metric(State_Vector);

    double metric_det = get_metric_det(s_Metric.Metric);
    double E_disk     = disk_Energy(r, Spacetimes);
    double L_disk     = disk_Angular_Momentum(r, Spacetimes);

    double Kepler    =    Keplerian_angular_velocity(r, Spacetimes);
    double dr_Kepler = dr_Keplerian_angular_velocity(r, Spacetimes);

    double Flux_coeff = -dr_Kepler / ((E_disk - Kepler * L_disk) * (E_disk - Kepler * L_disk)) / (4 * M_PI * sqrt(metric_det));

    if (e_metric == Wormhole) {

        r = sqrt(r * r + WH_R_THROAT * WH_R_THROAT);

    }

    double Flux_integral = solve_Flux_integral(r_in, r, INTEGRAL_ACCURACY, Spacetimes);

    Spacetimes[e_metric]->set_ignore_flag(false);

    return Flux_coeff * Flux_integral;

}

double Novikov_Thorne_Model::get_r_in()  { return r_in; };
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

        return ERROR;
    }

    if (NULL != p_Emission_Parameters) {

        this->s_Emission_Parameters = *p_Emission_Parameters;

    }
    else {

        std::cout << "p_Emission_Parameters is a NULL pointer! \n";

        return ERROR;
    }

    return OK;

}

Disk_model_parameters Optically_Thin_Toroidal_Model::get_disk_params() { return this->s_Disk_Parameters; };

/* Disk Temperature Functions */

int Optically_Thin_Toroidal_Model::update_disk_temperature(double State_vector[]) {

    double& r = State_vector[e_r];
    double Radial_Cutoff{};

    double T = T_ELECTRON_EXACT_CGS * this->s_Disk_Parameters.Power_law_radial_scale / r;

    if (r < this->s_Disk_Parameters.Disk_r_cutoff){

        Radial_Cutoff = (r - this->s_Disk_Parameters.Disk_r_cutoff) / this->s_Disk_Parameters.Disk_cutoff_scale;

        T *= exp(-Radial_Cutoff * Radial_Cutoff);

    }

    this->Disk_Temperature = T;

    if (T < 0) {
    
        return ERROR;

    }

    return OK;

}

double Optically_Thin_Toroidal_Model::get_disk_temperature(double State_vector[]) {

    int result = this->update_disk_temperature(State_vector);

    if (ERROR == result) {

        std::cout << "Invalid Disk Temperature: " << this->Disk_Temperature << "\n";

        return NULL;

    }

    return this->Disk_Temperature;

}

/* Disk Velocity Functions */

int Optically_Thin_Toroidal_Model::update_disk_velocity(double State_Vector[], Initial_conditions_type* s_Initial_Conditions) {

    double& r_source     = State_Vector[e_r];
    double& theta_source = State_Vector[e_theta];

    Metric_type s_Metric = s_Initial_Conditions->Spacetimes[e_metric]->get_metric(State_Vector);

    double rho = r_source * sin(theta_source);

    if (rho < 0.0) {

        rho *= -1.0;

    }

    double sqrt_rho = sqrt(rho);
    double ell      = sqrt_rho * sqrt_rho * sqrt_rho / (1 + rho);

    if (e_metric == Naked_Singularity) {

        ell *= pow(1 - JNW_R_SINGULARITY / r_source, JNW_GAMMA);

    }

    double u_t{}, u_phi{};

    double inv_metric[4][4]{};

    invert_metric(inv_metric, s_Metric.Metric);

    u_t = -1.0 / sqrt(-(inv_metric[0][0] - 2 * inv_metric[0][3] * ell + inv_metric[3][3] * ell * ell));
    u_phi = -u_t * ell;

    /*
    Convert U_source to contravariant components
    */

    this->Disk_velocity[0] = inv_metric[0][0] * u_t + inv_metric[0][3] * u_phi;
    this->Disk_velocity[1] = 0;
    this->Disk_velocity[2] = 0;
    this->Disk_velocity[3] = inv_metric[3][3] * u_phi + inv_metric[3][0] * u_t;

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

        return NULL;

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

    double hotspot_density  = exp(-(x_center - x_photon) * (x_center - x_photon) / HOTSPOT_SCALE / HOTSPOT_SCALE);
           hotspot_density *= exp(-(y_center - y_photon) * (y_center - y_photon) / HOTSPOT_SCALE / HOTSPOT_SCALE);
           hotspot_density *= exp(-(z_center - z_photon) * (z_center - z_photon) / HOTSPOT_SCALE / HOTSPOT_SCALE);
           hotspot_density *= HOTSPOT_REL_SCALE;

    this->Disk_hotspot_density = hotspot_density;

    if (Disk_hotspot_density < 0) {

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

        return ERROR;

    }

}

int Optically_Thin_Toroidal_Model::update_disk_density_profile(double State_Vector[]) {

    double& r  = State_Vector[e_r];
    double rho = sin(State_Vector[e_theta]);
    double h   = cos(State_Vector[e_theta]);

    double Height_Cutoff{};
    double Radial_Cutoff{}; 

    double electron_density_profile{};
    
    switch (e_disk_model){

        case Power_law:

            Height_Cutoff = h / (this->s_Disk_Parameters.Disk_opening_angle * rho);
            electron_density_profile = exp(-Height_Cutoff * Height_Cutoff / 2) / (r / R_0) / (r / R_0);

           if (r < this->s_Disk_Parameters.Disk_r_cutoff){
                
               Radial_Cutoff = (r - this->s_Disk_Parameters.Disk_r_cutoff) / this->s_Disk_Parameters.Disk_cutoff_scale;
               electron_density_profile *= exp(-Radial_Cutoff * Radial_Cutoff);

           }

            break;

        case Exponential_law:

            Height_Cutoff = h / this->s_Disk_Parameters.Exp_law_height_scale;
            Radial_Cutoff = r / this->s_Disk_Parameters.Power_law_radial_scale;

            electron_density_profile = exp(-Radial_Cutoff * Radial_Cutoff / 2 - Height_Cutoff * Height_Cutoff / 2);

            break;

        default:

            std::cout << "Wrong disk density model!" << '\n';

    }

    if (this->s_Disk_Parameters.Hotspot_scale != 0.0f) {

        this->Disk_density_profile = electron_density_profile + get_disk_hotspot(State_Vector);

    }
    else {

        this->Disk_density_profile =  electron_density_profile;

    }

    if (this->Disk_density_profile < 0) {

        return ERROR;

    }

    return OK;

}

double Optically_Thin_Toroidal_Model::get_disk_density_profile(double State_Vector[]) {

    int result = this->update_disk_density_profile(State_Vector);

    if (ERROR == result) {

        return ERROR;

    }

    return this->Disk_density_profile;

}

/* Disk Magnetic Field Functions */

double Optically_Thin_Toroidal_Model::get_magnetic_field(double B_field[3], double State_vector[]) {

    /*

    Everything is in GCS!

    */

    double electron_density = this->s_Disk_Parameters.Density_scale * this->get_disk_density_profile(State_vector);

    double B_CGS = sqrt(DISK_MAGNETIZATION * C_LIGHT_CGS * C_LIGHT_CGS * electron_density * M_PROTON_CGS * 4 * M_PI);

    B_field[x] = B_CGS * MAG_FIELD_GEOMETRY[x];
    B_field[y] = B_CGS * MAG_FIELD_GEOMETRY[y];
    B_field[z] = B_CGS * MAG_FIELD_GEOMETRY[z];

    return B_CGS;

}

double Optically_Thin_Toroidal_Model::get_electron_pitch_angle(double State_Vector[], double B_field_local[], Initial_conditions_type* s_Initial_Conditions) {

    double* U_source_coord = get_disk_velocity(State_Vector, s_Initial_Conditions);

    /*

    Transform U_source To The ZAMO Frame

    */

    Metric_type s_Metric = s_Initial_Conditions->Spacetimes[e_metric]->get_metric(State_Vector);

    double U_source_ZAMO[4]{};
    Contravariant_coord_to_ZAMO(s_Metric.Metric, U_source_coord, U_source_ZAMO);

    /*

    Boost U_source_ZAMO To The Fluid Frame

    */

    double Boost_matrix[4][4]{};

    Lorentz_boost_matrix(Boost_matrix, U_source_ZAMO);

    double U_source_Boosted[4]{};
    
    mat_vec_multiply_4D(Boost_matrix, U_source_ZAMO, U_source_Boosted);

    double  U_source_local[3] = { U_source_Boosted[1] / U_source_Boosted[0], U_source_Boosted[2] / U_source_Boosted[0], U_source_Boosted[3] / U_source_Boosted[0] };

    /*

    Get Sin Of The Angle Between The Local Fluid Velocity And Magnetic Field In The ZAMO Frame

    */

    double cos_angle = 1, sin_angle{};

    if (vector_norm(B_field_local, 3) > 1e-10 && dot_product(U_source_local, U_source_local) != 0) {

        cos_angle = dot_product(B_field_local, U_source_local) / sqrt(dot_product(B_field_local, B_field_local) * dot_product(U_source_local, U_source_local));
        sin_angle = sqrt(1 - cos_angle * cos_angle);

    }

    return sin_angle;

}

/* Disk Emission Functions */

double Optically_Thin_Toroidal_Model::get_emission_function_synchotron_exact(double State_vector[], Initial_conditions_type* s_Initial_conditions) {

    /* Electron Density in CGS */

    double electron_density = this->s_Disk_Parameters.Density_scale * this->get_disk_density_profile(State_vector);

    /* Dimentionless Electron Temperature */

    double T_electron     = this->get_disk_temperature(State_vector);
    double T_electron_dim = BOLTZMANN_CONST_CGS * T_electron / M_ELECTRON_CGS / C_LIGHT_CGS / C_LIGHT_CGS;
    double divisor        = std::cyl_bessel_k(2.0, 1.0 / T_electron_dim);

    /* Magnetic Field */

    double B_field_local[3];
    double B_CGS = this->get_magnetic_field(B_field_local, State_vector);

    /* Disk Coordinate Velocity */

    double* U_source_coord = this->get_disk_velocity(State_vector, s_Initial_conditions);

    /* Redshit */

    double redshift = Redshift(State_vector, U_source_coord);

    /* Cyclotron Frequency */

    double f_cyclo = Q_ELECTRON_CGS * B_CGS / (2 * M_PI * M_ELECTRON_CGS * C_LIGHT_CGS);

    /* Synchotron Frequency (without the sin(theta) term - that gets added on later from a pre-computed table) */

    double f_s_no_sin = 2. / 9 * f_cyclo * T_electron_dim * T_electron_dim;

    double X = 1e100;

    if (f_s_no_sin != 0) {

        X = OBS_FREQUENCY_CGS / f_s_no_sin / redshift;

    }

    /* The exact form of the emission function is from "Images and photon ring signatures of thick disks around black holes" */

    double X_1_2 = sqrt(X);
    double X_1_3 = cbrt(X);
    double X_1_6 = sqrt(X_1_3);

    /* Average the Emission Function over all electron pitch angles */

    double constant_coeff = M_SQRT2 * M_PI * Q_ELECTRON_CGS * Q_ELECTRON_CGS / 3 / C_LIGHT_CGS;
    double emission_function{};

    for (int averaging_idx = 0; averaging_idx <= NUM_SAMPLES_TO_AVG - 1; averaging_idx++) {

        double sin_pitch_angle = this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[averaging_idx];

        double X_1_2_angle_corrected = X_1_2 * this->s_Precomputed_e_pitch_angles.one_over_sqrt_sin[averaging_idx];
        double X_1_3_angle_corrected = X_1_3 * this->s_Precomputed_e_pitch_angles.one_over_cbrt_sin[averaging_idx];
        double X_1_6_angle_corrected = X_1_6 * this->s_Precomputed_e_pitch_angles.one_over_sin_to_1_6[averaging_idx];

        double X_term = X_1_2_angle_corrected + 1.887749 * X_1_6_angle_corrected; // The constant is = pow(2, 11.0 / 12.0)
        X_term *= X_term;

        if (divisor != 0.) {

            emission_function += constant_coeff * electron_density * (f_s_no_sin * sin_pitch_angle) * X_term * exp(-X_1_3_angle_corrected) / divisor * sin_pitch_angle * M_PI / NUM_SAMPLES_TO_AVG;

        }

    }

    /* Factor of 1 / 2 from the averaging */

    return emission_function / 2;

}

double Optically_Thin_Toroidal_Model::get_emission_function_synchotron_phenomenological(double State_vector[], Initial_conditions_type* s_Initial_conditions) {

    double& emission_power_law = this->s_Emission_Parameters.Emission_power_law;
    double& emission_scale     = this->s_Emission_Parameters.Emission_scale;

    double electron_density = emission_scale * get_disk_density_profile(State_vector);

    /* Disk Coordinate Velocity */

    double* U_source_coord = this->get_disk_velocity(State_vector, s_Initial_conditions);

    double redshift = Redshift(State_vector, U_source_coord);

    return electron_density * pow(redshift, emission_power_law);

}

double Optically_Thin_Toroidal_Model::get_absorbtion_function(double Emission_Function, double State_vector[], double redshift, double Frequency) {

    double Planck_function_CGS{};

    double& abs_coeff          = this->s_Emission_Parameters.Absorbtion_coeff;
    double& emission_scale     = this->s_Emission_Parameters.Emission_scale;
    double& source_f_power_law = this->s_Emission_Parameters.Source_f_power_law;
    double& emission_power_law = this->s_Emission_Parameters.Emission_power_law;

    switch (e_emission) {

    case Synchotron_exact: 

        Planck_function_CGS = get_planck_function_CGS(Frequency, this->get_disk_temperature(State_vector));

        if (Planck_function_CGS != 0.) {

            return Emission_Function / Planck_function_CGS;

        }
        else {

            return 0;

        }

    case Synchotron_phenomenological:

        return abs_coeff * emission_scale * this->get_disk_density_profile(State_vector) * pow(redshift, source_f_power_law + emission_power_law);

    default:

        std::cout << "Wrong emission model!" << '\n'; 

        return ERROR; 

    }
}

void Optically_Thin_Toroidal_Model::precompute_electron_pitch_angles() {

    for (int index = 0; index <= NUM_SAMPLES_TO_AVG - 1; index++) {

        double pitch_angle = double(index) / NUM_SAMPLES_TO_AVG * M_PI;
        this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index] = sin(pitch_angle);

        if (this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index] != 0) {

            this->s_Precomputed_e_pitch_angles.one_over_sqrt_sin[index]   = 1. / sqrt(this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index]);
            this->s_Precomputed_e_pitch_angles.one_over_cbrt_sin[index]   = 1. / cbrt(this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index]);
            this->s_Precomputed_e_pitch_angles.one_over_sin_to_1_6[index] = 1. / sqrt(this->s_Precomputed_e_pitch_angles.one_over_cbrt_sin[index]);

        }

    }

}

double get_planck_function_CGS(double Frequency, double Temperature) {

    return 2 * PLANCK_CONSTANT_CGS * Frequency * Frequency * Frequency / C_LIGHT_CGS / C_LIGHT_CGS / (exp(PLANCK_CONSTANT_CGS * Frequency / BOLTZMANN_CONST_CGS / Temperature) - 1);

}

