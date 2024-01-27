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

    double State_Vector[2] = { r, M_PI_2 };

    Spacetimes[e_metric]->set_ignore_flag(true);

    Metric_type s_Metric_source = Spacetimes[e_metric]->get_metric(State_Vector);

    double Kepler = this->Keplerian_angular_velocity(r, Spacetimes);

    double root = sqrt(-s_Metric_source.Metric[0][0] - 2 * s_Metric_source.Metric[0][3] * Kepler - s_Metric_source.Metric[3][3] * Kepler * Kepler);

    Spacetimes[e_metric]->set_ignore_flag(false);

    return  -(s_Metric_source.Metric[0][0] + s_Metric_source.Metric[0][3] * Kepler) / root;

}

double Novikov_Thorne_Model::disk_Angular_Momentum(double r, Spacetime_Base_Class* Spacetimes[]) {

    double State_Vector[2] = { r, M_PI_2 };
    
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

/* =============================================== Disk Transfer Functions Functions =============================================== */

void Optically_Thin_Toroidal_Model::get_faradey_functions(double State_vector[], 
                                                          double X, 
                                                          double X_1_2, 
                                                          double X_frac, 
                                                          double faradey_fucntions[STOKES_PARAM_NUM]) {

    /* The reference for this implementation is from Appendix B2 of https://arxiv.org/pdf/1602.03184.pdf */

    faradey_fucntions[I] = 0.0f; // This is zero by definition
    faradey_fucntions[U] = 0.0f; // This is zero by definition
    faradey_fucntions[Q] = 0.0f;
    faradey_fucntions[V] = 0.0f;

    if (X > std::numeric_limits<double>::min() && X < std::numeric_limits<double>::infinity()) {

        faradey_fucntions[Q] = 2.011 * exp(-X_frac / 4.7) - cos(X / 2) * exp(-pow(X,1.2) / 2.73) - 0.011 * exp(-X / 47.2) + 
                               (0.011 * exp(-X / 47.2) - pow(2, -1.0 / 3) / pow(3, 23. / 6) * 1e4 * M_PI * pow(X, -8.0 / 3)) / 2 * (1 + tanh(10 * log(X / 120)));
        faradey_fucntions[V] = 1.0 - 0.11 * log(1.0 + 0.035 * X);

    }

}

void Optically_Thin_Toroidal_Model::get_synchotron_emission_fit_function(Sync_emission_fit_functions e_fit_functions,
                                                                         double Emission_functions[4],
                                                                         double X,
                                                                         double X_1_2,
                                                                         double X_1_3) {

    Emission_functions[I] = 0;
    Emission_functions[U] = 0;
    Emission_functions[Q] = 0;
    Emission_functions[V] = 0;

    switch (e_fit_functions){

    case Leung_2011:

        /*

        The below expressions were originally designed to work with f_s as defined in the paper = 4 / 27 * f_crit. This is annoying for my implementation (I define X = f / f_crit),
        so I convert the expressions to work with f_crit explicitly.
        
        Ontop of that the whole function scales linearly with frequency. This is replaced by the variable "X" in the second line (in order to take the redshift into account without passing it in as an argument, because I don't like that).
        This is then corrected in the caller function, by multiplying with the critical frequency.

        Ref: https://www.aanda.org/articles/aa/pdf/2022/11/aa44339-22.pdf

        */

        Emission_functions[I] = 0.0;

        if (X_1_3 > std::numeric_limits<double>::min() && X > std::numeric_limits<double>::min() && X < std::numeric_limits<double>::infinity()) {

            Emission_functions[I] = M_SQRT2 * M_PI * Q_ELECTRON_CGS * Q_ELECTRON_CGS / 3 / C_LIGHT_CGS;
            Emission_functions[I] *= (4. / 27 * X);
            Emission_functions[I] *= ((2.5980762113) + (1.3747296369986) * 1.887749 / X_1_3); // The constants inside parentasies here are sqrt(27. / 4) and sqrt(cbrt(27. / 4)), respectfully
            Emission_functions[I] *= ((2.5980762113) + (1.3747296369986) * 1.887749 / X_1_3); // The other constant is pow (2, 11. / 12)
            Emission_functions[I] *= exp(-1.889881574 * X_1_3);                               // This constant is cbrt(27. / 4)
        }

        break;

    case Dexter_2016:

        /*
        
        Below for the I and Q functions, the frequency is replaced by the variable "X" (in order to take the redshift into account without passing it in as an argument, because I don't like that).
        This is then corrected in the caller function, by multiplying with the critical frequency.
        
        Ref: https://arxiv.org/pdf/1602.03184.pdf

        */

        if (X_1_3 > 10 * std::numeric_limits<double>::min() && X > 10 * std::numeric_limits<double>::min() && X_1_2 > 10 * std::numeric_limits<double>::min()) {

            double exponenet = exp(-1.8899 * X_1_3);

            Emission_functions[I] =     Q_ELECTRON_CGS * Q_ELECTRON_CGS / 1.73205080757 / C_LIGHT_CGS *     X * 2.5651 * (1 +  1.92 / X_1_3 + 0.9977 / X_1_3 / X_1_3) * exponenet;
            Emission_functions[Q] =     Q_ELECTRON_CGS * Q_ELECTRON_CGS / 1.73205080757 / C_LIGHT_CGS *     X * 2.5651 * (1 + 0.932 / X_1_3 + 0.4998 / X_1_3 / X_1_3) * exponenet;
            Emission_functions[V] = 4 * Q_ELECTRON_CGS * Q_ELECTRON_CGS / 1.73205080757 / C_LIGHT_CGS / 3 * X * (1.8138 / X + 3.423 / X_1_3 / X_1_3 + 0.02955 / X_1_2 + 2.0377 / X_1_3) * exponenet;
        }

        break;

    default:

        break;
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
    double T_electron     = this->get_disk_temperature(State_vector);
    double T_electron_dim = BOLTZMANN_CONST_CGS * T_electron / M_ELECTRON_CGS / C_LIGHT_CGS / C_LIGHT_CGS;
    double K0_Bessel      = std::cyl_bessel_k(0.0, 1.0 / T_electron_dim);
    double K1_Bessel      = std::cyl_bessel_k(1.0, 1.0 / T_electron_dim);
    double K2_Bessel      = std::cyl_bessel_k(2.0, 1.0 / T_electron_dim);

    /* Magnetic Field */
    double B_field_local[3]{};
    double B_CGS = this->get_magnetic_field(B_field_local, State_vector);

    /* Disk Coordinate Velocity */
    double* U_source_coord = this->get_disk_velocity(State_vector, s_Initial_conditions);

    /* Redshit */
    double redshift = Redshift(State_vector, U_source_coord);

    /* Cyclotron Frequency */
    double f_cyclo = Q_ELECTRON_CGS * B_CGS / (2 * M_PI * M_ELECTRON_CGS * C_LIGHT_CGS);

    /* The "averaged" critical frequency (without the sin(theta) term - that gets added on later from a pre-computed table) */
    double f_crit_no_sin = 3. / 2 * f_cyclo * T_electron_dim * T_electron_dim;

    /* Both the emission and faradey function expressions are in terms of an dimentionless variable X, but the definitions for X are different */
    double X_emission = 1e100;
    double X_faradey  = 1e100;

    if (f_crit_no_sin > std::numeric_limits<double>::min()) {

        X_emission = OBS_FREQUENCY_CGS / f_crit_no_sin / redshift;
        X_faradey  = T_electron_dim * sqrt(M_SQRT2 * 1e3 * f_cyclo / (OBS_FREQUENCY_CGS / redshift));
    }

    /* ==========================  Compute all the weird powers of X outside the averaging loop ========================== */

    double X_1_2_emission = sqrt(X_emission);
    double X_1_3_emission = cbrt(X_emission);

    double X_1_2_faradey  = sqrt(X_faradey);
    double X_frac_faradey =  pow(X_faradey, 1.035f);

    /* =================================================================================================================== */

    /* Average the Emission Function over all electron pitch angles */

    for (int averaging_idx = 23; averaging_idx <= 24 - 1; averaging_idx++) {

        double sin_pitch_angle = this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[averaging_idx];
        double cos_pitch_angle = this->s_Precomputed_e_pitch_angles.cos_electron_pitch_angles[averaging_idx];

        double X_1_2_emission_angle_corrected = X_1_2_emission * this->s_Precomputed_e_pitch_angles.one_over_sqrt_sin[averaging_idx];
        double X_1_3_emission_angle_corrected = X_1_3_emission * this->s_Precomputed_e_pitch_angles.one_over_cbrt_sin[averaging_idx];

        double X_faradey_angle_corrected      = X_faradey / this->s_Precomputed_e_pitch_angles.one_over_sqrt_sin[averaging_idx];;
        double X_1_2_faradey_angle_corrected  = X_1_2_faradey / this->s_Precomputed_e_pitch_angles.one_over_sin_to_1_4[averaging_idx];
        double X_frac_faradey_angle_corrected = X_frac_faradey / this->s_Precomputed_e_pitch_angles.one_over_sin_to_frac[averaging_idx];

        if (K2_Bessel > 1e10 * std::numeric_limits<double>::min()) {

            /* ================================================ The emission functions ================================================ */

            double current_emission_functions[4]{};

            this->get_synchotron_emission_fit_function(Dexter_2016, current_emission_functions, X_emission / sin_pitch_angle, X_1_2_emission_angle_corrected, X_1_3_emission_angle_corrected);

            current_emission_functions[I] *= electron_density * (f_crit_no_sin * sin_pitch_angle) / K2_Bessel; // Scale the emission by the remaining position-dependant factors
            //current_emission_functions[I] *= sin_pitch_angle * M_PI / NUM_SAMPLES_TO_AVG / 2;                  // Scale by the factors coming form averaging over all emission orientations

            current_emission_functions[Q] *= electron_density * (f_crit_no_sin * sin_pitch_angle) / K2_Bessel; // Scale the emission by the remaining position-dependant factors
            //current_emission_functions[Q] *= sin_pitch_angle * M_PI / NUM_SAMPLES_TO_AVG / 2;                  // Scale by the factors coming form averaging over all emission orientations

            if (sin_pitch_angle > std::numeric_limits<double>::min()) {

                current_emission_functions[V] *= electron_density * (f_crit_no_sin * sin_pitch_angle) / K2_Bessel; // Scale the emission by the remaining position-dependant factors
                current_emission_functions[V] *= (1. / T_electron_dim) * cos_pitch_angle / sin_pitch_angle;         // The V component has some extra angle dependance
                //current_emission_functions[V] *= sin_pitch_angle * M_PI / NUM_SAMPLES_TO_AVG / 2;                  // Scale by the factors coming form averaging over all emission orientations

            }

            // The U component is 0 by construction
            Emission_fucntions[I] += current_emission_functions[I];
            Emission_fucntions[Q] += current_emission_functions[Q];
            Emission_fucntions[V] += current_emission_functions[V];
            Emission_fucntions[U]  = 0.0f;

            /* ================================================ The faradey functions ================================================ */
            /* Originally derived in https://iopscience.iop.org/article/10.1086/592326/pdf - expressions 25, 26 and 33 */

            double current_faradey_functions[4]{};

            this->get_faradey_functions(State_vector, X_faradey_angle_corrected, X_1_2_faradey_angle_corrected, X_frac_faradey_angle_corrected, current_faradey_functions);

            double const omega_plasma_squared = 4 * M_PI * electron_density * Q_ELECTRON_CGS * Q_ELECTRON_CGS / M_ELECTRON_CGS;

            current_faradey_functions[Q] *= omega_plasma_squared * (2 * M_PI * f_cyclo) * (2 * M_PI * f_cyclo) * sin_pitch_angle * sin_pitch_angle * (K1_Bessel / K2_Bessel + 6 * T_electron_dim);
            current_faradey_functions[Q] /= 2 * C_LIGHT_CGS * (2 * M_PI * OBS_FREQUENCY_CGS / redshift) * (2 * M_PI * OBS_FREQUENCY_CGS / redshift) * (2 * M_PI * OBS_FREQUENCY_CGS / redshift);
            //current_faradey_functions[Q] *= sin_pitch_angle * M_PI / NUM_SAMPLES_TO_AVG / 2;    // Scale by the factors coming form averaging over all emission orientations

            current_faradey_functions[V] *= omega_plasma_squared * (2 * M_PI * f_cyclo) * cos_pitch_angle * (K0_Bessel) / K2_Bessel;
            current_faradey_functions[V] /= C_LIGHT_CGS * (2 * M_PI * OBS_FREQUENCY_CGS / redshift) * (2 * M_PI * OBS_FREQUENCY_CGS / redshift);
            //current_faradey_functions[V] *= sin_pitch_angle * M_PI / NUM_SAMPLES_TO_AVG / 2;    // Scale by the factors coming form averaging over all emission orientations

            // The I and U components are 0 by definition
            Faradey_functions[Q] = current_faradey_functions[Q];
            Faradey_functions[V] = current_faradey_functions[V];
            Faradey_functions[I] = 0.0f;
            Faradey_functions[U] = 0.0f;

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
    double& emission_scale     = this->s_Emission_Parameters.Emission_scale;

    double electron_density = emission_scale * get_disk_density_profile(State_vector);

    /* Disk Coordinate Velocity */

    double* U_source_coord = this->get_disk_velocity(State_vector, s_Initial_conditions);

    double redshift = Redshift(State_vector, U_source_coord);

    Emission_functions[I] = electron_density * pow(redshift, emission_power_law);

    for (int stokes_index = 1; stokes_index <= STOKES_PARAM_NUM - 1; stokes_index++) {

        Emission_functions[stokes_index] = 0;

    }

}

double Optically_Thin_Toroidal_Model::get_absorbtion_function_phenomenological(double Emission_Functions, double State_vector[], double redshift) {

    double& abs_coeff          = this->s_Emission_Parameters.Absorbtion_coeff;
    double& emission_scale     = this->s_Emission_Parameters.Emission_scale;
    double& source_f_power_law = this->s_Emission_Parameters.Source_f_power_law;
    double& emission_power_law = this->s_Emission_Parameters.Emission_power_law;

    return abs_coeff * emission_scale * this->get_disk_density_profile(State_vector) * pow(redshift, source_f_power_law + emission_power_law);

}

void Optically_Thin_Toroidal_Model::precompute_electron_pitch_angles() {

    for (int index = 0; index <= NUM_SAMPLES_TO_AVG - 1; index++) {

        double pitch_angle = double(index) / NUM_SAMPLES_TO_AVG * M_PI;
        this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index] = sin(pitch_angle);
        this->s_Precomputed_e_pitch_angles.cos_electron_pitch_angles[index] = cos(pitch_angle);

        if (this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index] != 0) {

            this->s_Precomputed_e_pitch_angles.one_over_sqrt_sin[index]    = 1. / sqrt(this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index]);
            this->s_Precomputed_e_pitch_angles.one_over_sin_to_1_4[index]  = 1. / sqrt(this->s_Precomputed_e_pitch_angles.one_over_sqrt_sin[index]);
            this->s_Precomputed_e_pitch_angles.one_over_cbrt_sin[index]    = 1. / cbrt(this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index]);
            this->s_Precomputed_e_pitch_angles.one_over_sin_to_frac[index] = 1. / pow(this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index], 1.035f);

        }

    }

}

double get_planck_function_CGS(double Frequency, double Temperature) {

    return 2 * PLANCK_CONSTANT_CGS * Frequency * Frequency * Frequency / C_LIGHT_CGS / C_LIGHT_CGS / (exp(PLANCK_CONSTANT_CGS * Frequency / BOLTZMANN_CONST_CGS / Temperature) - 1);

}

