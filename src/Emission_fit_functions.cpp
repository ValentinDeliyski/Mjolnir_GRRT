#include "General_math_functions.h"
#include "General_GR_functions.h"
#include "Disk_Models.h"
#include "Spacetimes.h"
#include "Constants.h"

#include "gsl/gsl_sf_hyperg.h"

/* =============================================== Thermal Synchrotron Transfer Functions =============================================== */

void Generic_Optically_Thin_Model::get_thermal_synchrotron_emission_fit_functions(const Thermal_transfer_f_arguments_type* const p_Transfer_arags,
                                                                                  double* const Emission_functions) const {

    /* The reference for the expressions below can be found in https://iopscience.iop.org/article/10.3847/1538-4357/ac1b28/pdf - equations (30). */

    /* Zero out the emission functions just in case. */
    memset(Emission_functions, 0, STOKES_PARAM_NUM * sizeof(double));

    /* Check weather the fit function p_Transfer_args are numerically OK to use in the expressions - they have problems at velry low densities where the emission/absorbtion is negligable.
       In such cases I directly return. */
    if (isnan(p_Transfer_arags->X) || isinf(p_Transfer_arags->X) || isinf(1.0 / p_Transfer_arags->cbrt_X)) { return; }

    /* This silly magic number pops up as a coefficient in the fit functions. */
    constexpr double TWO_TO_11_OVER_12 = 1.887749;

    /* The common exponential factor for each polarization component. */
    double exponent = exp(-p_Transfer_arags->cbrt_X);

    Emission_functions[I] = M_SQRT2 * M_PI / 27.0 * p_Transfer_arags->sin_pitch_angle * p_Transfer_arags->X * (1 + TWO_TO_11_OVER_12 / p_Transfer_arags->cbrt_X) * (1 + TWO_TO_11_OVER_12 / p_Transfer_arags->cbrt_X) * exponent;

    /* Return if the simulataion does not include polarization components */
    if (!this->Include_polarization) { return; }

    Emission_functions[Q] = -M_SQRT2 * M_PI / 27.0 * p_Transfer_arags->sin_pitch_angle * p_Transfer_arags->X
                            * (1 + (7. * p_Transfer_arags->T_electron_dim_to_24_25 + 35.) / (10. * p_Transfer_arags->T_electron_dim_to_24_25 + 75.) * TWO_TO_11_OVER_12 / p_Transfer_arags->cbrt_X)
                            * (1 + (7. * p_Transfer_arags->T_electron_dim_to_24_25 + 35.) / (10. * p_Transfer_arags->T_electron_dim_to_24_25 + 75.) * TWO_TO_11_OVER_12 / p_Transfer_arags->cbrt_X)
                            * exponent;

    Emission_functions[V] = p_Transfer_arags->cos_pitch_angle / p_Transfer_arags->T_electron_dim
                            * (M_PI / 3 + M_PI / 3 * p_Transfer_arags->cbrt_X + (2. / 300) * p_Transfer_arags->sqrt_X + (2 * M_PI / 19.) * p_Transfer_arags->sqrt_X * p_Transfer_arags->sqrt_X) * exponent;

}

void Generic_Optically_Thin_Model::get_thermal_synchrotron_absorbtion_fit_functions(const Thermal_transfer_f_arguments_type* const Transfer_arags,
                                                                                    const Emission_medium_state_type* p_Emission_medium_state,
                                                                                    const double* const Emission_function,
                                                                                    double* const Absorbtion_function) const {

    /* The reference for the expressions below can be found in https://iopscience.iop.org/article/10.3847/1538-4357/ac1b28/pdf - equation (32). */

    /* Zero out the absorbtion functions just in case. */
    memset(Absorbtion_function, 0, STOKES_PARAM_NUM * sizeof(double));

    const double& frequency   = Transfer_arags->frequency;
    const double exp_argument = PLANCK_CONSTANT_SI * frequency / (BOLTZMANN_CONST_SI * p_Emission_medium_state->Temperature);
    const double f_cyclo      = Q_ELECTRON_CGS * p_Emission_medium_state->Magnetic_fields.B_field_plasma_frame_norm / (2 * M_PI * M_ELECTRON_CGS * C_LIGHT_CGS);

    /* Check weather the exponent argument is numerically OK to use in the fit functions. */
    if (isnan(exp_argument) || isinf(exp_argument)) { return; }
    
    for (int stokes_idx = 0; stokes_idx <= STOKES_PARAM_NUM - 1; stokes_idx++) {

        Absorbtion_function[stokes_idx] = Emission_function[stokes_idx] * M_ELECTRON_CGS * C_LIGHT_CGS * C_LIGHT_CGS / 2 / PLANCK_CONSTANT_CGS / frequency / frequency * f_cyclo * (exp(exp_argument) - 1);

    }

}

void Generic_Optically_Thin_Model::get_thermal_synchrotron_faradey_fit_functions(const Thermal_transfer_f_arguments_type* const Transfer_arags,
                                                                                 double* const Faradey_fucntions) const {

    /* The reference for this implementation is from Appendix B2 of https://iopscience.iop.org/article/10.3847/1538-4357/ac1b28/pdf, expressions (33) to (37). */

    /* Zero out the Faradey functions just in case. */
    memset(Faradey_fucntions, 0, STOKES_PARAM_NUM * sizeof(double));

    /* Return if the simulataion does not include polarization components */
    if (!this->Include_polarization) { return; }

    double const K0_Bessel = std::cyl_bessel_k(0.0, 1.0 / Transfer_arags->T_electron_dim);
    double const K1_Bessel = std::cyl_bessel_k(1.0, 1.0 / Transfer_arags->T_electron_dim);
    double const K2_Bessel = std::cyl_bessel_k(2.0, 1.0 / Transfer_arags->T_electron_dim);

    if (isnan(Transfer_arags->X) || isinf(Transfer_arags->X) || isinf(1.0 / Transfer_arags->X) || isinf(1e10 / K2_Bessel)) { return; }

    double const common_exp_term = exp(-1.699 / Transfer_arags->sqrt_X);

    double const f_0 = 2.001 * exp(-19.78 / Transfer_arags->X_to_0_p_5175)
                     - cos(39.89 / Transfer_arags->sqrt_X) * exp(-70.16 / Transfer_arags->X_to_0_p_6)
                     - 0.011 * common_exp_term;

    double const f_m = f_0
                     + (0.011 * common_exp_term - 0.003135 * Transfer_arags->cbrt_X * Transfer_arags->cbrt_X * Transfer_arags->cbrt_X * Transfer_arags->cbrt_X)
                     * 0.5 * (1 + tanh(10 * log(0.6648 / Transfer_arags->sqrt_X)));

    double const delta_J_5 = 0.4379 * log(1 + 1.3414 / Transfer_arags->X_to_0_p_7515);

    Faradey_fucntions[Q] = -f_m * Transfer_arags->sin_pitch_angle * Transfer_arags->sin_pitch_angle * (K1_Bessel / K2_Bessel + 6 * Transfer_arags->T_electron_dim);

    Faradey_fucntions[V] = (K0_Bessel - delta_J_5) / K2_Bessel * Transfer_arags->cos_pitch_angle;

}

/* ========================================== Kappa Synchrotron Transfer Functions ========================================== */

void Generic_Optically_Thin_Model::get_kappa_synchrotron_emission_fit_functions(const Kappa_transfer_f_arguments_type* const p_Transfer_args,
                                                                                double* const Emission_functions) const {

    // The reference for these expressions is https://arxiv.org/pdf/1602.08749, equations (35), (36), (37) and (38).

    /* Zero out the emission functions just in case. */
    memset(Emission_functions, 0, STOKES_PARAM_NUM * sizeof(double));

    /* Check weather the fit function p_Transfer_args are numerically OK to use in the expressions - they have problems at velry low densities where the emission/absorbtion is negligable. 
       In such cases I directly return. */
    if (isnan(p_Transfer_args->sqrt_X) || isinf(p_Transfer_args->sqrt_X) || isnan(p_Transfer_args->X) || isinf(p_Transfer_args->X) || isnan(1. / p_Transfer_args->X_to_7_over_20) || isinf(1. / p_Transfer_args->X_to_7_over_20) || isinf(1. / p_Transfer_args->T_electron_dim) || isinf(1. / p_Transfer_args->sin_emission_angle)) {

        return;

    }

    // ------------------------------------------------------------------------ Low frequency fit coefficient ------------------------------------------------------------------------ //

    constexpr double THREE_TO_7_OVER_3 = 12.980246132766677;

    double Emission_functions_low[STOKES_PARAM_NUM]{};
    double Common_factor_low = p_Transfer_args->cbrt_X * p_Transfer_args->sin_emission_angle * (4 * M_PI / THREE_TO_7_OVER_3) * std::tgamma(p_Transfer_args->kappa - 4.0 / 3) / std::tgamma(p_Transfer_args->kappa - 2.0);

    // ----------------------------------------------------------------------- High frequency fit coefficient ------------------------------------------------------------------------ //

    double Emission_functions_high[STOKES_PARAM_NUM]{};
    double Common_factor_high = pow(p_Transfer_args->X, -(p_Transfer_args->kappa - 2.) / 2.) * p_Transfer_args->sin_emission_angle * pow(3.0, (p_Transfer_args->kappa - 1.) / 2) * (p_Transfer_args->kappa - 2.) * (p_Transfer_args->kappa - 1.) / 4 
                              * std::tgamma(p_Transfer_args->kappa / 4 - 1.0 / 3) * std::tgamma(p_Transfer_args->kappa / 4 + 4.0 / 3);

    // ------------------------------------------------------------------------ Bridging function ------------------------------------------------------------------------ //

    double power_I = 3 * pow(p_Transfer_args->kappa, -3. / 2);

    if (!isnan(Emission_functions_high[I] / Emission_functions_low[I])) {

        Emission_functions_low[I]  = Common_factor_low;
        Emission_functions_high[I] = Common_factor_high;

        Emission_functions[I] = Emission_functions_low[I] * pow(1. + pow(Emission_functions_high[I] / Emission_functions_low[I], -power_I), -1. / power_I);

        //Emission_functions[I] = Emission_functions_high[I];

    }

    /* Return if the simulataion does not include polarization components */
    if (!this->Include_polarization) { return; }

    double power_Q = 3.7 * pow(p_Transfer_args->kappa, -8. / 5);

    if (!isnan(Emission_functions_high[Q] / Emission_functions_low[Q])) {

        Emission_functions_low[Q]  = -Common_factor_low / 2;
        Emission_functions_high[Q] = -Common_factor_high * (16.0 / 25 + p_Transfer_args->kappa / 50);

        Emission_functions[Q] = Emission_functions_low[Q] * pow(1. + pow(Emission_functions_high[Q] / Emission_functions_low[Q], -power_Q), -1. / power_Q);
    }

    double power_V = 13. / 5 * pow(p_Transfer_args->kappa, -36. / 25);

    if (!isnan(Emission_functions_high[V] / Emission_functions_low[V])) {

        Emission_functions_low[V]  = -Common_factor_low * (9.0 / 16 * pow(pow(p_Transfer_args->sin_emission_angle, -12.0 / 5) - 1., 12.0 / 25)) * pow(p_Transfer_args->kappa, -66.0 / 125) / p_Transfer_args->T_electron_dim / p_Transfer_args->X_to_7_over_20;
        Emission_functions_high[V] = -Common_factor_high * (49.0 / 64 * pow(pow(p_Transfer_args->sin_emission_angle, -5.0 / 2) - 1, 11.0 / 25)) * pow(p_Transfer_args->kappa, -11.0 / 25) / p_Transfer_args->T_electron_dim / p_Transfer_args->sqrt_X;

        Emission_functions[V] = Emission_functions_low[V] * pow(1. + pow(Emission_functions_high[V] / Emission_functions_low[V], -power_V), -1. / power_V);

    }

}

void Generic_Optically_Thin_Model::get_kappa_synchrotron_absorbtion_fit_functions(const Kappa_transfer_f_arguments_type* const p_Transfer_args,
                                                                                  double* const Absorbtion_functions) const {

    // The reference for these expressions is https://arxiv.org/pdf/1602.08749, equations (39), (40), (41) and (42).
    
    /* Zero out the emission functions just in case. */
    memset(Absorbtion_functions, 0, STOKES_PARAM_NUM * sizeof(double));

    /* Check weather the fit function p_Transfer_args are numerically OK to use in the expressions - they have problems at velry low densities where the emission/absorbtion is negligable.
       In such cases I directly return. */
    if (isnan(p_Transfer_args->sqrt_X) || isinf(p_Transfer_args->sqrt_X) || isnan(p_Transfer_args->X) || isinf(p_Transfer_args->X) || isnan(1. / p_Transfer_args->X_to_7_over_20) || isinf(1. / p_Transfer_args->X_to_7_over_20) || isinf(1. / p_Transfer_args->T_electron_dim) || isinf(1. / p_Transfer_args->sin_emission_angle)) {

        return;

    }

    // ------------------------------------------------------------------------ Low frequency fit coefficients ------------------------------------------------------------------------ //

    constexpr double THREE_TO_1_OVER_6 = 1.2009369551760027;
    constexpr double GAMMA_OF_5_OVER_3 = 0.90274529295;

    // Below are the coefficients, present in the 2F1 hypergoemetric function from expression (39). Bcause in general |-kappa * T_electron_dim| > 1, I will use the algebraic relation 
    // from Abramowitz and Stegun 15.3.8. Note that our equivalent to the argument z in 15.3.8 is strictly real and negative, so taking fractional powers of it returns a real number.

    double a = p_Transfer_args->kappa - 1.0 / 3;
    double b = p_Transfer_args->kappa + 1.0;
    double c = p_Transfer_args->kappa + 2.0 / 3;
    double z = -p_Transfer_args->kappa * p_Transfer_args->T_electron_dim;

    double _2F1{};

    if (-z < 1.0) {

        _2F1 = gsl_sf_hyperg_2F1(a, b, c, z);

    }
    else {

        _2F1 = pow(1. - z, -a) * std::tgamma(c) / std::tgamma(b) * std::tgamma(b - a) / std::tgamma(c - a) * gsl_sf_hyperg_2F1(a, c - b, a - b + 1., 1. / (1. - z))
             + pow(1. - z, -b) * std::tgamma(c) / std::tgamma(a) * std::tgamma(a - b) / std::tgamma(c - b) * gsl_sf_hyperg_2F1(b, c - a, b - a + 1., 1. / (1 - z));

    }

    double Absorbtion_functions_low[STOKES_PARAM_NUM]{};
    double Common_factor_low = 1.0 / p_Transfer_args->cbrt_X / p_Transfer_args->cbrt_X * THREE_TO_1_OVER_6 * 10.0 / 41 * 2 * M_PI / pow(p_Transfer_args->T_electron_dim * p_Transfer_args->kappa, 10.0 / 3 - p_Transfer_args->kappa) * (p_Transfer_args->kappa - 2) * (p_Transfer_args->kappa - 1) * p_Transfer_args->kappa / (3 * p_Transfer_args->kappa - 1)
                             * GAMMA_OF_5_OVER_3 * _2F1;

    // ----------------------------------------------------------------------- High frequency fit coefficients ------------------------------------------------------------------------ //

    double Absorbtion_functions_high[STOKES_PARAM_NUM]{};
    double Common_factor_high = pow(p_Transfer_args->X, -(1 + p_Transfer_args->kappa) / 2) * M_PI * (2 / M_2_SQRTPI) / 3 * (p_Transfer_args->kappa - 2) * (p_Transfer_args->kappa - 1) * p_Transfer_args->kappa / (p_Transfer_args->kappa * p_Transfer_args->T_electron_dim) / (p_Transfer_args->kappa * p_Transfer_args->T_electron_dim) / (p_Transfer_args->kappa * p_Transfer_args->T_electron_dim)
                              * (2 * std::tgamma(2 + p_Transfer_args->kappa / 2) / (2 + p_Transfer_args->kappa) - 1.0);


    // ------------------------------------------------------------------------ Bridging function ------------------------------------------------------------------------ //

    double power_I = pow(-7.0 / 4 + 8.0 / 5 * p_Transfer_args->kappa, -43.0 / 50);

    if (!isnan(Absorbtion_functions_high[I] / Absorbtion_functions_low[I])) {

        Absorbtion_functions_low[I]  = Common_factor_low;
        Absorbtion_functions_high[I] = Common_factor_high * (pow(3.0 / p_Transfer_args->kappa, 19.0 / 4) + 3.0 / 5);

        Absorbtion_functions[I] = Absorbtion_functions_low[I] * pow(1. + pow(Absorbtion_functions_high[I] / Absorbtion_functions_low[I], -power_I), -1.0 / power_I);

        //Absorbtion_functions[I] = Absorbtion_functions_high[I];
    }

    /* Return if the simulataion does not include polarization components */
    if (!this->Include_polarization) { return; }

    double power_Q = 7.0 / 5 * pow(p_Transfer_args->kappa, -23.0 / 20);

    if (!isnan(Absorbtion_functions_high[Q] / Absorbtion_functions_low[Q])) {


        Absorbtion_functions_low[Q]  = -Common_factor_low * 25.0 / 48;
        Absorbtion_functions_high[Q] = -Common_factor_high * (441 * pow(p_Transfer_args->kappa, -144.0 / 25) + 11.0 / 20);

        Absorbtion_functions[Q] = Absorbtion_functions_low[Q] * pow(1. + pow(Absorbtion_functions_high[Q] / Absorbtion_functions_low[Q], -power_Q), -1.0 / power_Q);

    }

    double power_V = 61.0 / 50 * pow(p_Transfer_args->kappa, -142.0 / 125) + 7.0 / 1000;

    if (!isnan(Absorbtion_functions_high[V] / Absorbtion_functions_low[V])) {

        Absorbtion_functions_low[V]  = -Common_factor_low * pow((pow(p_Transfer_args->sin_emission_angle, -114.0 / 50) - 1), 223.0 / 500) / p_Transfer_args->X_to_7_over_20 * pow(p_Transfer_args->kappa, -7.0 / 10);
        Absorbtion_functions_high[V] = -Common_factor_high * 143.0 / 10 * pow(p_Transfer_args->T_electron_dim, -116.0 / 125) * sqrt(pow(p_Transfer_args->sin_emission_angle, -41.0 / 20) - 1) * (169 * pow(p_Transfer_args->kappa, -8) + 13.0 / 2500 * p_Transfer_args->kappa - 1.0 / 200 + 47.0 / 200 / p_Transfer_args->kappa) / p_Transfer_args->sqrt_X;

        Absorbtion_functions[V] = Absorbtion_functions_low[V] * pow(1. + pow(Absorbtion_functions_high[V] / Absorbtion_functions_low[V], -power_V), -1.0 / power_V);

    }

}

/* ========================================== Phenomenological Synchrotron Transfer Functions ========================================== */

void Generic_Optically_Thin_Model::get_phenomenological_synchrotron_fit_functions(const Phenomenological_transfer_f_arguments_type* const p_Transfer_args,
                                                                                  Transfer_functions_type* const p_Transfer_functions) const {

    /* The reference for this implementation is https://iopscience.iop.org/article/10.3847/1538-4357/ab96c6 - expressions (9) and (11). */

    /* === Zero out the transfer functions just in case === */
    memset(p_Transfer_functions, 0, sizeof(Transfer_functions_type));

    /* This is common factor infront of all emission functions (up to a density - that gets added in the caller function), which carrires most of the units.
       For the sake of unifying all emission models I want to add on this factor in a single caller function, for all emission models (Thermal, kappa and so on).
       This phenomenological model does not have this factor included for the sake of simplicity (unlike the other models). For this reason I divide by it
       in this function, while the caller multiplies by it afterwards. This amounts to multiplying by one, but it allows me to have only one function (the caller)
       that adds on the units to the transfer functions. This makes reducing the radiative transfer to a dimensionless problem a lot more straightforward. */
    double common_factor_emission = Q_ELECTRON_CGS * Q_ELECTRON_CGS / C_LIGHT_CGS * p_Transfer_args->f_cyclo;

    /* The common factor for the absorbtion function is different, but it serves the same purpose. */
    double common_factor_absorbtion = Q_ELECTRON_CGS * Q_ELECTRON_CGS / M_ELECTRON_CGS / C_LIGHT_CGS / p_Transfer_args->frequency;

    /* References for the sake of readability. */
    const double& emission_power_law = this->s_Emission_params.Phenomenological_emission_power_law;
    const double& source_f_power_law = this->s_Emission_params.Phenomenological_source_f_power_law;
    const double& emission_coeff     = this->s_Emission_params.Phenomenological_emission_coeff;
    const double& abs_coeff          = this->s_Emission_params.Phenomenological_absorbtion_coeff;

    p_Transfer_functions->Emission_functions[I] = (emission_coeff / this->s_Disk_params.Electron_density_scale / common_factor_emission) * pow(p_Transfer_args->redshift, emission_power_law);

    p_Transfer_functions->Absorbtion_functions[I] = (abs_coeff * emission_coeff / this->s_Disk_params.Electron_density_scale / common_factor_absorbtion) * pow(p_Transfer_args->redshift, source_f_power_law + emission_power_law);

}
