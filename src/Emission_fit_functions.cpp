#include "General_math_functions.h"
#include "General_GR_functions.h"
#include "Disk_Models.h"
#include "Spacetimes.h"
#include "Constants.h"

#include "gsl/gsl_sf_hyperg.h"

/* =============================================== Thermal synchrotron Transfer Functions =============================================== */


void Generic_Optically_Thin_Model::get_thermal_synchrotron_fit_functions(double Emission_fucntions[STOKES_PARAM_NUM],
                                                                        double Faradey_functions[STOKES_PARAM_NUM],
                                                                        Thermal_emission_f_arguments* Emission_args,
                                                                        Thermal_faradey_f_arguments* Faradey_args) {

    if (this->Include_polarization) {

        this->get_thermal_synchrotron_emission_fit_functions(Dexter_2016, Emission_fucntions, Emission_args->X, Emission_args->sqrt_X, Emission_args->cbrt_X);
        this->get_thermal_synchrotron_faradey_fit_functions(Faradey_args->X, Faradey_args->X_to_1_point_2, Faradey_args->X_to_1_point_035, Faradey_functions);

    }
    else {

        this->get_thermal_synchrotron_emission_fit_functions(Leung_2011, Emission_fucntions, Emission_args->X, Emission_args->sqrt_X, Emission_args->cbrt_X);

    }
}


//! Evaluates the thermal sychrotron fit functions.
/*! Evaluates the thermal sychrotron fit functions (emission and Fatadey), based on two sources:
*       1) https://www.aanda.org/articles/aa/pdf/2022/11/aa44339-22.pdf for the unpolarized case.
*       2) https://arxiv.org/pdf/1602.03184.pdf for the polarized case
*   
*   \param [in] e_fit_functions - Enum that selects which fit functions to evaluate.
*   \param [out] Emission_functions - Pointer to the array that holds the fit funct ions.
*   \param [in] X - Dimentionless parameter that the fit functions depend on.
*   \param [in] sqrt_X - Square root of X - computed ouside this function because it saves clock cycles when averaging over pitch angles.
*   \param [in] cbrt_X - Cube root of X - computed ouside this function because it saves clock cycles when averaging over pitch angles.
*/
void Generic_Optically_Thin_Model::get_thermal_synchrotron_emission_fit_functions(Thermal_Syncotron_fit_selector e_fit_functions,
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

        /*!< The below expressions were originally designed to work with f_s as defined in the paper = 4 / 27 * f_crit. This is annoying for my implementation (I define X = f / f_crit),
         *   so I convert the expressions to work with f_crit explicitly.
         *
         *   Ontop of that the whole function scales linearly with frequency. This is replaced by the variable "X" in the second line (in order to take the redshift into account without passing it in as an argument, because I don't like that).
         *   This is then corrected in the caller function, by multiplying with the critical frequency.
         *   
         *   Ref: https://www.aanda.org/articles/aa/pdf/2022/11/aa44339-22.pdf
         */
        if (!isnan(cbrt_X) && !isinf(cbrt_X) && !isnan(X) && !isinf(X) && !isnan(1.0 / cbrt_X) && !isinf(1.0 / cbrt_X)) {

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

        if (!isnan(cbrt_X) && !isinf(cbrt_X) && !isnan(X) && !isinf(X) && !isnan(sqrt_X) && !isinf(sqrt_X) && !isnan(1.0 / cbrt_X) && !isinf(1.0 / cbrt_X)) {

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

void Generic_Optically_Thin_Model::get_thermal_synchrotron_faradey_fit_functions(double X,
                                                                                double X_to_1_point_2,
                                                                                double X_frac,
                                                                                double faradey_fucntions[STOKES_PARAM_NUM]) {

    /* The reference for this implementation is from Appendix B2 of https://arxiv.org/pdf/1602.03184.pdf */

    faradey_fucntions[I] = 0.0f; // This is zero by definition
    faradey_fucntions[U] = 0.0f; // This is zero by definition
    faradey_fucntions[Q] = 0.0f;
    faradey_fucntions[V] = 0.0f;

    constexpr double TWO_TO_MINUS_ONE_THIRD = 0.7937005259840997373758528196361;
    constexpr double THREE_TO_23_OVER_6 = 67.44733739;

    double exp_term = exp(-X / 47.2);

    if (!isnan(X) && !isinf(X) && !isinf(exp_term)) {

        faradey_fucntions[Q] = 2.011 * exp(-X_frac / 4.7) - cos(X / 2) * exp(-X_to_1_point_2 / 2.73) - 0.011 * exp_term
            + (0.011 * exp_term - TWO_TO_MINUS_ONE_THIRD / THREE_TO_23_OVER_6 * 1e4 * M_PI * pow(X, -8.0 / 3)) / 2 * (1 + tanh(10 * log(X / 120)));
        faradey_fucntions[V] = 1.0 - 0.11 * log(1.0 + 0.035 * X);

    }

}

/* ========================================== Kappa synchrotron Transfer Functions ========================================== */

void Generic_Optically_Thin_Model::get_kappa_synchrotron_fit_functions(double Emission_fucntions[STOKES_PARAM_NUM],
                                                                      double Faradey_functions[STOKES_PARAM_NUM], 
                                                                      double Absorbtion_fucntions[STOKES_PARAM_NUM],
                                                                      Kappa_transfer_f_arguments* args){

    this->get_kappa_synchrotron_emission_fit_functions(Emission_fucntions, args->X, args->sqrt_X, args->cbrt_X, args->X_to_7_over_20, args->kappa, args->sin_emission_angle, args->T_electron_dim);
    this->get_kappa_synchrotron_absorbtion_fit_functions(Absorbtion_fucntions, args->X, args->sqrt_X, args->cbrt_X, args->X_to_7_over_20, args->kappa, args->sin_emission_angle, args->T_electron_dim);

    for (int index = 0; index <= STOKES_PARAM_NUM - 1; index++) {

        Faradey_functions[index] = 0.0;

    }

}

void Generic_Optically_Thin_Model::get_kappa_synchrotron_emission_fit_functions(double Emission_functions[STOKES_PARAM_NUM],
                                                                               double X,
                                                                               double sqrt_X,
                                                                               double cbrt_X,
                                                                               double X_to_7_over_20,
                                                                               double kappa,
                                                                               double sin_emission_angle,
                                                                               double T_electron_dim) {

    // The reference for these expressions is https://arxiv.org/pdf/1602.08749, equations (35), (36), (37) and (38).

    Emission_functions[I] = 0.0;
    Emission_functions[U] = 0.0;
    Emission_functions[Q] = 0.0;
    Emission_functions[V] = 0.0;

    if (!isnan(sqrt_X) && !isinf(sqrt_X) && !isnan(X) && !isinf(X) && !isnan(1. / X_to_7_over_20) && !isinf(1. / X_to_7_over_20) && !isinf(1. / T_electron_dim) && !isinf(1. / sin_emission_angle)) {

        // ------------------------------------------------------------------------ Low frequency fit ------------------------------------------------------------------------ //

        constexpr double THREE_TO_7_OVER_3 = 12.980246132766677;

        double Emission_functions_low[STOKES_PARAM_NUM]{};
        double Common_factor_low = cbrt_X * sin_emission_angle * (4 * M_PI / THREE_TO_7_OVER_3) * std::tgamma(kappa - 4.0 / 3) / std::tgamma(kappa - 2.0);

        Emission_functions_low[I] = Common_factor_low;
        Emission_functions_low[Q] = -Common_factor_low / 2;
        Emission_functions_low[U] = 0;
        Emission_functions_low[V] = -Common_factor_low * (9.0 / 16 * pow(pow(sin_emission_angle, -12.0 / 5) - 1., 12.0 / 25)) * pow(kappa, -66.0 / 125) / T_electron_dim / X_to_7_over_20;

        // ----------------------------------------------------------------------- High frequency fit ------------------------------------------------------------------------ //

        double Emission_functions_high[STOKES_PARAM_NUM]{};
        double Common_factor_high = pow(X, -(kappa - 2.) / 2.) * sin_emission_angle * pow(3, (kappa - 1) / 2) * (kappa - 2) * (kappa - 1) / 4 * std::tgamma(kappa / 4 - 1.0 / 3) * std::tgamma(kappa / 4 + 4.0 / 3);

        Emission_functions_high[I] = Common_factor_high;
        Emission_functions_high[Q] = -Common_factor_high * (16.0 / 25 + kappa / 50);
        Emission_functions_high[U] = 0;
        Emission_functions_high[V] = -Common_factor_high * (49.0 / 64 * pow(pow(sin_emission_angle, -5.0 / 2) - 1, 11.0 / 25)) * pow(kappa, -11.0 / 25) / T_electron_dim / sqrt_X;

        // ------------------------------------------------------------------------ Bridging function ------------------------------------------------------------------------ //

        double power_I = 3 * pow(kappa, -3. / 2);

        if (!isnan(Emission_functions_high[I] / Emission_functions_low[I])) {

            Emission_functions[I] = Emission_functions_low[I] * pow(1. + pow(Emission_functions_high[I] / Emission_functions_low[I], -power_I), -1. / power_I);

        }

        double power_Q = 3.7 * pow(kappa, -8. / 5);

        if (!isnan(Emission_functions_high[Q] / Emission_functions_low[Q])) {

            Emission_functions[Q] = Emission_functions_low[Q] * pow(1. + pow(Emission_functions_high[Q] / Emission_functions_low[Q], -power_Q), -1. / power_Q);
        }

        Emission_functions[U] = 0;

        double power_V = 13. / 5 * pow(kappa, -36. / 25);

        if (!isnan(Emission_functions_high[V] / Emission_functions_low[V])) {

            Emission_functions[V] = Emission_functions_low[V] * pow(1. + pow(Emission_functions_high[V] / Emission_functions_low[V], -power_V), -1. / power_V);

        }

       //Emission_functions[I] = Emission_functions_high[I];

    }

}

void Generic_Optically_Thin_Model::get_kappa_synchrotron_absorbtion_fit_functions(double Absorbtion_functions[STOKES_PARAM_NUM],
                                                                                 double X,
                                                                                 double sqrt_X,
                                                                                 double cbrt_X,
                                                                                 double X_to_7_over_20,
                                                                                 double kappa,
                                                                                 double sin_emission_angle,
                                                                                 double T_electron_dim) {

    // The reference for these expressions is https://arxiv.org/pdf/1602.08749, equations (39), (40), (41) and (42).

    Absorbtion_functions[I] = 0.0;
    Absorbtion_functions[U] = 0.0;
    Absorbtion_functions[Q] = 0.0;
    Absorbtion_functions[V] = 0.0;
    double Absorbtion_functions_high[STOKES_PARAM_NUM]{};
    if (!isnan(sqrt_X) && !isinf(sqrt_X) && !isnan(X) && !isinf(X) && !isnan(1.0 / X_to_7_over_20) && !isinf(1.0 / X_to_7_over_20)) {

        // ------------------------------------------------------------------------ Low frequency fit ------------------------------------------------------------------------ //

        constexpr double THREE_TO_1_OVER_6 = 1.2009369551760027;
        constexpr double GAMMA_OF_5_OVER_3 = 0.90274529295;

        // Below are the coefficients, present in the 2F1 hypergoemetric function from expression (39). Bcause in general |(-kappa * T_electron_dim)| > 1, I will use the algebraic relation 
        // from Abramowitz and Stegun 15.3.8. Note that our equivalent to the argument z in 15.3.8 is strictly real and negative, so taking fractional powers of it returns a real number.

        double a = kappa - 1.0 / 3;
        double b = kappa + 1.0;
        double c = kappa + 2.0 / 3;
        double z = -kappa * T_electron_dim;

        double _2F1{};

        if (-z < 1) {

            _2F1 = gsl_sf_hyperg_2F1(a, b, c, z);

        }
        else {

            _2F1 = pow(1. - z, -a) * std::tgamma(c) / std::tgamma(b) * std::tgamma(b - a) / std::tgamma(c - a) * gsl_sf_hyperg_2F1(a, c - b, a - b + 1., 1. / (1. - z))
                 + pow(1. - z, -b) * std::tgamma(c) / std::tgamma(a) * std::tgamma(a - b) / std::tgamma(c - b) * gsl_sf_hyperg_2F1(b, c - a, b - a + 1., 1. / (1 - z));

        }

        double Absorbtion_functions_low[STOKES_PARAM_NUM]{};
        double Common_factor_low = 1.0 / cbrt_X / cbrt_X * THREE_TO_1_OVER_6 * 10.0 / 41 * 2 * M_PI / pow(T_electron_dim * kappa, 10.0 / 3 - kappa) * (kappa - 2) * (kappa - 1) * kappa / (3 * kappa - 1)
            * GAMMA_OF_5_OVER_3 * _2F1;

        Absorbtion_functions_low[I] = Common_factor_low;
        Absorbtion_functions_low[Q] = -Common_factor_low * 25.0 / 48;
        Absorbtion_functions_low[U] = 0;
        Absorbtion_functions_low[V] = -Common_factor_low * pow((pow(sin_emission_angle, -114.0 / 50) - 1), 223.0 / 500) / X_to_7_over_20 * pow(kappa, -7.0 / 10);

        // ----------------------------------------------------------------------- High frequency fit ------------------------------------------------------------------------ //

        double Common_factor_high = pow(X, -(1 + kappa) / 2) * M_PI * (2 / M_2_SQRTPI) / 3 * (kappa - 2) * (kappa - 1) * kappa / (kappa * T_electron_dim) / (kappa * T_electron_dim) / (kappa * T_electron_dim)
            * (2 * std::tgamma(2 + kappa / 2) / (2 + kappa) - 1.0);

        Absorbtion_functions_high[I] = Common_factor_high * (pow(3.0 / kappa, 19.0 / 4) + 3.0 / 5);
        Absorbtion_functions_high[Q] = -Common_factor_high * (441 * pow(kappa, -144.0 / 25) + 11.0 / 20);
        Absorbtion_functions_high[U] = 0;
        Absorbtion_functions_high[V] = -Common_factor_high * 143.0 / 10 * pow(T_electron_dim, -116.0 / 125) * sqrt(pow(sin_emission_angle, -41.0 / 20) - 1) * (169 * pow(kappa, -8) + 13.0 / 2500 * kappa - 1.0 / 200 + 47.0 / 200 / kappa) / sqrt_X;

        // ------------------------------------------------------------------------ Bridging function ------------------------------------------------------------------------ //

        double power_I = pow(-7.0 / 4 + 8.0 / 5 * kappa, -43.0 / 50);
        if (!isnan(Absorbtion_functions_high[I] / Absorbtion_functions_low[I])) {

            Absorbtion_functions[I] = Absorbtion_functions_low[I] * pow(1. + pow(Absorbtion_functions_high[I] / Absorbtion_functions_low[I], -power_I), -1.0 / power_I);

        }

        double power_Q = 7.0 / 5 * pow(kappa, -23.0 / 20);
        if (!isnan(Absorbtion_functions_high[Q] / Absorbtion_functions_low[Q])) {

            Absorbtion_functions[Q] = Absorbtion_functions_low[Q] * pow(1. + pow(Absorbtion_functions_high[Q] / Absorbtion_functions_low[V], -power_Q), -1.0 / power_Q);

        }

        Absorbtion_functions[U] = 0.0;

        double power_V = 61.0 / 50 * pow(kappa, -142.0 / 125) + 7.0 / 1000;
        if (!isnan(Absorbtion_functions_high[V] / Absorbtion_functions_low[V])) {

            Absorbtion_functions[V] = Absorbtion_functions_low[V] * pow(1. + pow(Absorbtion_functions_high[V] / Absorbtion_functions_low[V], -power_V), -1.0 / power_V);

        }

    }

   //Absorbtion_functions[I] = Absorbtion_functions_high[I];

}