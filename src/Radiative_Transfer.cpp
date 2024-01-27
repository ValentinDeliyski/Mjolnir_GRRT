#define _USE_MATH_DEFINES

#include "Enumerations.h"
#include "Constants.h"
#include "General_math_functions.h"

void Implicit_Trapezoid_Radiative_Transfer(double const Emission_Functions[STOKES_PARAM_NUM],
                                           double const Absorbtion_Functions[STOKES_PARAM_NUM],
                                           double const Faradey_Functions[STOKES_PARAM_NUM],
                                           double const step,
                                           double Stokes_Vector[STOKES_PARAM_NUM]) {

    /* ==================================================================================================|
    |                                                                                                    |
    |    This function applies the implicit trapezoidal rule to solve the radiative transfer             |
    |    equations. The equations are linear, so the method gives an explicit expression for             |
    |    the Stokes vector. The reference for this implementation is from the RAPTOR code/paper          |
    |    (the paper screwed up the explanation of the method tough - they skipped defining the           |
    |    variables x):                                                                                   |
    |    https://arxiv.org/pdf/2007.03045.pdf, https://github.com/tbronzwaer/raptor/tree/polarization    |
    |                                                                                                    |
    |    The variables u_ij and ell_ij are the components of the lower and upper triangular matricies    |
    |    in the LU decomposition of matrix A from (B.1)                                                  |
    |                                                                                                    |
    |=================================================================================================== */
    
    /* Here I define these for brevity, because it gets hairy otherwise - definitions follow (B.4) */

    const double*& alpha = Absorbtion_Functions;
    const double*& rho = Faradey_Functions;

    double u_11 = 1 + step * alpha[I] / 2;
    double u_12 = step * alpha[Q] / 2;
    double u_13 = 0;
    double u_14 = step * alpha[V] / 2;

    double ell_21 = step * alpha[Q] / 2 / u_11;

    double u_22 = 1 + step * alpha[I] / 2 - ell_21 * u_12;
    double u_23 = step * rho[V] / 2;
    double u_24 = -ell_21 * u_14;

    double ell_31 = 0;
    double ell_32 = -step * rho[V] / 2 / u_22;

    double u_33 = 1 + step * alpha[I] / 2 - ell_32 * u_23;
    double u_34 = step * rho[Q] / 2 - ell_32 * u_24;

    double ell_41 = step * alpha[V] / 2 / u_11;
    double ell_42 = -ell_41 * u_12 / u_22;
    double ell_43 = -(step * rho[Q] / 2 + ell_42 * u_23) / u_33;

    double u_44 = 1 + step * alpha[I] / 2 - ell_41 * u_14 - ell_42 * u_24 - ell_43 * u_34;

    double M[4][4] = { {alpha[I], alpha[Q], alpha[U], alpha[V]},
                       {alpha[Q], alpha[I],   rho[V],  -rho[U]},
                       {alpha[U],  -rho[V], alpha[I],   rho[Q]},
                       {alpha[V],   rho[U],  -rho[Q], alpha[I]} };

    double M_dot_Stokes_Vector[4]{};

    mat_vec_multiply_4D(M, Stokes_Vector, M_dot_Stokes_Vector);

    double b_1 = Stokes_Vector[I] + step / 2 * (2 * Emission_Functions[I] - M_dot_Stokes_Vector[I]);
    double b_2 = Stokes_Vector[Q] + step / 2 * (2 * Emission_Functions[Q] - M_dot_Stokes_Vector[Q]);
    double b_3 = Stokes_Vector[U] + step / 2 * (2 * Emission_Functions[U] - M_dot_Stokes_Vector[U]);
    double b_4 = Stokes_Vector[V] + step / 2 * (2 * Emission_Functions[V] - M_dot_Stokes_Vector[V]);

    double y_1 = b_1;
    double y_2 = b_2 - ell_21 * y_1;
    double y_3 = b_3 - ell_32 * y_2;
    double y_4 = b_4 - (ell_41 * y_1 + ell_42 * y_2 + ell_43 * y_3);

    double x_4 = y_4 / u_44;
    double x_3 = (y_3 - u_34 * x_4) / u_33;
    double x_2 = (y_2 - u_23 * x_3 - u_24 * x_4) / u_22;
    double x_1 = (y_1 - u_12 * x_2 - u_14 * x_4) / u_11;

    Stokes_Vector[I] = x_1;
    Stokes_Vector[Q] = x_2;
    Stokes_Vector[U] = x_3;
    Stokes_Vector[V] = x_4;

}

static void Get_radiative_transfer_matrix(double const absorbtion_functions[STOKES_PARAM_NUM],
                                          double const faradey_functions[STOKES_PARAM_NUM],
                                          double const step,
                                          double Transfer_Operator[STOKES_PARAM_NUM][STOKES_PARAM_NUM],
                                          double Integrated_Transfer_Operator[STOKES_PARAM_NUM][STOKES_PARAM_NUM]) {

    /* The reference for this implementation is from appendix D in https://arxiv.org/pdf/1602.03184.pdf, originally derived in https://doi.org/10.1007/BF00165988 */

    for (int row_idx = 0; row_idx <= STOKES_PARAM_NUM - 1; row_idx++) {

        for (int colum_idx = 0; colum_idx <= STOKES_PARAM_NUM - 1; colum_idx++) {

            Transfer_Operator[row_idx][colum_idx] = 0;
            Integrated_Transfer_Operator[row_idx][colum_idx] = 0;

        }

    }


    // Here I define a bunch of references, because its going to get hairy if I don't...
    const double*& alpha = absorbtion_functions;
    const double*& rho   = faradey_functions;

    /* These are the variables defined in D8 - D13, used in calculating the M matricies */

    double const alpha_squared = alpha[Q] * alpha[Q] +
        alpha[U] * alpha[U] +
        alpha[V] * alpha[V];

    double const rho_squared = rho[Q] * rho[Q] +
        rho[U] * rho[U] +
        rho[V] * rho[V];

    double const alpha_rho = alpha[Q] * rho[Q] +
        alpha[U] * rho[U] +
        alpha[V] * rho[V];

    // sigma is the sign of the variable alpha_rho
    double const sigma = copysign(1.0, alpha_rho);

    // These quantities can go ever so sligtly negative, which physically should not happen, but nmerically it does.
    // This breaks the sqrt() functions, and so guards have to be put in place
    double Theta = (alpha_squared - rho_squared) * (alpha_squared - rho_squared) / 4 + alpha_rho * alpha_rho;

    if (Theta >= 0.0) {

        Theta = 2 * sqrt(Theta);

    }
    else {

        Theta = std::numeric_limits<double>::min();

    }

    double Lambda[2] = { (Theta / 2 + (alpha_squared - rho_squared) / 2),
                         (Theta / 2 - (alpha_squared - rho_squared) / 2) };

    if (Lambda[0] >= 0.0 && Lambda[1] >= 0.0) {

        Lambda[0] = sqrt(Lambda[0]);
        Lambda[1] = sqrt(Lambda[1]);

    }
    else {

        Lambda[0] = 0.0;
        Lambda[1] = 0.0;

    }


    /* Thesse are used in the "scaling factors" infront of the M matricies */

    double const exp_I = exp(-alpha[I] * step);

    double const cosh_term = cosh(Lambda[0] * step);
    double const cos_term = cos(Lambda[1] * step);

    double const sinh_term = sinh(Lambda[0] * step);
    double const sin_term = sin(Lambda[1] * step);

    /* ========================== M_1 Matrix calculation ========================== */

    double const M_1_scale_factor = exp_I * (cosh_term + cos_term) / 2;

    double const M_1[4][4] = { {1.0, 0.0, 0.0, 0.0},
                               {0.0, 1.0, 0.0, 0.0},
                               {0.0, 0.0, 1.0, 0.0},
                               {0.0, 0.0, 0.0, 1.0} };

    /* ========================== M_2 Matrix calculation ========================== */

    double M_2_scale_factor{};

    if (Theta > std::numeric_limits<double>::min()) {

        M_2_scale_factor = -exp_I * sin_term / Theta;
    }

    double const M_2[4][4] = { {                         0,                          (Lambda[1] * alpha[Q] - sigma * Lambda[0] * rho[Q]), (Lambda[1] * alpha[U] - sigma * Lambda[0] * rho[U]), (Lambda[1] * alpha[V] - sigma * Lambda[0] * rho[V])},
                               {(Lambda[1] * alpha[Q] - sigma * Lambda[0] * rho[Q]),                           0,                          (sigma * Lambda[0] * alpha[V] + Lambda[1] * rho[V]), (-sigma * Lambda[0] * alpha[U] - Lambda[1] * rho[U])},
                               {(Lambda[1] * alpha[U] - sigma * Lambda[0] * rho[U]), (-sigma * Lambda[0] * alpha[V] - Lambda[1] * rho[V]),                           0,                          (sigma * Lambda[0] * alpha[Q] + Lambda[1] * rho[Q])},
                               {(Lambda[1] * alpha[V] - sigma * Lambda[0] * rho[V]), (sigma * Lambda[0] * alpha[U] + Lambda[1] * rho[U]), (-sigma * Lambda[0] * alpha[Q] - Lambda[1] * rho[Q]),                           0                         } };

    /* ========================== M_3 Matrix calculation ========================== */

    double M_3_scale_factor{};

    if (Theta > std::numeric_limits<double>::min()) {

        M_3_scale_factor = -exp_I * sinh_term / Theta;
    }

    double const M_3[4][4] = { {						 0,							 (Lambda[0] * alpha[Q] + sigma * Lambda[1] * rho[Q]), (Lambda[0] * alpha[U] + sigma * Lambda[1] * rho[Q]), (Lambda[0] * alpha[V] + sigma * Lambda[1] * rho[V])},
                               {(Lambda[0] * alpha[Q] + sigma * Lambda[1] * rho[Q]),	 		               0,                          (-sigma * Lambda[1] * alpha[V] + Lambda[0] * rho[V]), (sigma * Lambda[1] * alpha[U] - Lambda[0] * rho[U])},
                               {(Lambda[0] * alpha[U] + sigma * Lambda[1] * rho[U]), (sigma * Lambda[1] * alpha[V] - Lambda[0] * rho[V]),	                         0,	                         (-sigma * Lambda[1] * alpha[Q] + Lambda[0] * rho[Q])},
                               {(Lambda[0] * alpha[V] + sigma * Lambda[1] * rho[V]), (-sigma * Lambda[1] * alpha[U] + Lambda[0] * rho[U]), (sigma * Lambda[1] * alpha[Q] - Lambda[0] * rho[Q]),                           0                         } };

    /* ========================== M_4 Matrix calculation ========================== */

    double M_4_scale_factor{};

    if (Theta > std::numeric_limits<double>::min()) {

        M_4_scale_factor = exp_I * (cosh_term - cos_term) / Theta;
    }

    double const M_4[4][4] = { {   (alpha_squared + rho_squared) / 2,                      (alpha[V] * rho[U] - alpha[U] * rho[V]),                                     (alpha[Q] * rho[V] - alpha[V] * rho[Q]),                                     (alpha[U] * rho[Q] - alpha[Q] * rho[U])},
                               {(alpha[U] * rho[V] - alpha[V] * rho[U]), (alpha[Q] * alpha[Q] + rho[Q] * rho[Q] - (alpha_squared + rho_squared) / 2),                   (alpha[Q] * alpha[U] + rho[Q] * rho[U]),                                     (alpha[V] * alpha[Q] + rho[V] * rho[Q])},
                               {(alpha[V] * rho[Q] - alpha[Q] * rho[V]),                   (alpha[Q] * alpha[U] + rho[Q] * rho[U]),                   (alpha[U] * alpha[U] + rho[U] * rho[U] - (alpha_squared + rho_squared) / 2),                   (alpha[U] * alpha[V] + rho[U] * rho[V])},
                               {(alpha[Q] * rho[U] - alpha[U] * rho[Q]),                   (alpha[V] * alpha[Q] + rho[V] * rho[Q]),                                     (alpha[U] * alpha[V] + rho[U] * rho[V]),                   (alpha[V] * alpha[V] + rho[V] * rho[V] - (alpha_squared + rho_squared) / 2)} };

    /* ========================== This is the formal operator O(s,s') - the solution to D1 ========================== */

    for (int row_idx = 0; row_idx <= STOKES_PARAM_NUM - 1; row_idx++) {

        for (int colum_idx = 0; colum_idx <= STOKES_PARAM_NUM - 1; colum_idx++) {

            Transfer_Operator[row_idx][colum_idx] = M_1_scale_factor * M_1[row_idx][colum_idx] +
                M_2_scale_factor * M_2[row_idx][colum_idx] +
                M_3_scale_factor * M_3[row_idx][colum_idx] +
                M_4_scale_factor * M_4[row_idx][colum_idx];

        }

    }

    /* ========================== The intergral of O(s,s') for constant M matricies ========================== */

    /* This part of the implementation is adapted from equation (24) of https://academic.oup.com/mnras/article/475/1/43/4712230 */

    bool math_guard_1 = fabs(alpha[I] * alpha[I] - Lambda[0] * Lambda[0]) > std::numeric_limits<double>::min();
    bool math_guard_2 = fabs(alpha[I] * alpha[I] + Lambda[1] * Lambda[1]) > std::numeric_limits<double>::min();

    if (math_guard_1 && math_guard_2) {

        double f_1 = 1.0 / (alpha[I] * alpha[I] - Lambda[0] * Lambda[0]);
        double f_2 = 1.0 / (alpha[I] * alpha[I] + Lambda[1] * Lambda[1]);

        for (int row_idx = 0; row_idx <= STOKES_PARAM_NUM - 1; row_idx++) {

            for (int colum_idx = 0; colum_idx <= STOKES_PARAM_NUM - 1; colum_idx++) {


                Integrated_Transfer_Operator[row_idx][colum_idx] = -Lambda[0] * f_1 * M_3[row_idx][colum_idx] + alpha[I] * f_1 / 2 * (M_1[row_idx][colum_idx] + M_4[row_idx][colum_idx]) +
                                                                    Lambda[1] * f_2 * M_2[row_idx][colum_idx] + alpha[I] * f_2 / 2 * (M_1[row_idx][colum_idx] - M_4[row_idx][colum_idx]) -
                                                                    exp_I * ((-Lambda[0] * f_1 * M_3[row_idx][colum_idx] +  alpha[I] * f_1 / 2 * (M_1[row_idx][colum_idx] + M_4[row_idx][colum_idx])) * cosh_term +
                                                                             (-Lambda[1] * f_2 * M_2[row_idx][colum_idx] +  alpha[I] * f_2 / 2 * (M_1[row_idx][colum_idx] - M_4[row_idx][colum_idx])) * cos_term +
                                                                             ( -alpha[I] * f_2 * M_2[row_idx][colum_idx] - Lambda[1] * f_2 / 2 * (M_1[row_idx][colum_idx] - M_4[row_idx][colum_idx])) * sin_term -
                                                                             (  alpha[I] * f_1 * M_3[row_idx][colum_idx] - Lambda[0] * f_1 / 2 * (M_1[row_idx][colum_idx] + M_4[row_idx][colum_idx])) * sinh_term);

            }

        }

    }

}

void Analytic_Radiative_Transfer(double const emission_functions[STOKES_PARAM_NUM],
                                 double const absorbtion_functions[STOKES_PARAM_NUM],
                                 double const faradey_functions[STOKES_PARAM_NUM],
                                 double const step,
                                 double Intensity[STOKES_PARAM_NUM]){

    double transfer_operator[4][4]{};
    double integrated_transfer_operator[4][4];

    Get_radiative_transfer_matrix(absorbtion_functions, faradey_functions, step, transfer_operator, integrated_transfer_operator);

    double Transfered_emission_vector[STOKES_PARAM_NUM]{};
    double temp_Intensity[STOKES_PARAM_NUM]{};

    // Placeholder vector for use in the mat_vec_multiply_4D() function
    for (int index = 0; index <= STOKES_PARAM_NUM - 1; index++) {

        temp_Intensity[index] = Intensity[index];

    }

    mat_vec_multiply_4D(integrated_transfer_operator, emission_functions, Transfered_emission_vector);
    mat_vec_multiply_4D(transfer_operator, temp_Intensity, Intensity);

    for (int index = 0; index <= STOKES_PARAM_NUM - 1; index++) {

        Intensity[index] += Transfered_emission_vector[index];

    }

}