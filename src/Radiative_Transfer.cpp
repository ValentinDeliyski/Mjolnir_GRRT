#include "Radiative_Transfer.h"

void static get_M_matrix(double M_matrix[e_Stokes_param_num][e_Stokes_param_num],
                         const double* absorbtion_function, 
                         const double* faradey_function){

    const double*& alpha = absorbtion_function;
    const double*& rho = faradey_function;

    M_matrix[0][0] = M_matrix[1][1] = M_matrix[2][2] = M_matrix[3][3] = alpha[I];
    M_matrix[0][1] = M_matrix[1][0] = alpha[Q];
    M_matrix[0][2] = M_matrix[2][0] = alpha[U];
    M_matrix[0][3] = M_matrix[3][0] = alpha[V];
    M_matrix[1][2] =  rho[V];
    M_matrix[2][1] = -rho[V];
    M_matrix[1][3] = -rho[U];
    M_matrix[3][1] =  rho[U];
    M_matrix[2][3] =  rho[Q];
    M_matrix[3][2] = -rho[Q];

}

void Implicit_Trapezoid_Radiative_Transfer(double* const Emission_Functions,
                                           double* const Absorbtion_Functions,
                                           double* const Faradey_Functions,
                                           const double step,
                                           double* const Stokes_Vector) {

    /* ==================================================================================================|
    |                                                                                                    |
    |    This function applies the implicit trapezoidal rule to solve the radiative transfer             |
    |    equations. The equations are linear, so the method gives an explicit expression for             |
    |    the Stokes vector. The reference for this implementation is from the RAPTOR code/paper          |
    |    (the paper screwed up the explanation of the method - they skipped defining the variables x):   |                                                                              
    |    https://arxiv.org/pdf/2007.03045.pdf, https://github.com/tbronzwaer/raptor/tree/polarization    |
    |                                                                                                    |
    |    The variables u_ij and ell_ij are the components of the lower and upper triangular matricies    |
    |    in the LU decomposition of matrix A from (B.1)                                                  |
    |                                                                                                    |
    |=================================================================================================== */
    
    /* Here I define these for brevity, because it gets hairy otherwise - definitions follow (B.4) */

    const double* const& alpha = Absorbtion_Functions;
    const double* const& rho = Faradey_Functions;

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

    double M_matrix[4][4]{};

    get_M_matrix(M_matrix, Absorbtion_Functions, Faradey_Functions);

    double M_dot_Stokes_Vector[4]{};

    mat_vec_multiply_4D(M_matrix, Stokes_Vector, M_dot_Stokes_Vector);

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

static Return_Values Get_radiative_transfer_matrix(double* const absorbtion_functions,
                                                   double* const faradey_functions,
                                                   const double step,
                                                   double Transfer_Operator[e_Stokes_param_num][e_Stokes_param_num],
                                                   double Integrated_Transfer_Operator[e_Stokes_param_num][e_Stokes_param_num]) {

    /* The reference for this implementation is from appendix D in https://arxiv.org/pdf/1602.03184.pdf, originally derived in https://doi.org/10.1007/BF00165988 */

    for (int row_idx = 0; row_idx < e_Stokes_param_num ; row_idx++) {

        for (int colum_idx = 0; colum_idx < e_Stokes_param_num ; colum_idx++) {

            Transfer_Operator[row_idx][colum_idx] = 0;
            Integrated_Transfer_Operator[row_idx][colum_idx] = 0;

        }

    }


    // Here I define a bunch of references, because its going to get hairy if I don't...
    const double* const& alpha = absorbtion_functions;
    const double* const& rho   = faradey_functions;

    /* These are the variables defined in D8 - D13, used in calculating the M matricies */

    const double alpha_squared = alpha[Q] * alpha[Q] + alpha[U] * alpha[U] + alpha[V] * alpha[V];
    const double rho_squared = rho[Q] * rho[Q] + rho[U] * rho[U] + rho[V] * rho[V];
    const double alpha_rho = alpha[Q] * rho[Q] + alpha[U] * rho[U] + alpha[V] * rho[V];

    // sigma is the sign of the variable alpha_rho
    const double sigma = copysign(1.0, alpha_rho);

    // This is identically zero if no polarization is included
    double Theta = 2 * sqrt((alpha_squared - rho_squared) * (alpha_squared - rho_squared) / 4 + alpha_rho * alpha_rho) + 1e-40;

    if (isnan(1.0 / Theta) || isinf(1.0 / Theta)) { return ERROR; }

    double Lambda[2] = { sqrt(Theta / 2 + (alpha_squared - rho_squared) / 2) + 1e-40,
                         sqrt(Theta / 2 - (alpha_squared - rho_squared) / 2) + 1e-40 };

    if (isnan(Lambda[0]) || isnan(Lambda[1])) { return ERROR; }

    /* Thesse are used in the "scaling factors" infront of the M matricies */
    const double exp_I = exp(-alpha[I] * step);

    const double cosh_term = cosh(Lambda[0] * step);
    const double cos_term = cos(Lambda[1] * step);

    const double sinh_term = sinh(Lambda[0] * step);
    const double sin_term = sin(Lambda[1] * step);

    /* ========================== M_1 Matrix calculation ========================== */

    const double M_1_scale_factor = exp_I * (cosh_term + cos_term) / 2;

    const double M_1[4][4] = { {1.0, 0.0, 0.0, 0.0},
                               {0.0, 1.0, 0.0, 0.0},
                               {0.0, 0.0, 1.0, 0.0},
                               {0.0, 0.0, 0.0, 1.0} };

    /* ========================== M_2 Matrix calculation ========================== */

    const double M_2_scale_factor = -exp_I * sin_term / Theta;

    const double M_2[4][4] = { {                         0,                          (Lambda[1] * alpha[Q] - sigma * Lambda[0] * rho[Q]), ( Lambda[1] * alpha[U] - sigma * Lambda[0] * rho[U]), ( Lambda[1] * alpha[V] - sigma * Lambda[0] * rho[V])},
                               {(Lambda[1] * alpha[Q] - sigma * Lambda[0] * rho[Q]),                           0,                         ( sigma * Lambda[0] * alpha[V] + Lambda[1] * rho[V]), (-sigma * Lambda[0] * alpha[U] - Lambda[1] * rho[U])},
                               {(Lambda[1] * alpha[U] - sigma * Lambda[0] * rho[U]), (-sigma * Lambda[0] * alpha[V] - Lambda[1] * rho[V]),                           0,                         ( sigma * Lambda[0] * alpha[Q] + Lambda[1] * rho[Q])},
                               {(Lambda[1] * alpha[V] - sigma * Lambda[0] * rho[V]), ( sigma * Lambda[0] * alpha[U] + Lambda[1] * rho[U]), (-sigma * Lambda[0] * alpha[Q] - Lambda[1] * rho[Q]),                           0                         } };

    /* ========================== M_3 Matrix calculation ========================== */

    const double M_3_scale_factor = -exp_I * sinh_term / Theta;

    const double M_3[4][4] = { {						 0,							 ( Lambda[0] * alpha[Q] + sigma * Lambda[1] * rho[Q]), ( Lambda[0] * alpha[U] + sigma * Lambda[1] * rho[U]), ( Lambda[0] * alpha[V] + sigma * Lambda[1] * rho[V])},
                               {(Lambda[0] * alpha[Q] + sigma * Lambda[1] * rho[Q]),	 		               0,                          (-sigma * Lambda[1] * alpha[V] + Lambda[0] * rho[V]), ( sigma * Lambda[1] * alpha[U] - Lambda[0] * rho[U])},
                               {(Lambda[0] * alpha[U] + sigma * Lambda[1] * rho[U]), ( sigma * Lambda[1] * alpha[V] - Lambda[0] * rho[V]),	                         0,	                         (-sigma * Lambda[1] * alpha[Q] + Lambda[0] * rho[Q])},
                               {(Lambda[0] * alpha[V] + sigma * Lambda[1] * rho[V]), (-sigma * Lambda[1] * alpha[U] + Lambda[0] * rho[U]), ( sigma * Lambda[1] * alpha[Q] - Lambda[0] * rho[Q]),                           0                          } };

    /* ========================== M_4 Matrix calculation ========================== */

    const double M_4_scale_factor = exp_I * (cosh_term - cos_term) / Theta;

    const double M_4[4][4] = { {   (alpha_squared + rho_squared) / 2,                      (alpha[V] * rho[U] - alpha[U] * rho[V]),                                     (alpha[Q] * rho[V] - alpha[V] * rho[Q]),                                     (alpha[U] * rho[Q] - alpha[Q] * rho[U])},
                               {(alpha[U] * rho[V] - alpha[V] * rho[U]), (alpha[Q] * alpha[Q] + rho[Q] * rho[Q] - (alpha_squared + rho_squared) / 2),                   (alpha[Q] * alpha[U] + rho[Q] * rho[U]),                                     (alpha[V] * alpha[Q] + rho[V] * rho[Q])},
                               {(alpha[V] * rho[Q] - alpha[Q] * rho[V]),                   (alpha[Q] * alpha[U] + rho[Q] * rho[U]),                   (alpha[U] * alpha[U] + rho[U] * rho[U] - (alpha_squared + rho_squared) / 2),                   (alpha[U] * alpha[V] + rho[U] * rho[V])},
                               {(alpha[Q] * rho[U] - alpha[U] * rho[Q]),                   (alpha[V] * alpha[Q] + rho[V] * rho[Q]),                                     (alpha[U] * alpha[V] + rho[U] * rho[V]),                   (alpha[V] * alpha[V] + rho[V] * rho[V] - (alpha_squared + rho_squared) / 2)} };

    /* ========================== This is the formal operator O(s,s') - the solution to D1 ========================== */

    for (int row_idx = 0; row_idx < e_Stokes_param_num ; row_idx++) {

        for (int colum_idx = 0; colum_idx < e_Stokes_param_num ; colum_idx++) {

            Transfer_Operator[row_idx][colum_idx] = M_1_scale_factor * M_1[row_idx][colum_idx] +
                                                    M_2_scale_factor * M_2[row_idx][colum_idx] +
                                                    M_3_scale_factor * M_3[row_idx][colum_idx] +
                                                    M_4_scale_factor * M_4[row_idx][colum_idx];

        }

    }

    /* ========================== The intergral of O(s,s') for constant M matricies ========================== */

    /* This part of the implementation is adapted from equation (24) of https://academic.oup.com/mnras/article/475/1/43/4712230 */

    const double f_1 = 1.0 / (alpha[I] * alpha[I] - Lambda[0] * Lambda[0]);
    const double f_2 = 1.0 / (alpha[I] * alpha[I] + Lambda[1] * Lambda[1]);

    if (isinf(f_1) || isinf(f_2)) { return ERROR; }

    for (int row_idx = 0; row_idx < e_Stokes_param_num ; row_idx++) {

        for (int colum_idx = 0; colum_idx < e_Stokes_param_num ; colum_idx++) {


            Integrated_Transfer_Operator[row_idx][colum_idx] = -Lambda[0] * f_1 * M_3[row_idx][colum_idx] + alpha[I] * f_1 / 2. * (M_1[row_idx][colum_idx] + M_4[row_idx][colum_idx]) +
                                                               -Lambda[1] * f_2 * M_2[row_idx][colum_idx] + alpha[I] * f_2 / 2. * (M_1[row_idx][colum_idx] - M_4[row_idx][colum_idx]) -
                                                                exp_I * ((-Lambda[0] * f_1 * M_3[row_idx][colum_idx] +  alpha[I] * f_1 / 2 * (M_1[row_idx][colum_idx] + M_4[row_idx][colum_idx])) * cosh_term +
                                                                         (-Lambda[1] * f_2 * M_2[row_idx][colum_idx] +  alpha[I] * f_2 / 2 * (M_1[row_idx][colum_idx] - M_4[row_idx][colum_idx])) * cos_term +
                                                                         ( -alpha[I] * f_2 * M_2[row_idx][colum_idx] - Lambda[1] * f_2 / 2 * (M_1[row_idx][colum_idx] - M_4[row_idx][colum_idx])) * sin_term -
                                                                         (  alpha[I] * f_1 * M_3[row_idx][colum_idx] - Lambda[0] * f_1 / 2 * (M_1[row_idx][colum_idx] + M_4[row_idx][colum_idx])) * sinh_term);

        }

    }

    return OK;

}

static void run_Decoupled_Radiative_Transfer(const double* const Emission_Functions,
                                             const double* const Absorbtion_Functions,
                                             const double step, 
                                             double* const Stokes_Vector) {

    if (Absorbtion_Functions[I] > 0) {

        for (int idx = 0; idx < e_Stokes_param_num; idx++) {

            Stokes_Vector[idx] = Stokes_Vector[idx] * exp(-Absorbtion_Functions[I] * step) + (1 - exp(-Absorbtion_Functions[I] * step)) * Emission_Functions[idx] / Absorbtion_Functions[I];

        }

    }
    else {

        for (int idx = 0; idx < e_Stokes_param_num; idx++) {

            Stokes_Vector[idx] += Emission_Functions[idx] * step;

        }

    }

}

void Analytic_Radiative_Transfer(double* const Emission_Functions,
                                 double* const Absorbtion_Functions,
                                 double* const Faradey_Functions,
                                 const double step,
                                 double* const Stokes_Vector){

    double transfer_operator[4][4]{};
    double integrated_transfer_operator[4][4];

    if (OK != Get_radiative_transfer_matrix(Absorbtion_Functions, Faradey_Functions, step, transfer_operator, integrated_transfer_operator)) {

        run_Decoupled_Radiative_Transfer(Emission_Functions, Absorbtion_Functions, step, Stokes_Vector);

        return;

    }

    double Transfered_emission_vector[e_Stokes_param_num]{};
    double temp_Stokes_Vector[e_Stokes_param_num]{};

    memcpy(temp_Stokes_Vector, Stokes_Vector, 4 * sizeof(double));

    mat_vec_multiply_4D(integrated_transfer_operator, Emission_Functions, Transfered_emission_vector);
    mat_vec_multiply_4D(transfer_operator, temp_Stokes_Vector, Stokes_Vector);

    for (int index = 0; index < e_Stokes_param_num; index++) {

        Stokes_Vector[index] += Transfered_emission_vector[index];

    }

}

void static Radiative_transfer_RHS(const double* const Emission_Functions, 
                                   const double* const Absorbtion_Functions, 
                                   const double* const Faradey_Functions, 
                                   const double* const Stokes_Vector,
                                   double* const RHS) {

    double M_matrix[4][4]{};
    get_M_matrix(M_matrix, Absorbtion_Functions, Faradey_Functions);

    double M_dot_Stokes[4]{};
    mat_vec_multiply_4D(M_matrix, Stokes_Vector, M_dot_Stokes);

    for (int idx = 0; idx < e_Stokes_param_num ; idx++) {

        RHS[idx] = Emission_Functions[idx] - M_dot_Stokes[idx];

    }

}

void RK5_radiative_transfer(double* const Emission_Functions,
                            double* const Absorbtion_Functions,
                            double* const Faradey_Functions,
                            double* const State_Vector,
                            const Simulation_Context_type* p_Sim_Context,
                            double* const Stokes_Vector) {

    double RHS[Nyström_size * e_Stokes_param_num]{};
    double EOM[Nyström_size * e_Dynamic_state_size]{};

    double Temp_Stokes_Vector[e_Stokes_param_num]{}, Temp_State_Vector[e_Dynamic_state_size]{};

    for (int RK5_stage = 0; RK5_stage < Nyström_size; RK5_stage++) {

        memcpy(Temp_Stokes_Vector, Stokes_Vector, e_Stokes_param_num * sizeof(double));
        memcpy(Temp_State_Vector, State_Vector, e_Dynamic_state_size * sizeof(double));

        for (int derivative_indexer = 0; derivative_indexer < RK5_stage; derivative_indexer++) {

            for (int idx = 0; idx < 4; idx++) {

                Temp_Stokes_Vector[idx] += Nyström_Deriv_coeffs[RK5_stage][derivative_indexer] * RHS[idx + derivative_indexer * e_Stokes_param_num] * State_Vector[e_step];

            }

            for (int idx = 0; idx < e_Dynamic_state_size; idx++) {

                Temp_State_Vector[idx] += Nyström_Deriv_coeffs[RK5_stage][derivative_indexer] * EOM[idx + derivative_indexer * e_Dynamic_state_size] * State_Vector[e_step];

            }
        }

        Transfer_functions_type Total_Transfer_Functions{};

        for (int emission_medium = Disk; emission_medium <= Hotspot; emission_medium++) {

            Transfer_functions_type Temp_Transfer_functions{};

            p_Sim_Context->p_Emission_Model->get_radiative_transfer_functions(Temp_State_Vector,
                                                                              p_Sim_Context,
                                                                              static_cast<Emission_medium_enums>(emission_medium),
                                                                              &Temp_Transfer_functions);

            add_vectors(Temp_Transfer_functions.Absorbtion_functions, Total_Transfer_Functions.Absorbtion_functions, e_Stokes_param_num, Total_Transfer_Functions.Absorbtion_functions);
            add_vectors(Temp_Transfer_functions.Emission_functions, Total_Transfer_Functions.Emission_functions, e_Stokes_param_num, Total_Transfer_Functions.Emission_functions);
            add_vectors(Temp_Transfer_functions.Faradey_functions, Total_Transfer_Functions.Faradey_functions, e_Stokes_param_num, Total_Transfer_Functions.Faradey_functions);

        }

        Radiative_transfer_RHS(Total_Transfer_Functions.Emission_functions, 
                               Total_Transfer_Functions.Absorbtion_functions,
                               Total_Transfer_Functions.Faradey_functions, 
                               Temp_Stokes_Vector, 
                               RHS + RK5_stage * e_Stokes_param_num);

        p_Sim_Context->p_Spacetime->get_EOM(Temp_State_Vector, EOM + RK5_stage * e_Dynamic_state_size);

    }

    for (int idx = 0; idx < e_Stokes_param_num; idx++) {

        for (int deriv_idx = 0; deriv_idx < Nyström_size; deriv_idx++) {

            Stokes_Vector[idx] += State_Vector[e_step] * Nyström_Coeff_sol[deriv_idx] * RHS[idx + deriv_idx * e_Stokes_param_num];

        }
    }

}