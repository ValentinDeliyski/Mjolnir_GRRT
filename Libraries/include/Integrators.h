#pragma once
#include "General_math_functions.h"
#include "Enumerations.h"
#include "Spacetimes.h"
#include "Constants.h"
#include "Structs.h"
#include "gsl\gsl_multiroots.h"
#include <functional>
#include <cmath>

class Step_controller_class {

public:

    Step_controller_class(const Integrator_parameters_type Integrator_parameters);

    //! Updates the integration step, based on the previous State Error estimates, and the current State Vector.
    /*! Updates the integration step, based on the previous State Error estimates, and the current State Vector.
     *   Currently the following step controllers are implemented. The reference is https://arxiv.org/pdf/1806.08693:
     *      1) PID controller
     *      2) Gustafsson controller
     *
     *   \param [in] State_Vector - Pointer to the array that holds the photon State Vector.
     *   \return Nothing
     */
    void update_step(const double* State_Vector);

    void update_state_errors(const double* State_Vector, const double* State_Error_Vector);

    Integrator_parameters_type Parameters;

    double step;
    double previous_step;

    double current_err;
    double prev_err;
    double sec_prev_err;

};

class Integrator_class {

private:

    /* ====================================================== The adaptive RK7(8) Butcher table ====================================================== */

    // The coefficients are from https://ntrs.nasa.gov/api/citations/19680027281/downloads/19680027281.pdf, table X
    double RK78_Coeff_deriv[RK78_size][RK78_size - 1] =
    {
        {       0,         0,         0,          0,            0,            0,          0,          0,        0,         0,     0,  0  },
        {    2. / 27,      0,         0,          0,            0,            0,          0,          0,        0,         0,     0,  0  },
        {    1. / 36,   1. / 12,      0,          0,            0,            0,          0,          0,        0,         0,     0,  0  },
        {    1. / 24,      0,      1. / 8,        0,            0,            0,          0,          0,        0,         0,     0,  0  },
        {    5. / 12,      0,    -25. / 16,   25. / 16,         0,            0,          0,          0,        0,         0,     0,  0  },
        {    1. / 20,      0,         0,       1. / 4,       1. / 5,          0,          0,          0,        0,         0,     0,  0  },
        {  -25. / 108,     0,         0,     125. / 108,   -65. / 27,    125. / 54,       0,          0,        0,         0,     0,  0  },
        {   31. / 300,     0,         0,          0,        61. / 225,    -2. / 9,    13. / 900,      0,        0,         0,     0,  0  },
        {      2.,         0,         0,     -53. / 6,     704. / 45,   -107. / 9,    67. / 90,       3.,       0,         0,     0,  0  },
        {  -91. / 108,     0,         0,      23. / 108,  -976. / 135,   311. / 54,  -19. / 60,   17. / 6,  -1. / 12,      0,     0,  0  },
        { 2383. / 4100,    0,         0,    -341. / 164,  4496. / 1025, -301. / 82, 2133. / 4100, 45. / 82, 45. / 164, 18. / 41,  0,  0  },
        {    3. / 205,     0,         0,          0,            0,        -6. / 41,   -3. / 205,  -3. / 41,  3. / 41,   6. / 41,  0,  0  },
        {-1777. / 4100,    0,         0,    -341. / 164,  4496. / 1025, -289. / 82, 2193. / 4100, 51. / 82, 33. / 164, 12. / 41,  0,  1. }
    };

    double RK78_Coeff_sol[13]      = { 41. / 840, 0, 0, 0, 0, 34. / 105, 9. / 35, 9. / 35, 9. / 280, 9. / 280, 41. / 840,     0,         0     };
    double RK78_Coeff_test_sol[13] = {     0,     0, 0, 0, 0, 34. / 105, 9. / 35, 9. / 35, 9. / 280, 9. / 280,     0,     41. / 840, 41. / 840 };

    double RK78_Stability_pol_coeffs[13] = { 1., 0.5, 0.16666667, 0.04166667, 0.00833333, 0.00138889, 0.00019841, 2.31653715e-05, 2.36714395e-06, 5.18294488e-08, -4.31912073e-08, 0., 0. };


    /* ====================================================== The adaptive DP RK8(7) Butcher table ====================================================== */

    // The coefficients are from https://ntrs.nasa.gov/api/citations/19680027281/downloads/19680027281.pdf, table X
    double RK78DP_Coeff_deriv[RK78_size][RK78_size - 1] =
    {
        {       0,         0,         0,          0,            0,            0,          0,          0,        0,         0,     0,  0  },
        {    2. / 27,      0,         0,          0,            0,            0,          0,          0,        0,         0,     0,  0  },
        {    1. / 36,   1. / 12,      0,          0,            0,            0,          0,          0,        0,         0,     0,  0  },
        {    1. / 24,      0,      1. / 8,        0,            0,            0,          0,          0,        0,         0,     0,  0  },
        {    5. / 12,      0,    -25. / 16,   25. / 16,         0,            0,          0,          0,        0,         0,     0,  0  },
        {    1. / 20,      0,         0,       1. / 4,       1. / 5,          0,          0,          0,        0,         0,     0,  0  },
        {  -25. / 108,     0,         0,     125. / 108,   -65. / 27,    125. / 54,       0,          0,        0,         0,     0,  0  },
        {   31. / 300,     0,         0,          0,        61. / 225,    -2. / 9,    13. / 900,      0,        0,         0,     0,  0  },
        {      2.,         0,         0,     -53. / 6,     704. / 45,   -107. / 9,    67. / 90,       3.,       0,         0,     0,  0  },
        {  -91. / 108,     0,         0,      23. / 108,  -976. / 135,   311. / 54,  -19. / 60,   17. / 6,  -1. / 12,      0,     0,  0  },
        { 2383. / 4100,    0,         0,    -341. / 164,  4496. / 1025, -301. / 82, 2133. / 4100, 45. / 82, 45. / 164, 18. / 41,  0,  0  },
        {    3. / 205,     0,         0,          0,            0,        -6. / 41,   -3. / 205,  -3. / 41,  3. / 41,   6. / 41,  0,  0  },
        {-1777. / 4100,    0,         0,    -341. / 164,  4496. / 1025, -289. / 82, 2193. / 4100, 51. / 82, 33. / 164, 12. / 41,  0,  1. }
    };

    double RK78DP_Coeff_sol[13] = { 41. / 840, 0, 0, 0, 0, 34. / 105, 9. / 35, 9. / 35, 9. / 280, 9. / 280, 41. / 840,     0,         0 };
    double RK78DP_Coeff_test_sol[13] = { 0,     0, 0, 0, 0, 34. / 105, 9. / 35, 9. / 35, 9. / 280, 9. / 280,     0,     41. / 840, 41. / 840 };

    double RK78DP_Stability_pol_coeffs[13] = { 1., 0.5, 0.16666667, 0.04166667, 0.00833333, 0.00138889, 0.00019841, 2.31653715e-05, 2.36714395e-06, 5.18294488e-08, -4.31912073e-08, 0., 0. };

    /* ====================================================== The fixed step RK5 Butcher table ====================================================== */

    // The coefficients are from this tage https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods
    double Nyström_Deriv_coeffs[Nyström_size][Nyström_size] =
    { 
        {   0.,      0.,         0.,       0.,   0., 0.},
        {1. / 3,     0.,         0.,       0.,   0., 0.},
        {4. / 25, 6. / 25,       0.,       0.,   0., 0.},
        {1. / 4,    -3.,     15. / 4,      0.,   0., 0.},
        {2. / 27, 10. / 9,  -50. / 81,  8. / 81, 0., 0.},
        {2. / 25, 12. / 25,   2. / 15,  8. / 75, 0., 0.} 
    };

    double Nyström_Coeff_sol[6] = { 23. / 192 , 0., 125. / 192, 0., -27. / 64, 125. / 192 };

    /* ====================================================== The fixed step Adams-Bashforth coefficients ====================================================== */
    
    // The first index corresponds to the order of the method, while the second index (idx_2) corresponds to the coefficient infront of the f_{n + idx_2) term.
    // The reference for this is table 1 from https://math.iit.edu/~fass/478578_Chapter_2.pdf.
    double AB_history_coefficients[5][5] =
    {
        {      1.,            0.,          0.,           0.,          0.,   },
        {    -0.5,         3. / 2,         0.,           0.,          0.,   },
        {   5. / 12,     -16. / 12,    23. / 12,         0.,          0.,   },
        {  -9. / 24,      37. / 24,   -59. / 24,     55. / 24,        0.,   },
        { 251. / 720,  -1274. / 720, 2616. / 720, -2774. / 720, 1901. / 720,}
    };

    // The first index corresponds to the order of the method, while the second index (idx_2) corresponds to the coefficient infront of the y_{n + idx_2) term.
    // The reference for this is table 3 from https://math.iit.edu/~fass/478578_Chapter_2.pdf.
    double BDF_history_coefficients[5][6] =
    {
            {     1.,         0.,          0.,          0.,          0.,   },
            {  -1. / 3,     4. / 2,        0.,          0.,          0.,   },
            {   2. / 11,   -9. / 11,    18. / 11,       0.,          0.,   },
            {  -3. / 25,   16. / 25,   -36. / 25,    48. / 25,       0.,   },
            {  12. / 137, -75. / 137,  200. / 137, -300. / 137, 300. / 137,}
    };

    // The index corresponds to the order of the method, while value at that index to the coefficient infront of the f_{n + idx + 1) term.
    // The reference for this is table 3 from https://math.iit.edu/~fass/478578_Chapter_2.pdf.
    double BDF_RHS_coefficients[5] = { 1., 2. / 3, 6. / 11, 12. / 25, 60. / 137 };

    bool In_stiff_region;
    bool Implicit_method_init_status;

    Geodesic_Integrator_enums e_Active_integrator;
    Initial_conditions_type* p_Init_conditions;

    gsl_vector* gsl_trial_State_Vector;
    gsl_multiroot_fsolver *Root_finder;
    gsl_multiroot_function Function_to_solve;

    void Run_RK78();

    Stability_return_type Check_method_stability();

    void Check_integration_complete_status();

    void Init_BDF();

    int Run_BDF();

    void Update_ray_log(const double* const New_State_vector);


public:

    Integrator_class(const Simulation_Context_type* const p_Sim_Context, Results_type* p_Ray_results);

    s_Ray_log_type* p_Ray_log_struct;
    RHS_wrapper_struct RHS_Wrapper_params;

    Spacetime_Base_Class* p_Spacetime;
    Step_controller_class* p_Step_controller;

    bool continue_integration;
    bool integration_complete;

    const double* const get_current_State_Vector() const;
    const double* const get_previous_State_Vector() const;

    int get_implicit_method_system(const gsl_vector* State_Vector, void* Params, gsl_vector* System_to_solve);

    void Propagate_ray();

};