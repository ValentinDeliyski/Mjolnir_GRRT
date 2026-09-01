#pragma once

#include "Enumerations.h"
#include "Constants.h"
#include "Integrators.h"
#include "Structs.h"

#include "gsl/gsl_interp2d.h"
#include "gsl/gsl_spline.h"
#include "gsl/gsl_spline2d.h"

class Emission_Integrator_class {

private:

    static constexpr int RK54_size = 7;
    static constexpr int RK78_size = 13;

    /* ====================================================== The adaptive DP RK5(4) Butcher table ====================================================== */

    const double RK54_Coeff_deriv[RK78_size][RK78_size - 1] =
    {
        {       0,               0,              0,             0,            0,           0   , 0., 0., 0., 0., 0., 0.},
        {    1. / 5,             0,              0,             0,            0,           0   , 0., 0., 0., 0., 0., 0.},
        {    3. / 40,         9. / 40,           0,             0,            0,           0   , 0., 0., 0., 0., 0., 0.},
        {   44. / 45,       -56. / 15,       32. / 9,           0,            0,           0   , 0., 0., 0., 0., 0., 0.},
        {19372. / 6561,  -25360. / 2187,  64448. / 6561,  -212. / 729,        0,           0   , 0., 0., 0., 0., 0., 0.},
        { 9017. / 3168,    -355. / 33,    46732. / 5247,    49. / 176, -5103. / 18656,     0   , 0., 0., 0., 0., 0., 0.},
        {   35. / 384,           0,         500. / 1113,   125. / 192, -2187. / 6784,  11. / 84, 0., 0., 0., 0., 0., 0.},
        {       0.,              0,              0.,            0.,           0.,          0.,   0., 0., 0., 0., 0., 0.},
        {       0.,              0,              0.,            0.,           0.,          0.,   0., 0., 0., 0., 0., 0.},
        {       0.,              0,              0.,            0.,           0.,          0.,   0., 0., 0., 0., 0., 0.},
        {       0.,              0,              0.,            0.,           0.,          0.,   0., 0., 0., 0., 0., 0.},
        {       0.,              0,              0.,            0.,           0.,          0.,   0., 0., 0., 0., 0., 0.},
        {       0.,              0,              0.,            0.,           0.,          0.,   0., 0., 0., 0., 0., 0.}
    };

    const double RK54_Coeff_sol_main[RK78_size] = { 35. / 384, 0, 500. / 1113, 125. / 192, -2187. / 6784, 11. / 84, 0., 0., 0.,0., 0., 0., 0. };
    const double RK54_Coeff_test_embeded[RK78_size] = { 5179. / 57600, 0, 7571. / 16695, 393. / 640, -92097. / 339200, 187. / 2100, 1. / 40, 0., 0., 0., 0., 0., 0. };
    const double RK54_Coeff_affine_param[RK78_size] = { 0.0, 1. / 5, 3. / 10, 4. / 5, 8. / 9, 1.0, 1.0, 0., 0., 0., 0., 0., 0. };

    /* ====================================================== The adaptive RK-Fehlberg 7(8) Butcher table ====================================================== */

    // The coefficients are from https://ntrs.nasa.gov/api/citations/19680027281/downloads/19680027281.pdf, table X
    const double RK78_Fhelberg_Coeff_deriv[RK78_size][RK78_size - 1] =
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

    const double RK78_Fhelberg_Coeff_sol_main[RK78_size] = { 41. / 840, 0, 0, 0, 0, 34. / 105, 9. / 35, 9. / 35, 9. / 280, 9. / 280, 41. / 840, 0., 0. };
    const double RK78_Fhelberg_Coeff_sol_embeded[RK78_size] = { 0, 0, 0, 0, 0, 34. / 105, 9. / 35, 9. / 35, 9. / 280, 9. / 280, 0, 41. / 840, 41. / 840 };
    const double RK78_Fhelberg_Coeff_affine_param[RK78_size] = { 0.0, 2. / 27, 1. / 9, 1. / 6, 5. / 12, 1. / 2, 5. / 6, 1. / 6, 2. / 3, 1. / 3, 1., 0., 1. };

    /* ====================================================== The adaptive DP RK8(7) Butcher table ====================================================== */

    // The coefficients are from https://www.sciencedirect.com/science/article/pii/0771050X81900103?via%3Dihub, table 2
    const double RK78_DP_Coeff_deriv[RK78_size][RK78_size - 1] =
    {
        {              0,              0,         0,                0,                         0,                          0,                          0,                         0,                          0,                        0,                        0,            0 },
        {           1. / 18,           0,         0,                0,                         0,                          0,                          0,                         0,                          0,                        0,                        0,            0 },
        {           1. / 48,        1. / 16,      0,                0,                         0,                          0,                          0,                         0,                          0,                        0,                        0,            0 },
        {           1. / 32,           0,      3. / 32,             0,                         0,                          0,                          0,                         0,                          0,                        0,                        0,            0 },
        {           5. / 16,           0,    -75. / 64,         75. / 64,                      0,                          0,                          0,                         0,                          0,                        0,                        0,            0 },
        {           3. / 80,           0,         0,             3. / 16,                   3. / 20,                       0,                          0,                         0,                          0,                        0,                        0,            0 },
        {    29443841. / 614563906,    0,         0,      77736538. / 692538347,    -28693883. / 1125000000,     23124283. / 1800000000,               0,                         0,                          0,                        0,                        0,            0 },
        {    16016141. / 946692911,    0,         0,      61564180. / 158732637,     22789713. / 633445777,     545815736. / 2771057229,   -180193667. / 1043307555,              0,                          0,                        0,                        0,            0 },
        {    39632708. / 573591083,    0,         0,    -433636366. / 683701615,   -421739975. / 2616292301,    100302831. / 723423059,     790204164. / 839813087,    800635310. / 3783071287,               0,                        0,                        0,            0 },
        {   246121993. / 1340847787,   0,         0,  -37695042795. / 15268766246, -309121744. / 1061227803,    -12992083. / 490766935,    6005943493. / 2108947869,   393006217. / 1396673457,    123872331. / 1001029789,             0,                        0,            0 },
        { -1028468189. / 846180014,    0,         0,    8478235783. / 508512852,   1311729495. / 1432422823, -10304129995. / 1701304382, -48777925059. / 3047939560, 15336726248. / 1032824649, -45442868181. / 3398467696, 3065993473. / 597172653,              0,            0 },
        {   185892177. / 718116043,    0,         0,   -3185094517. / 667107341,   -477755414. / 1098053517,   -703635378. / 230739211,    5731566787. / 1027545527,  5232866602. / 850066563,   -4093664535. / 808688257,  3962137247. / 1805957418,   65686358. / 487910083,  0 },
        {   403863854. / 491063109,    0,         0,   -5068492393. / 434740067,   -411421997. / 543043805,     652783627. / 914296604,   11173962825. / 925320556, -13158990841. / 6184727034,   3936647629. / 1978049680, -160528059. / 685178525,   248638103. / 1413531060, 0 }
    };

    const double RK78_DP_Coeff_sol_main[RK78_size] = { 13451932. / 455176623, 0, 0, 0, 0, -808719846. / 976000145, 1757004468. / 5645159321, 656045339. / 265891186, -3867574721. / 1518517206, 465885868. / 322736535,   53011238. / 667516719,          2. / 45,            0 };
    const double RK78_DP_Coeff_sol_embeded[RK78_size] = { 14005451. / 335480064, 0, 0, 0, 0,  -59238493. / 1068277825, 181606767. / 758867731,  561292985. / 797845732, -1041891430. / 1371343529, 760417239. / 1151165299, 118820643. / 751138087, -528747749. / 2220607170, 1. / 4 };
    const double RK78_DP_Coeff_affine_param[RK78_size] = { 0.0, 1. / 18, 1. / 12, 1. / 8, 5. / 16, 3. / 8, 59. / 400, 93. / 200, 5490023248. / 9719169821, 13. / 20, 1201146811. / 1299019798, 1.0, 1.0 };

    Integrator_enums e_Active_rad_transfer_integrator;
    Integrator_enums e_Active_parallel_transport_integrator;

    /* ---------------------- Counters ---------------------- */

    size_t Current_emission_log_idx;
    size_t Current_polarization_log_idx;

    /* ------ Holds a shared pointer to the main results log struct ------ */
    std::shared_ptr<Ray_log_type> p_Ray_log_struct;

    /* ------ Holds a shared pointer to the polarization debug log struct ------ */
    std::shared_ptr<Polarization_debug_type> p_Polarization_debug_struct;

    /* ------ Holds a shared pointer to the sim context ----- */
    const Simulation_Context_type* p_Sim_Context;


    /* ---------------------------- Pointers to the geodesic spline instance ----------------------------- */

    // The +2 is because I need a spline of the geodesic in both local and global coordinates (so far only the radial coordinate and momentum differs)
    gsl_spline* p_Ray_spline_instance[e_Dynamic_state_size + 2];
    gsl_interp_accel* p_Spline_accelerator;

    /* ------------- These are convenient storage for dynamical states and their derivatives ------------- */

    std::complex<double> Parallel_Transport_RHS_log[RK78_size * e_Stokes_param_num]{};
    double Rad_Transfer_RHS_log[RK78_size * e_Stokes_param_num]{};

    std::complex<double> Current_Pol_Vector[e_Stokes_param_num]{};
    double Current_Stokes_Vector[e_Stokes_param_num]{};
    double Current_State_Vector[e_Dynamic_state_size]{};
    double Current_Affine_Param{};

    std::complex<double> Temp_Pol_Vector[e_Stokes_param_num]{};
    double Temp_Stokes_Vector[e_Stokes_param_num]{};


    /* --------------------------------- Internal functions --------------------------------- */

    void Update_emission_log();
    void Update_polarization_log();

    const double* const get_ray_Local_State_Vector(const double Affine_param);

    const double* const get_ray_Global_State_Vector(const double Affine_param);

    void get_Radiative_transfer_RHS(const double* const Emission_Functions,
                                    const double* const Absorbtion_Functions,
                                    const double* const Faraday_Functions,
                                    const double* const Stokes_Vector,
                                    double* const RHS);

    void get_Parallel_Transport_RHS(const double* const Global_State_Vector,
                                    std::complex<double>* Vector_to_transport,
                                    Tensor_type_enums e_Vec_type,
                                    std::complex<double>* const RHS);

    void Get_radiative_transfer_operators(const double* const Absorbtion_functions,
                                          const double* const Faraday_functions,
                                          double const CGS_Step,
                                          double Transfer_Operator[e_Stokes_param_num][e_Stokes_param_num],
                                          double Integrated_Transfer_Operator[e_Stokes_param_num][e_Stokes_param_num]);

    void Run_Analytic_Stokes_Vector_Propagator(const double Start_Affine_Param, const double End_Affine_Param);

    void __Run_Analytic_No_Faraday_conversion_propagator(const double CGS_Step, const Transfer_functions_type* p_Transfer_Functions);

    void __Run_Analytic_Pure_Polarized_Emission_propagator(const double CGS_Step, const Transfer_functions_type* p_Transfer_Functions);

    void __Run_Analytic_No_Absorbtion_propagator(const double CGS_Step, const Transfer_functions_type* p_Transfer_Functions);

    void __Run_Analytic_Unpolarized_Emission_propagator(const double CGS_Step, const Transfer_functions_type* p_Transfer_Functions);

    void Run_Runge_Kutta_Stokes_Vector(const double Start_Affine_Param, const double End_Affine_Param);

public:

    Emission_Integrator_class(const Simulation_Context_type* p_Sim_Context, Results_type &p_Ray_results);
    ~Emission_Integrator_class();

    template<typename Vec_type>
    void set_Polarization_Vector(const Vec_type* const Polarization_Vector) {

        /* NOTE: Do not replace this with a memcpy call! */

        for (int idx = 0; idx < e_Stokes_param_num; idx++) {

            this->Current_Pol_Vector[idx] = Polarization_Vector[idx];

        }

    };

    void Map_Stokes_to_Polarization_Vector(const double Stokes_Tetrad[e_Stokes_param_num][e_Stokes_param_num], bool Map_Between_Intermediate);

    void Map_Polarization_Vector_to_Stokes(const double inv_Stokes_Tetrad[e_Stokes_param_num][e_Stokes_param_num], bool Map_Between_Intermediate);

    void normalize_polarization_vector(const double Affine_param);

    void Propagate_Stokes_Vector(const double Start_Affine_Param, const double End_Affine_Param);

    void Propagate_Polarization_Vector(const double Start_Affine_Param, const double End_Affine_Param, const Tensor_type_enums Vec_type);

    const double* const get_current_Stokes_Vector() const;

    const std::complex<double>* const get_current_Polarization_Vector() const;

};