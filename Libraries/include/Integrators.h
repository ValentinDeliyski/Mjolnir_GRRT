#pragma once
#include "Enumerations.h"
#include "Spacetimes.h"
#include "Structs.h"
#include "Step_Controller.h"
#include "gsl\gsl_multiroots.h"

class Geodesic_Integrator_class {

private:

    static constexpr int RK54_size = 7;
    static constexpr int RK78_size = 13;
    static constexpr int ESDIRK54_size = 7;

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

    const double RK78_Fhelberg_Coeff_sol_main[RK78_size]    = { 41. / 840, 0, 0, 0, 0, 34. / 105, 9. / 35, 9. / 35, 9. / 280, 9. / 280, 41. / 840,     0,         0     };
    const double RK78_Fhelberg_Coeff_sol_embeded[RK78_size] = {     0,     0, 0, 0, 0, 34. / 105, 9. / 35, 9. / 35, 9. / 280, 9. / 280,     0,     41. / 840, 41. / 840 };

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

    const double RK78_DP_Coeff_sol_main[RK78_size]    = { 13451932. / 455176623, 0, 0, 0, 0, -808719846. / 976000145, 1757004468. / 5645159321, 656045339. / 265891186, -3867574721. / 1518517206, 465885868. / 322736535,   53011238. / 667516719,          2. / 45,            0 };
    const double RK78_DP_Coeff_sol_embeded[RK78_size] = { 14005451. / 335480064, 0, 0, 0, 0,  -59238493. / 1068277825, 181606767. / 758867731,  561292985. / 797845732, -1041891430. / 1371343529, 760417239. / 1151165299, 118820643. / 751138087, -528747749. / 2220607170, 1. / 4 };

    /* =================================================== The adaptive ESDIRK54 Butcher table ==================================================== */

    // The coefficients are from https://ntrs.nasa.gov/api/citations/20160005923/downloads/20160005923.pdf, table 25
    const double ESDIRK54_Coeff_deriv[ESDIRK54_size][ESDIRK54_size] =
    {
        {                 0.,                                0.,                               0.,                              0.,                              0.,                              0.,                  0.},
        {              23. / 125,                        23. / 125,                            0.,                              0.,                              0.,                              0.,                  0.},
        {   -121529886477. / 3189120653983,   -121529886477. / 3189120653983,              23. / 125,                           0.,                              0.,                              0.,                  0.},
        {    186345625210. / 8596203768457,    186345625210. / 8596203768457,   3681435451073. / 12579882114497,            23. / 125,                           0.,                              0.,                  0.},
        {  -9898129553915. / 11630542248213  -9898129553915. / 11630542248213, 19565727496993. / 11159348038501, 2073446517052. / 4961027473423,             23. / 125,                           0.,                  0.},
        { -39752543191591. / 7894275939720, -39752543191591. / 7894275939720,  52228808998390. / 5821762529307,  2756378382725. / 8748785577174, 17322065038796. / 10556643942083,            23. / 125,               0.},
        {  -1319096626979. / 17356965168099  -1319096626979. / 17356965168099,  4356877330928. / 10268933656267,  922991294344. / 3350617878647,  4729382008034. / 14755765856909, -308199069217. / 5897303561678, 23. / 125}
    };

    const double ESDIRK54_Coeff_sol_main[ESDIRK54_size]    = {  -1319096626979. / 17356965168099,    -1319096626979. / 17356965168099,   4356877330928. / 10268933656267,   922991294344. / 3350617878647,    4729382008034. / 14755765856909,   -308199069217. / 5897303561678,               23. / 125 };
    const double ESDIRK54_Coeff_sol_embeded[ESDIRK54_size] = { -12068858301481. / 111697653055985,  -12068858301481. / 111697653055985, 30204157393951. / 62440428688139, 26156819792768. / 110856972047457, 33531609809941. / 89326307438822, -18686091006953. / 578397443530870, 10582397456777. / 69011126173064 };

    bool Force_scatter;

    double Max_affine_param;
    double Max_integration_count;

    /* ---------------------- Debug flags ---------------------- */

    int N_steps_rejected;
    int NaN_checker_count;

    /* ------------- Integration termination flags ------------- */

    bool Normal_termination_condition;
    bool Max_affine_param_reached;
    bool Max_integration_count_reached;
    bool Step_too_small;

    /* --------------------------------------------------------- */

    double Intermediate_RHS_log[RK78_size * e_Dynamic_state_size]{};
    double Current_Dynamic_state[e_Dynamic_state_size];

    Integrator_enums e_Active_integrator;
    Initial_conditions_type* p_Init_conditions;

    Ray_log_type* p_Ray_log_struct;

    Spacetime_Base_Class* p_Spacetime;
    Step_controller_class* p_Step_controller;

    /* ----------- Root finder environment variables ----------- */

    gsl_vector* gsl_trial_State_Vector;
    gsl_multiroot_fsolver* Root_finder;
    gsl_multiroot_function Function_to_solve;

    /* --------------------------------- Internal functions --------------------------------- */

    Return_Values Run_NaN_checker(const double* const New_State, const double* const New_State_Embeded);

    void Run_Explicit_Runge_Kutta();

    void Run_ESDIRK54();

    void Check_integration_complete_status();

    void Update_ray_log(const double* const New_State_vector);
    void Update_debug_log();

public:

    Geodesic_Integrator_class(const Simulation_Context_type* const p_Sim_Context, Results_type* p_Ray_results);
    ~Geodesic_Integrator_class();

    RHS_wrapper_struct RHS_Wrapper_params;

    /* ------ This gets filled only in simulation mode 3 ----- */
    Adaptive_RK_Integrator_debug_type RK_Integrator_debug_log;

    bool continue_integration;
    bool integration_complete;

    bool Locate_event(Event_detection_enums e_Event, double* const State_at_Event_Global, double* const State_at_Event_Local);

    const double get_dense_output(const double Param, const State_enums idx, bool Is_current_RHS_evaluated) const;
    const double* const get_current_State_Vector_global() const;
    const double* const get_current_State_Vector_local() const;

    const double* const get_previous_State_Vector_global() const;
    const double* const get_previous_State_Vector_local() const;

    /* -------------------- This gets called by the RHS wrapper in the ESDIRK54 routine -------------------- */
    int get_implicit_method_system(const gsl_vector* State_Vector, void* Params, gsl_vector* System_to_solve);

    void Propagate_ray();

};