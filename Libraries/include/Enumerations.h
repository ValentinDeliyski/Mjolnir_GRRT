#pragma once

#include <type_traits>

enum Spacetime_enums: std::underlying_type_t <std::byte> {

    Kerr			      = 0,
    Wormhole		      = 1,
    Reg_Black_Hole        = 2,
    Janis_Newman_Winicour = 3,
    Einstein_Gauss_Bonnet = 4,
    BH_w_Dark_Matter      = 5,
    Numerical             = 6,
    Minkowski             = 7,
    e_Spacetime_number    = 8

};

enum Magnetic_field_geometry_enums: std::underlying_type_t <std::byte> {

    Toroidal = 0,
    Radial   = 1,
    Vertical = 2,
    Constant = 3

};

enum Magnetic_field_magnitude_enums: std::underlying_type_t <std::byte>{

    Power_law_based     = 0,
    Magnetization_based = 1,

};

enum Derivative_selector_enums : std::underlying_type_t <std::byte> {

    None = 0,
    First_radial_derivative = 1,
    Second_radial_derivative = 2,
    First_theta_derivative = 3,
    Second_mixed_derivative = 4,
    Second_theta_derivative = 5

};

enum Step_controller_type_enums: std::underlying_type_t <std::byte> {

    PID = 0,
    Gustafsson = 1

};

enum Tensor_type_enums: std::underlying_type_t <std::byte> {

    Covariant = 0,
    Contravariant = 1

};

enum Emission_medium_enums: std::underlying_type_t <std::byte> {

    Disk = 0,
    Hotspot = 1,
    Jet = 2

};

enum Ensamble_enums: std::underlying_type_t <std::byte> {

    e_Thermal_ensamble          = 0,
    e_Power_law_ensamble        = 1,
    e_Kappa_ensamble            = 2,
    e_Phenomenological_ensamble = 3,
    e_Debug_constant_functions  = 4

};

enum Profile_enums: std::underlying_type_t <std::byte> {

    e_Power_law             = 0,
    e_Gaussian              = 1,
    e_Hybrid_power_gaussian = 2,
    e_Spherical             = 3,

};

enum Disk_model_enums : std::underlying_type_t <std::byte> {

    /* This is the model used in https://arxiv.org/pdf/2206.12066, with an added cutoff exponential. */
    e_Phenom_RIAF_1 = 0,

    /* This is the model used in https://arxiv.org/pdf/2209.09931, with an added cutoff exponential. */
    e_Phenom_RIAF_2 = 1,

    /* This is the model used in https://iopscience.iop.org/article/10.3847/1538-4357/ab96c6/pdf */
    e_Colab_test_1 = 2,

    /* This is the model from https://articles.adsabs.harvard.edu/pdf/1974ApJ...191..499P */
    e_Page_Thorne = 3,

    /* This model exists for testing purposes (see the plasma integration tests in https://www.aanda.org/articles/aa/pdf/2020/09/aa38573-20.pdf )*/
    e_Debug_constant_density = 4

};

enum Numerical_Anzatz_enums : std::underlying_type_t <std::byte> {

    /* This is the anzatz used in Stoycho / Galin's paper - https://arxiv.org/pdf/2402.08469 */
    e_Anzatz_1 = 0,

    /* This is the anzatz used in the Heirdeiro paper - https://arxiv.org/pdf/1501.04319 */
    e_Anzatz_2 = 1

};

enum Velocity_enums: std::underlying_type_t <std::byte> {

    e_Keplarian = 0,
    e_Theta_dependant = 1,
    e_Circular_fixed_rate = 2

};

enum State_enums: std::underlying_type_t <std::byte> {

    e_t     = 0,
    e_r     = 1,
    e_theta = 2,
    e_phi   = 3,

    e_p_t     = 4,
    e_p_r     = 5,
    e_p_theta = 6,
    e_p_phi   = 7,

    e_Dynamic_state_size = 8,

    e_step         = 8,
    e_affine_param = 9,

    e_Full_state_size = 10,

};

enum XYZ_enums: std::underlying_type_t <std::byte> {

    e_x = 0,
    e_y = 1,
    e_z = 2

};

enum Image_Orders: std::underlying_type_t <std::byte> {

    e_direct       = 0,
    e_first        = 1,
    e_second       = 2,
    e_third        = 3,
    e_order_number = 4,

};

enum Return_Values: std::underlying_type_t <std::byte> {

    OK    =  0,
    ERROR =  1

};

enum Orbit_select: std::underlying_type_t <std::byte> {

    Inner = 0,
    Outer = 1,

};

enum Stokes_parameters: std::underlying_type_t <std::byte> { 

    I = 0,
    Q = 1,
    U = 2,
    V = 3,
    e_Stokes_param_num = 4

};

enum Radiative_Transfer_Integrator: std::underlying_type_t <std::byte> {

    Analytic = 0,
    Implicit_Trapezoid = 1,
    RK5 = 2,

};

enum Geodesic_Integrator_enums: std::underlying_type_t <std::byte> {

    RK5_fixed_step = 0,
    RK78_adaptive_step = 1,
    BDF_w_AB_predictor_order_1 = 2,
    BDF_w_AB_predictor_order_2 = 3,
    BDF_w_AB_predictor_order_3 = 4,
    BDF_w_AB_predictor_order_4 = 5,
    BDF_w_AB_predictor_order_5 = 6,

};