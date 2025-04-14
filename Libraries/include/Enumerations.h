#pragma once

#include <type_traits>

enum Spacetime_enums: std::underlying_type_t <std::byte> {

    Kerr			      = 0,
    Wormhole		      = 1,
    Reg_Black_Hole        = 2,
    Janis_Newman_Winicour = 3,
    Einstein_Gauss_Bonnet = 4,
    BH_w_Dark_Matter      = 5,
    SPACETIME_NUMBER      = 6

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

};

enum Profile_enums: std::underlying_type_t <std::byte> {

    e_Power_law_profile       = 0,
    e_Exponential_law_profile = 1,
    e_Gaussian_profile        = 2,
    e_Spherical_profile       = 3,

};

enum Velocity_enums: std::underlying_type_t <std::byte> {

    e_Keplarian = 0,
    e_Theta_dependant = 1

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

    e_step    = 8,

    e_State_Number = 9,

};

enum XYZ_enums: std::underlying_type_t <std::byte> {

    x = 0,
    y = 1,
    z = 2

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

};

enum Metric_Parameter_Selector: std::underlying_type_t <std::byte> {

    Spin = 0,
    WH_Redshift = 1,
    JNW_Gamma = 2,
    GB_Gamma = 3,
    RBH_Param = 4,
    BH_w_DM_Halo_Compactness = 5,
    BH_w_DM_Halo_M_Halo = 6

};

