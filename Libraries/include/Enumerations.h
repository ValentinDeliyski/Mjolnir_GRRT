#pragma once

#ifndef ENUMS

    #define ENUMS

    enum Spacetime_enums {

        Kerr			  = 0,
        Wormhole		  = 1,
        Reg_Black_Hole    = 2,
        Naked_Singularity = 3,
        Gauss_Bonnet      = 4,
        BH_w_Dark_Matter  = 5,
        SPACETIME_NUMBER  = 6

    };

    enum Emission_model_enums {

        Synchotron_exact            = 0,
        Synchotron_phenomenological = 1

    };

    enum Disk_density_model_enums {

        Power_law       = 0,
        Exponential_law = 1

    };

    enum State_enums {

        e_r = 0,
        e_theta = 1,
        e_phi = 2,
        e_p_phi = 3,
        e_p_theta = 4,
        e_p_r = 5,
        e_Intensity = 6,
        e_Optical_Depth = 7,
        e_State_Number = 8

    };

    enum Ray_path_log_enums {

        e_path_log_r = 0,
        e_path_log_theta = 1,
        e_path_log_phi = 2,
        e_path_log_p_phi = 3,
        e_path_log_p_theta = 4,
        e_path_log_p_r = 5,
        e_path_log_step = 6,
        e_path_log_number = 7

    };

    enum Spacetime_coords {

        e_t_coord     = 0,
        e_r_coord     = 1,
        e_theta_coord = 2,
        e_phi_coord   = 3

    };

    enum XYZ_enums {

        x = 0,
        y = 1,
        z = 2

    };

    enum Image_Orders {

        direct    = 0,
        first     = 1,
        second    = 2,
        third     = 3,
        ORDER_NUM = 4,

    };

    enum Return_Values {

        OK    = 0,
        ERROR = 255

    };

    enum Disk_Model{

        Novikov_Thorne          = 0,
        Optically_Thin_Toroidal = 1,
        DISK_MODEL_NUM          = 2

    };

    enum Orbit_select{

        Inner = 0,
        Outer = 1,

    };

    enum Emission_Interpolation_points {

        Current  = 0,
        Next     = 1,
        INTERPOLATION_NUM = 2,

    };

    enum Sync_emission_fit_functions {

        Leung_2011  = 0,
        Dexter_2016 = 1

    };

    enum Stokes_parameters {

        I = 0,
        Q = 1,
        U = 2,
        V = 3,
        STOKES_PARAM_NUM = 4

    };

    enum Radiative_Transfer_Integrator {

        Analytic = 0,
        Implicit_Trapezoid = 1

    };

#endif
