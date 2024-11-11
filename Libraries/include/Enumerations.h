#pragma once

#ifndef ENUMS

    #define ENUMS

    enum Spacetime_enums {

        Kerr			      = 0,
        Wormhole		      = 1,
        Reg_Black_Hole        = 2,
        Janis_Newman_Winicour = 3,
        Einstein_Gauss_Bonnet = 4,
        BH_w_Dark_Matter      = 5,
        SPACETIME_NUMBER      = 6

    };

    enum Ensamble_enums {

        e_Thermal_ensamble          = 0,
        e_Power_law_ensamble        = 1,
        e_Kappa_ensamble            = 2,
        e_Phenomenological_ensamble = 3,

    };

    enum Profile_enums {

        e_Power_law_profile       = 0,
        e_Exponential_law_profile = 1,
        e_Gaussian_profile        = 2

    };

    enum State_enums {

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

        OK    =  0,
        ERROR = -1

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

    enum Thermal_Syncotron_fit_selector {

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
        Implicit_Trapezoid = 1,
        RK4 = 2

    };

    enum Metric_Parameter_Selector {

        Spin = 0,
        WH_Redshift = 1,
        JNW_Gamma = 2,
        GB_Gamma = 3,
        RBH_Param = 4,
        BH_w_DM_Halo_Compactness = 5,
        BH_w_DM_Halo_M_Halo = 6

    };

#endif
