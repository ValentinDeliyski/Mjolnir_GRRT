#pragma once

#ifndef ENUMS

    #define ENUMS

    #include <vector>

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

        e_r             = 0,
        e_theta         = 1,
        e_phi           = 2,
        e_p_phi         = 3,
        e_p_theta       = 4,
        e_p_r           = 5,
        e_Intensity     = 6,
        e_Optical_Depth = 7,
        e_State_Number  = 8

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

    struct Initial_conditions_type {

        double init_metric[4][4];
        double init_metric_Redshift_func;
        double init_metric_Shitft_func;

        double init_Pos[3];
        double init_Three_Momentum[3];

    };

    struct Results_type {

        double Flux_NT[ORDER_NUM]{};
        double Redshift_NT[ORDER_NUM]{};

        double Intensity[ORDER_NUM]{};
        double Optical_Depth{};

        double Source_Coords[3][ORDER_NUM]{};
        double Three_Momentum[3][ORDER_NUM]{};

        double Image_Coords[2]{};

        std::vector<std::vector<double>> Ray_log;

        double Parameters[SPACETIME_NUMBER]{};

    };

    struct Trace_Results_type {

        double Intensity;
        double Optical_Depth;

        double Coordinates[3];
        double Three_Momentum[3];

        double Image_Coords[2];

    };

#endif
