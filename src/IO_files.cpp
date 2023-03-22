#pragma once

#define _USE_MATH_DEFINES

#include <iostream>

#include "IO_files.h"
#include "Constants.h"

extern const bool truncate;
extern double hotspot_phi_angle;

std::string File_Names[] = {

    "Kerr_n=0",
    "Kerr_n=1",
    "Kerr_n=2",
    "Kerr_n=3",

    "Kerr_Momenta_n=0",
    "Kerr_Momenta_n=1",
    "Kerr_Momenta_n=2",
    "Kerr_Momenta_n=3",

    "Wormhole_n=0",
    "Wormhole_n=1",
    "Wormhole_n=2",
    "Wormhole_n=3",

    "Wormhole_Momenta_n=0",
    "Wormhole_Momenta_n=1",
    "Wormhole_Momenta_n=2",
    "Wormhole_Momenta_n=3",

    "RBH_n=0",
    "RBH_n=1",
    "RBH_n=2",
    "RBH_n=3",

    "RBH_Momenta_n=0",
    "RBH_Momenta_n=1",
    "RBH_Momenta_n=2",
    "RBH_Momenta_n=3",

    "JNW_n=0",
    "JNW_n=1",
    "JNW_n=2",
    "JNW_n=3",

    "JNW_Momenta_n=0",
    "JNW_Momenta_n=1",
    "JNW_Momenta_n=2",
    "JNW_Momenta_n=3",

    "Gauss_Bonnet_n=0",
    "Gauss_Bonnet_n=1",
    "Gauss_Bonnet_n=2",
    "Gauss_Bonnet_n=3",

    "Gauss_Bonnet_Momenta_n=0",
    "Gauss_Bonnet_Momenta_n=1",
    "Gauss_Bonnet_Momenta_n=2",
    "Gauss_Bonnet_Momenta_n=3",

    "BH_Dark_Matter_Halo_n=0",
    "BH_Dark_Matter_Halo_n=1",
    "BH_Dark_Matter_Halo_n=2",
    "BH_Dark_Matter_Halo_n=3",

    "BH_Dark_Matter_Halo_Momenta_n=0",
    "BH_Dark_Matter_Halo_Momenta_n=1",
    "BH_Dark_Matter_Halo_Momenta_n=2",
    "BH_Dark_Matter_Halo_Momenta_n=3"

};

int open_output_files(std::ofstream data[], std::ofstream momentum_data[]) {

    const int File_number = SPACETIME_NUMBER * 8;

    std::filesystem::path dir(std::filesystem::current_path() / "Sim_Results"); // Main results directory
    std::filesystem::create_directories(dir / "Hotspot_animation_frames");      // Only relevant for simulation mode 3

    std::filesystem::path file[File_number];
    std::filesystem::path full_path[File_number];
    std::filesystem::path file_extention(".txt");

    for (int File_Index = 0; File_Index <= File_number - 1; File_Index += 1) {

        if (Active_Sim_Mode == 3) {

            file[File_Index] = File_Names[File_Index] + "_" + std::to_string(hotspot_phi_angle * 10);

            file[File_Index].replace_extension(file_extention);

            full_path[File_Index] = dir / "Hotspot_animation_frames" / file[File_Index];

        }
        else {

            file[File_Index] = File_Names[File_Index];

            file[File_Index].replace_extension(file_extention);

            full_path[File_Index] = dir / file[File_Index];
        }

    }

    auto open_type = std::ios::app;

    if (truncate) {

        open_type = std::ios::trunc;

    }

    for (int File_Index = 0; File_Index <= ORDER_NUM - 1; File_Index += 1) {

        data[File_Index].open(full_path[File_Index + 2 * e_metric * 4], open_type);

        if (Active_Sim_Mode == 2) {

            momentum_data[File_Index].open(full_path[File_Index + (2 * e_metric + 1) * 4], open_type);

        }

    }

    for (int Image_order = direct; Image_order <= ORDER_NUM - 1; Image_order += 1) {

        data[Image_order] << "Observer Radial Position [M]: " << r_obs << '\n'
                          << "Observer Inclination [Deg]: " << theta_obs * 180 / M_PI
                          << '\n'
                          << "Observation Window Dimentions (-X,+X,-Y,+Y) [Rad]: "
                          << H_angle_min << ","
                          << H_angle_max << ","
                          << V_angle_min << ","
                          << H_angle_max
                          << '\n'
                          << "Simulation Resolutoin: "
                          << RESOLUTION
                          << " x "
                          << RESOLUTION
                          << '\n';

        switch (e_metric) {

            case Kerr:

                data[Image_order] << "Spin Parameter: " << SPIN
                                  << '\n';

                break;

            case Wormhole:

                data[Image_order] << "Spin Parameter: " << SPIN << ", Redshift Parameter: " << WH_REDSHIFT
                                  << '\n';
                break;

            case Reg_Black_Hole:

                data[Image_order] << "Parameter: " << RBH_PARAM
                                  << '\n';
                break;

            case Naked_Singularity:

                data[Image_order] << "Gamma: " << JNW_GAMMA
                                  << '\n';
                break;

            case Gauss_Bonnet:

                data[Image_order] << "Gamma: " << GAUSS_BONNET_GAMMA
                                  << '\n';
            break;

            case BH_w_Dark_Matter:

                data[Image_order] << "Mass_Halo: " << M_HALO << ", Halo Length Scale: " << A_0
                                  << '\n';
           break;

        }

        data[Image_order] << "Image X Coord [M],"
                          << " "
                          << "Image Y Coord [M],"
                          << " "
                          << "Novikov-Thorne Disk Redshift [-],"
                          << " "
                          << "Novikov-Thorne Flux [M^-2],"
                          << " "
                          << "Optically Thin Disk Intensity [Jy/sRad],"
                          << " "
                          << "Novikov-Thorne Source r Coord [M],"
                          << " "
                          << "Novikov-Thorne Source phi Coord [M] [Rad]"
                          << '\n';


        if (Active_Sim_Mode == 2) {

            momentum_data[Image_order] << "Radial Momentum (covariant),"
                << " "
                << "Theta Momentum (covariant),"
                << " "
                << "Phi Momentum (covariant),"
                << " "
                << "Metric Parameter"
                << '\n';

        }
    }

    return OK;

}

int close_output_files(std::ofstream data[], std::ofstream momentum_data[]) {

    for (int i = 0; i <= 3; i++) {

        data[i].close();

        if (Active_Sim_Mode == 2) {

            momentum_data[i].close();
        }

    }

    return OK;

}

int get_geodesic_data(double J_data[], double p_theta_data[], int* Data_number) {

    std::ifstream geodesic_data;

    geodesic_data.open(input_file_path, std::ios::in);

    while (true) {

        geodesic_data >> J_data[*Data_number] >> p_theta_data[*Data_number];

        if (geodesic_data.eof()) {

            break;

        }

        *Data_number += 1;

    }

    geodesic_data.close();

    return OK;
}

int write_to_file(Results_type s_Ray_results, std::ofstream data[], std::ofstream momentum_data[]) {

    for (int Image_order = direct; Image_order <= ORDER_NUM - 1; Image_order += 1) {

        data[Image_order] << s_Ray_results.Image_Coords[x]
                          << " "
                          << s_Ray_results.Image_Coords[y]
                          << " "
                          << s_Ray_results.Redshift_NT[Image_order]
                          << " "
                          << s_Ray_results.Flux_NT[Image_order]
                          << " "
                          << s_Ray_results.Intensity[Image_order] * CGS_TO_JANSKY
                          << " "
                          << s_Ray_results.Source_Coords[e_r][Image_order]
                          << " "
                          << s_Ray_results.Source_Coords[e_phi][Image_order]
                          << '\n';

        if (Active_Sim_Mode == 2) {

            momentum_data[Image_order] << s_Ray_results.Three_Momentum[e_r][Image_order]
                                       << " "
                                       << s_Ray_results.Three_Momentum[e_theta][Image_order]
                                       << " "
                                       << s_Ray_results.Three_Momentum[e_phi][Image_order]
                                       << " "
                                       << s_Ray_results.Parameters[e_metric]
                                       << '\n';

        }
    }

    return OK;

}

