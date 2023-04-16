#pragma once

#define _USE_MATH_DEFINES

#include <iostream>

#include "IO_files.h"
#include "Constants.h"

extern Spacetime_Base_Class* Spacetimes[];
extern Optically_Thin_Toroidal_Model OTT_Model;

File_manager_class::File_manager_class(std::string input_file_path, bool truncate) {

    Input_file_path_sim_mode_2 = input_file_path;
    Truncate_files = truncate;

}

void File_manager_class::get_geodesic_data(double J_data[], double p_theta_data[], int* Data_number) {

    std::ifstream geodesic_data;

    geodesic_data.open(Input_file_path_sim_mode_2, std::ios::in);

    while (true) {

        geodesic_data >> J_data[*Data_number] >> p_theta_data[*Data_number];

        if (geodesic_data.eof()) {

            break;

        }

        *Data_number += 1;

    }

    geodesic_data.close();

}

void File_manager_class::write_simulation_metadata() {

    for (int Image_order = direct; Image_order <= ORDER_NUM - 1; Image_order += 1) {

        Image_Output_files[Image_order] << "============================================================ SIMULATION METADATA ============================================================"
                                        << "\n"
                                        << "Spacetime: "
                                        << Image_File_Names[e_metric]
                                        << "\n"
                                        << "Image Order [-]: "
                                        << Image_order
                                        << "\n"
                                        << "Observer Radial Position [M]: " << r_obs << '\n'
                                        << "Observer Inclination [Deg]: " << theta_obs * 180 / M_PI
                                        << '\n'
                                        << "------------------------------------------------------- Optically Thin Disk Metadata -------------------------------------------------------"
                                        << "\n"
                                        << "Active disk profile: ";

        switch (e_disk_model) {

        case Power_law:

            Image_Output_files[Image_order] << "Power law (rho ~ 1 / r^2)"
                                            << "\n"
                                            << "Disk Opening Angle [tan(angle)]: "
                                            << DISK_OPENING_ANGLE
                                            << "\n"
                                            << "R_0[M]: "
                                            << R_0
                                            << "\n";

            break;

        case Exponential_law:

            Image_Output_files[Image_order] << "Exponential law (rho ~ exp(-r^2))"
                                            << "\n"
                                            << "Disk Height Scale [M]: "
                                            << DISK_HEIGHT_SCALE
                                            << "\n"
                                            << "Disk Radial Scale [M]: "
                                            << DISK_RADIAL_SCALE
                                            << "\n";

            break;

        }

        switch (e_emission) {

        case Synchotron_phenomenological:

            Image_Output_files[Image_order] << "Emission Model: " 
                                            << "Phenomenological"
                                            << "\n"
                                            << "Emission Power Law Exponent [-]: "
                                            << EMISSION_POWER_LAW
                                            << "\n"
                                            << "Absorbtion Coefficient [-]: "
                                            << DISK_ABSORBTION_COEFF
                                            << "\n"
                                            << "Source Function Power Law Exponent [-]: "
                                            << SOURCE_F_POWER_LAW
                                            << "\n"
                                            << "Emission Scale [Jy / sRad]: "
                                            << EMISSION_SCALE_PHENOMENOLOGICAL
                                            << "\n";

            break;

        case Synchotron_exact:

            Image_Output_files[Image_order] << "Emission Model: "
                                            << "Exact"
                                            << "\n"
                                            << "Disk Magnetization [-]: "
                                            << DISK_MAGNETIZATION
                                            << "\n"
                                            << "Magnetic Field Geometry: "
                                            << "[" << MAG_FIELD_GEOMETRY[0]
                                            << ", " << MAG_FIELD_GEOMETRY[1]
                                            << ", " << MAG_FIELD_GEOMETRY[2] << "]"
                                            << "\n"
                                            << "Max Electron Temperature [K]: "
                                            << T_ELECTRON_EXACT_CGS
                                            << "\n"
                                            << "Max Disk Density [g / cm^3]: "
                                            << N_ELECTRON_EXACT_CGS
                                            << "\n";

            break;

        }

        Image_Output_files[Image_order] << "------------------------------------------------------------- Hotpost Metadata -------------------------------------------------------------"
                                        << "\n"
                                        << "Hotspot Relative Scale [-]: "
                                        << HOTSPOT_REL_SCALE
                                        << "\n"
                                        << "Hotspot Radial Scale [M]: "
                                        << HOTSPOT_SCALE
                                        << "\n"
                                        << "Hotspot r_center [M]: "
                                        << HOTSPOT_R_COORD
                                        << "\n";


        Image_Output_files[Image_order] << "------------------------------------------------------- Novikov-Thorne Disk Metadata -------------------------------------------------------"
                                        << "\n"
                                        << "Inner Disk Radius [M]: "
                                        << r_in * (r_in != NULL) + Spacetimes[e_metric]->get_ISCO()[Inner] * (r_in == NULL)
                                        << "\n"
                                        << "Outer Disk Radius [M]: "
                                        << r_out
                                        << "\n";


        switch (Active_Sim_Mode) {

        case 2:

            Image_Output_files[Image_order] << "Image X Coord [M],"
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
                                            << "Novikov-Thorne Source phi Coord [Rad],"
                                            << " "
                                            << "Radial Momentum (covariant),"
                                            << " "
                                            << "Theta Momentum (covariant),"
                                            << " "
                                            << "Phi Momentum (covariant),"
                                            << " "
                                            << "Metric Parameter"
                                            << '\n';

            break;

        default:

            Image_Output_files[Image_order] << "Observation Window Dimentions (-X,+X,-Y,+Y) [Rad]: "
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

                Image_Output_files[Image_order] << "Spin Parameter: " << SPIN
                    << '\n';
                break;

            case Wormhole:

                Image_Output_files[Image_order] << "Spin Parameter: " << SPIN << ", Redshift Parameter: " << WH_REDSHIFT
                    << '\n';
                break;

            case Reg_Black_Hole:

                Image_Output_files[Image_order] << "Parameter: " << RBH_PARAM
                    << '\n';
                break;

            case Naked_Singularity:

                Image_Output_files[Image_order] << "Gamma: " << JNW_GAMMA
                    << '\n';
                break;

            case Gauss_Bonnet:

                Image_Output_files[Image_order] << "Gamma: " << GAUSS_BONNET_GAMMA
                    << '\n';
                break;

            case BH_w_Dark_Matter:

                Image_Output_files[Image_order] << "Mass_Halo: " << M_HALO << ", Halo Length Scale: " << A_0
                    << '\n';
                break;

            }

            Image_Output_files[Image_order] << "Image X Coord [M],"
                                            << " "
                                            << "Image Y Coord [M],"
                                            << " "
                                            << "Novikov-Thorne Disk Redshift [-],"
                                            << " "
                                            << "Novikov-Thorne Flux [M^-2],"
                                            << " "
                                            << "Optically Thin Disk Intensity [Jy/sRad]"
                                            << '\n';

        }

    }
}


void File_manager_class::open_image_output_files() {

    const int File_number = SPACETIME_NUMBER;

    std::filesystem::path dir(std::filesystem::current_path() / "Sim_Results"); // Main results directory
    std::filesystem::create_directories(dir / "Hotspot_animation_frames");      // Only relevant for simulation mode 3

    std::filesystem::path Image_file_names[File_number];
    std::filesystem::path Momenta_file_names[File_number];

    std::filesystem::path Image_full_path[ORDER_NUM];
    std::filesystem::path Momenta_full_path[ORDER_NUM];
    std::filesystem::path file_extention(".txt");

    auto open_type = std::ios::app;

    if (Truncate_files) {

        open_type = std::ios::trunc;

    }

    for (int File_Index = 0; File_Index <= ORDER_NUM - 1; File_Index += 1) {

        if (Active_Sim_Mode == 3) {

            Image_file_names[File_Index] = Image_File_Names[e_metric]
                + "_n"
                + std::to_string(File_Index)
                + "_"
                + std::to_string(OTT_Model.HOTSPOT_PHI_COORD * 10);

            Image_file_names[File_Index].replace_extension(file_extention);
            Image_full_path[File_Index] = dir / "Hotspot_animation_frames" / Image_file_names[File_Index];

        }
        else {

            Image_file_names[File_Index] = Image_File_Names[e_metric]
                + "_n"
                + std::to_string(File_Index);

            Image_file_names[File_Index].replace_extension(file_extention);
            Image_full_path[File_Index] = dir / Image_file_names[File_Index];
        }

        Image_Output_files[File_Index].open(Image_full_path[File_Index], open_type);

    }

    File_manager_class::write_simulation_metadata();

}

void File_manager_class::open_log_output_file() {

    std::filesystem::path dir(std::filesystem::current_path() / "Sim_Results"); // Main results directory

    std::filesystem::path Log_File_Name;
    std::filesystem::path Log_File_full_path;
    std::filesystem::path file_extention(".txt");

    auto open_type = std::ios::trunc;

    Log_File_Name = "Photon_Log";
    Log_File_Name.replace_extension(file_extention);
    Log_File_full_path = dir / Log_File_Name;
     
    Log_Output_File.open(Log_File_full_path, open_type);

    File_manager_class::write_simulation_metadata();

}

void File_manager_class::write_image_data_to_file(Results_type s_Ray_results) {

    for (int Image_order = direct; Image_order <= ORDER_NUM - 1; Image_order += 1) {

        Image_Output_files[Image_order] << s_Ray_results.Image_Coords[x]
                                  << " "
                                  << s_Ray_results.Image_Coords[y]
                                  << " "
                                  << s_Ray_results.Redshift_NT[Image_order]
                                  << " "
                                  << s_Ray_results.Flux_NT[Image_order]
                                  << " "
                                  << s_Ray_results.Intensity[Image_order] * CGS_TO_JANSKY;

        if (Active_Sim_Mode == 2) {

            Image_Output_files[Image_order] << " "
                                            << s_Ray_results.Source_Coords[e_r][Image_order]
                                            << " "
                                            << s_Ray_results.Source_Coords[e_phi][Image_order]
                                            << " "
                                            << s_Ray_results.Three_Momentum[e_r][Image_order]
                                            << " "
                                            << s_Ray_results.Three_Momentum[e_theta][Image_order]
                                            << " "
                                            << s_Ray_results.Three_Momentum[e_phi][Image_order]
                                            << " "
                                            << s_Ray_results.Parameters[e_metric];

        }

        Image_Output_files[Image_order] << '\n';
    }

};

void File_manager_class::log_photon_path(Results_type s_Ray_results) {

    for (int log_index = 0; log_index <= s_Ray_results.Ray_log.size() - 1; log_index++) {

        for (int state_index = e_r; state_index <= e_State_Number - 1; state_index++) {

            Log_Output_File << s_Ray_results.Ray_log[log_index][state_index] << " ";
          
        }

        Log_Output_File << '\n';
    }
   
};

void File_manager_class::close_image_output_files() {

    for (int File_Index = 0; File_Index <= ORDER_NUM - 1; File_Index++) {

        Image_Output_files[File_Index].close();

    }

}

void File_manager_class::close_log_output_file() {

    Log_Output_File.close();

}