#pragma once

#define _USE_MATH_DEFINES

#include "IO_files.h"
#include "Constants.h"
#include "Disk_models.h"
#include "Spacetimes.h"
#include <iostream>

File_manager_class::File_manager_class(Initial_conditions_type *p_Initial_Conditions, std::string input_file_path, bool truncate) {

    this->p_Initial_Conditions       = p_Initial_Conditions;
    this->Input_file_path_sim_mode_2 = input_file_path;
    this->Truncate_files             = truncate;
}

void File_manager_class::get_geodesic_data(double J_data[], double p_theta_data[], int* Data_number) {

    std::ifstream geodesic_data;
    std::string line;

    double J_input{};
    double P_input{};

    geodesic_data.open(Input_file_path_sim_mode_2, std::ios::in);

    while (true) {

        if (geodesic_data >> J_input >> P_input) {

            J_data[*Data_number] = J_input;
            p_theta_data[*Data_number] = P_input;
            
            *Data_number += 1;
        }

        if (geodesic_data.eof()) {

            break;
        }
    }

    geodesic_data.close();
}

void File_manager_class::write_simulation_metadata(int Sim_mode_2_number) {

    for (int Image_order = direct; Image_order <= ORDER_NUM - 1; Image_order += 1) {

        Image_Output_files[Image_order] << "============================================================ SIMULATION METADATA ============================================================"
                                        << "\n"
                                        << "Spacetime: "
                                        << Base_File_Names[e_metric]
                                        << "\n"
                                        << "Image Order [-]: "
                                        << Image_order
                                        << "\n"
                                        << "Observer Radial Position [M]: " << p_Initial_Conditions->init_Pos[e_r] << '\n'
                                        << "Observer Inclination [Deg]: " << theta_obs * 180 / M_PI << '\n'
                                        << "Observation Frequency [Hz]: " << OBS_FREQUENCY_CGS << '\n'
                                        << "Active Simulation Mode: " << Active_Sim_Mode
                                        << '\n';

        if (Active_Sim_Mode == 2) {


            Image_Output_files[Image_order] << "Number Of Photons Per Param Value: " << Sim_mode_2_number
                                            << '\n'
                                            << "Number Of Param Values: " << PARAM_SWEEP_NUMBER
                                            << '\n';
        }


        Image_Output_files[Image_order] << "------------------------------------------------------- Optically Thin Disk Metadata -------------------------------------------------------"
                                        << "\n"
                                        << "Active disk profile: ";

        switch (e_disk_model) {

        case Power_law:

            Image_Output_files[Image_order] << "Power law (rho ~ 1 / r^2)"
                                            << "\n"
                                            << "Disk Opening Angle [tan(angle)]: "
                                            << DISK_OPENING_ANGLE
                                            << "\n"
                                            << "R_0 [M]: "
                                            << DISK_R_0
                                            << "\n"
                                            << "R_Cutoff [M]: ";

            Image_Output_files[Image_order] << DISK_R_Cutoff << "\n";     
            Image_Output_files[Image_order] << "Cutoff Scale [M]: " << DISK_CUTOFF_SCALE << "\n";

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

        switch (e_disk_emission) {

        case Phenomenological_synchotron:

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

        case Thermal_dist_synchotron:

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
                                            << DISK_T_ELECTRON_SCALE_CGS
                                            << "\n"
                                            << "Max Disk Density [g / cm^3]: "
                                            << DISK_N_ELECTRON_SCALE_CGS
                                            << "\n";

            break;
        }

        Image_Output_files[Image_order] << "------------------------------------------------------------- Hotpost Metadata -------------------------------------------------------------"
                                        << "\n"
                                        << "Max Hotspot Density[g / cm ^ 3] "
                                        << HOTSPOT_N_ELECTRON_SLACE_CGS
                                        << "\n"
                                        << "Hotspot Radial Scale [M]: "
                                        << HOTSPOT_SPREAD
                                        << "\n"
                                        << "Hotspot r_center [M]: "
                                        << HOTSPOT_R_COORD
                                        << "\n";


        Image_Output_files[Image_order] << "------------------------------------------------------- Novikov-Thorne Disk Metadata -------------------------------------------------------"
                                        << "\n"
                                        << "Inner Disk Radius [M]: "
                                        << r_in
                                        << "\n"
                                        << "Outer Disk Radius [M]: "
                                        << r_out
                                        << "\n";


        if (Active_Sim_Mode != 2) {

            Image_Output_files[Image_order] << "Observation Window Dimentions (-X,+X,-Y,+Y) [Rad]: "
                                            << H_angle_min << ","
                                            << H_angle_max << ","
                                            << V_angle_min << ","
                                            << V_angle_max
                                            << '\n'
                                            << "Simulation Resolutoin: "
                                            << NUM_RAYS_X
                                            << " x "
                                            << NUM_RAYS_Y
                                            << '\n';

            Metric_Parameters_type &Parameters = this->p_Initial_Conditions->Metric_Parameters;

            switch (e_metric) {

            case Kerr:

                Image_Output_files[Image_order] << "Spin Parameter: " << Parameters.Spin
                    << '\n';
                break;

            case Wormhole:

                Image_Output_files[Image_order] << "Spin Parameter: " << Parameters.Spin << ", Redshift Parameter: " << Parameters.Redshift_Parameter
                    << '\n';
                break;

            case Reg_Black_Hole:

                Image_Output_files[Image_order] << "Parameter: " << Parameters.RBH_Parameter
                    << '\n';
                break;

            case Naked_Singularity:

                Image_Output_files[Image_order] << "Gamma: " << Parameters.JNW_Gamma_Parameter
                    << '\n';
                break;

            case Gauss_Bonnet:

                Image_Output_files[Image_order] << "Gamma: " << Parameters.GB_Gamma_Parameter
                    << '\n';
                break;

            case BH_w_Dark_Matter:

                Image_Output_files[Image_order] << "Mass_Halo: " << Parameters.Halo_Mass << ", Halo Compactness: " << Parameters.Compactness
                    << '\n';
                break;
            }
        }

            Image_Output_files[Image_order] << "Image X Coord [M],"
                                            << " "
                                            << "Image Y Coord [M],"
                                            << " "
                                            << "Novikov-Thorne Disk Redshift [-],"
                                            << " "
                                            << "Novikov-Thorne Flux [M^-2],"
                                            << " "
                                            << "Synchotron Intensity I [Jy/sRad],"
                                            << " "
                                            << "Synchotron Intensity Q [Jy/sRad],"
                                            << " "
                                            << "Synchotron Intensity U [Jy/sRad],"
                                            << " "
                                            << "Synchotron Intensity V [Jy/sRad],"
                                            << " ";
            if (Active_Sim_Mode == 2) {

                Image_Output_files[Image_order] << "Source r Coord [M],"
                                                << " "
                                                << "Source phi Coord [Rad],"
                                                << " " 
                                                << "Radial Momentum (covariant),"
                                                << " "
                                                << "Theta Momentum (covariant),"
                                                << " "
                                                << "Phi Momentum (covariant),"
                                                << " ";

                switch (e_metric) {

                case Kerr:

                    Image_Output_files[Image_order] << "Spin Parameter";
                    break;

                case Wormhole:

                    Image_Output_files[Image_order] << "Spin Parameter,"
                        << " "
                        << "Redshift Parameter";
                    break;

                case Reg_Black_Hole:

                    Image_Output_files[Image_order] << "Parameter";
                    break;

                case Naked_Singularity:

                    Image_Output_files[Image_order] << "Gamma";
                    break;

                case Gauss_Bonnet:

                    Image_Output_files[Image_order] << "Gamma";
                    break;

                case BH_w_Dark_Matter:

                    Image_Output_files[Image_order] << "Halo Mass,"
                                                    << " " 
                                                    << "Halo Compactness";
                    break;
                }
            }

            Image_Output_files[Image_order] << '\n';
    }
}

void File_manager_class::open_image_output_files(int Sim_mode_2_number) {

    const int File_number = SPACETIME_NUMBER;

    // Create the path to the main results directory
    std::filesystem::path dir(std::filesystem::current_path() / "Sim_Results"); 

    // Init the std::path variables where we will store the names of the output files
    std::filesystem::path Image_file_names[File_number];
    std::filesystem::path Momenta_file_names[File_number];

    // Init the std::path variables of the full file paths
    std::filesystem::path Image_full_path[ORDER_NUM]{};
    std::filesystem::path Momenta_full_path[ORDER_NUM]{};

    // Specify the output file extention
    std::filesystem::path file_extention(".txt");

    // Set weather we truncate the file upon opening or not
    auto open_type = std::ios::app;

    if (Truncate_files) {

        open_type = std::ios::trunc;

    }

    // Loop over all the files and populate the (so far empty) 
    for (int File_Index = 0; File_Index <= ORDER_NUM - 1; File_Index += 1) {

        Image_file_names[File_Index] = Base_File_Names[e_metric]
                                     + "_n"
                                     + std::to_string(File_Index);

        Image_file_names[File_Index].replace_extension(file_extention);
        Image_full_path[File_Index] = dir / Image_file_names[File_Index];
       

        Image_Output_files[File_Index].open(Image_full_path[File_Index], open_type);

    }

    if (Truncate_files) {

        // If we are truncating the file, we should write the metadata to it
        File_manager_class::write_simulation_metadata(Sim_mode_2_number);

    }

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

    File_manager_class::write_simulation_metadata(int(0));
}

void File_manager_class::write_image_data_to_file(Results_type* s_Ray_results) {

    for (int Image_order = direct; Image_order <= ORDER_NUM - 1; Image_order += 1) {

        Image_Output_files[Image_order] << s_Ray_results->Image_Coords[x]
                                        << " "
                                        << s_Ray_results->Image_Coords[y]
                                        << " "
                                        << s_Ray_results->Redshift_NT[Image_order]
                                        << " "
                                        << s_Ray_results->Flux_NT[Image_order]
                                        << " "
                                        << s_Ray_results->Intensity[Image_order][I] * CGS_TO_JANSKY
                                        << " "
                                        << s_Ray_results->Intensity[Image_order][Q] * CGS_TO_JANSKY
                                        << " "
                                        << s_Ray_results->Intensity[Image_order][U] * CGS_TO_JANSKY
                                        << " "
                                        << s_Ray_results->Intensity[Image_order][V] * CGS_TO_JANSKY;

        if (Active_Sim_Mode == 2) {

            Image_Output_files[Image_order] << " "
                                            << s_Ray_results->Source_Coords[e_r][Image_order]
                                            << " "
                                            << s_Ray_results->Source_Coords[e_phi][Image_order]
                                            << " "
                                            << s_Ray_results->Three_Momentum[e_r][Image_order]
                                            << " "
                                            << s_Ray_results->Three_Momentum[e_theta][Image_order]
                                            << " "
                                            << s_Ray_results->Three_Momentum[e_phi][Image_order]
                                            << " ";
                switch (e_metric) {

                case Kerr:
                    Image_Output_files[Image_order] << s_Ray_results->Parameters.Spin;
                    break;

                case Wormhole:
                    Image_Output_files[Image_order] << s_Ray_results->Parameters.Spin
                                                    << " " 
                                                    << s_Ray_results->Parameters.Redshift_Parameter;
                    break;

                case Reg_Black_Hole:
                    Image_Output_files[Image_order] << s_Ray_results->Parameters.RBH_Parameter;
         
                    break;

                case Naked_Singularity:
                    Image_Output_files[Image_order] << s_Ray_results->Parameters.JNW_Gamma_Parameter;
                    break;

                case Gauss_Bonnet:
                    Image_Output_files[Image_order] << s_Ray_results->Parameters.GB_Gamma_Parameter;
                    break;

                case BH_w_Dark_Matter:
                    Image_Output_files[Image_order] << s_Ray_results->Parameters.Halo_Mass << " " << s_Ray_results->Parameters.Compactness;
                    break;
                }
        }

        Image_Output_files[Image_order] << '\n';

    }
}

void File_manager_class::log_photon_path(Results_type* s_Ray_results) {

    //for (int log_index = 0; log_index <= s_Ray_results->Ray_log.size() - 1; log_index++) {

    //    for (int state_index = e_r; state_index <= e_State_Number - 1; state_index++) {

    //        Log_Output_File << s_Ray_results->Ray_log[log_index][state_index] << " ";
    //      
    //    }

    //    Log_Output_File << '\n';
    //}
   
};

void File_manager_class::close_image_output_files() {

    for (int File_Index = 0; File_Index <= ORDER_NUM - 1; File_Index++) {

        Image_Output_files[File_Index].close();
    }
}

void File_manager_class::close_log_output_file() {

    Log_Output_File.close();

}