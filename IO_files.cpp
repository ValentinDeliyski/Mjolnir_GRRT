#pragma once

#define _USE_MATH_DEFINES

#include <filesystem>
#include <string>
#include <fstream>
#include <iostream>

#include "Enumerations.h"
#include "IO_files.h"
#include "Constants.h"

extern e_Spacetimes e_metric;
extern Const_bool lens_from_file;
extern Const_bool truncate;

extern Const_Float r_obs;
extern Const_Float theta_obs;

std::string File_Names[] = {

	"Kerr_Data0",
	"Kerr_Data1",
	"Kerr_Data2",
	"Kerr_Data3",

	"Kerr_Momentum_Data0",
	"Kerr_Momentum_Data1",
	"Kerr_Momentum_Data2",
	"Kerr_Momentum_Data3",

	"Wormhole_Data0",
	"Wormhole_Data1",
	"Wormhole_Data2",
	"Wormhole_Data3",

	"Wormhole_Momentum_Data0",
	"Wormhole_Momentum_Data1",
	"Wormhole_Momentum_Data2",
	"Wormhole_Momentum_Data3",

	"RBH_Data0",
	"RBH_Data1",
	"RBH_Data2",
	"RBH_Data3",

	"RBH_Momentum_Data0",
	"RBH_Momentum_Data1",
	"RBH_Momentum_Data2",
	"RBH_Momentum_Data3",

    "JNW_Data0",
	"JNW_Data1",
	"JNW_Data2",
	"JNW_Data3",

	"JNW_Momentum_Data0",
	"JNW_Momentum_Data1",
	"JNW_Momentum_Data2",
	"JNW_Momentum_Data3"

};

int open_output_files(std::ofstream data[], std::ofstream momentum_data[]) {

	int const File_number = SPACETIME_NUMBER * 8;

	std::filesystem::path dir("C:\\Users\\Valentin\\Documents\\University stuff\\General Relativity");
	std::filesystem::path file[File_number];
	std::filesystem::path full_path[File_number];
	std::filesystem::path file_extention(".txt");

	for (int File_Index = 0; File_Index <= File_number - 1; File_Index += 1) {

		file[File_Index] = File_Names[File_Index];

		file[File_Index].replace_extension(file_extention);

		full_path[File_Index] = dir / file[File_Index];

	}

	auto open_type = std::ios::app;

	if (truncate) {

		open_type = std::ios::trunc;

	}

	for (int File_Index = 0; File_Index <= ORDER_NUM - 1; File_Index += 1) {

		switch (e_metric) {

		case Kerr:

			data[File_Index].open(full_path[File_Index + 0 * 4], open_type);

			momentum_data[File_Index].open(full_path[File_Index + 1 * 4], open_type);

			break;

		case Wormhole:

			data[File_Index].open(full_path[File_Index + 2 * 4], open_type);

			momentum_data[File_Index].open(full_path[File_Index + 3 * 4], open_type);

			break;

		case Reg_Black_Hole:

			data[File_Index].open(full_path[File_Index + 4 * 4], open_type);

			momentum_data[File_Index].open(full_path[File_Index + 5 * 4], open_type);

			break;

		case Naked_Singularity:

			data[File_Index].open(full_path[File_Index + 6 * 4], open_type);

			momentum_data[File_Index].open(full_path[File_Index + 7 * 4], open_type);

			break;

		}
	}

	for (int Image_order = direct; Image_order <= ORDER_NUM - 1; Image_order += 1) {

		data[Image_order] << "Observer Radial Position [M]: " << r_obs << '\n'
						  << "Observer Inclination [Deg]: " << theta_obs * 180 / M_PI << '\n'
						  << "Image X Coord [M],"
						  << " "
  		   				  << "Image Y Coord [M],"
	  					  << " "
						  << "Novikov-Thorne Disk Redshift [-],"
						  << " "
						  << "Novikov-Thorne Flux Redshift [?],"
						  << " "
						  << "Optically Thin Disk Intensity [Jy/sRad],"
						  << " "
						  << "Novikov-Thorne Disk Phi Coord [Rad],"
						  << " "
						  << "Novikov-Thorne Disk Radial Coord [M]"
						  << '\n';

		if (lens_from_file) {

			momentum_data[Image_order] << "Radial Momentum (covariant),"
				<< " "
				<< "Theta Momentum (covariant),"
				<< " "
				<< "Phi Momentum (covariant),"
				<< " "
				<< "Metric Parameter,"
				<< '\n';

		}
	}

	return OK;

}

int close_output_files(std::ofstream data[], std::ofstream momentum_data[]) {

	for (int i = 0; i <= 3; i++) {

		data[i].close();
		momentum_data[i].close();

	}

	return OK;

}

int get_geodesic_data(double J_data[], double p_theta_data[], int* Data_number) {

	std::ifstream geodesic_data;

	geodesic_data.open("C:\\Users\\Valentin\\Documents\\University stuff\\General Relativity\\Polarization\\Schwarzschild_Impact_parameters\\First_relativistic\\geodesic_data_20_deg_Sch_r6_198_photons.txt", std::ios::in);

	while (true) {

		geodesic_data >> J_data[*Data_number] >> p_theta_data[*Data_number];

		if (geodesic_data.eof() == true) {

			break;

		}

		*Data_number += 1;

	}

	geodesic_data.close();

	return OK;
}

int write_to_file(results Ray_results, std::ofstream data[], std::ofstream momentum_data[]) {


	for (int Image_order = direct; Image_order <= ORDER_NUM - 1; Image_order += 1) {

		data[Image_order] << Ray_results.Image_Coords[0]
						  << " "
						  << Ray_results.Image_Coords[1]
						  << " "
						  << Ray_results.Redshift_NT[Image_order]
						  << " "
						  << Ray_results.Flux_NT[Image_order]
						  << " "
						  << Ray_results.Intensity[Image_order]
						  << " "
						  << std::setprecision(6) << Ray_results.Source_Coords[e_phi][Image_order]
						  << " "
						  << Ray_results.Source_Coords[e_r][Image_order]
						  << '\n';

		if (lens_from_file) {

			momentum_data[Image_order] << Ray_results.Three_Momentum[e_r][Image_order]
									   << " "
									   << Ray_results.Three_Momentum[e_theta][Image_order]
									   << " "
									   << Ray_results.Three_Momentum[e_phi][Image_order]
									   << " "
									   << Ray_results.Parameters[e_metric]
									   << '\n';

		}
	}
	return OK;

}

