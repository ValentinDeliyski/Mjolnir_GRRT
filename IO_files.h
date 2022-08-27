#pragma once

#include <filesystem>
#include <string>

int open_output_files(Spacetimes e_metric, std::ofstream data[], std::ofstream momentum_data[], bool truncate) {

	int const File_number = sizeof(File_Names) / sizeof(std::string);

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

int write_to_file(double Image_coordiantes[], double redshift, double Flux, double State_vector[], double parameter, double J,
				  int Image_oder, bool lens_from_file, std::ofstream data[], std::ofstream momentum_data[]) {


	switch (Image_oder) {

		case direct:

			data[0] << -Image_coordiantes[0] << " " << Image_coordiantes[2] << " " << redshift << " " << std::fixed << std::setprecision(15) << Flux << " " << std::setprecision(6) << State_vector[2] << " " << State_vector[0] << '\n';

			if (lens_from_file) {

				momentum_data[0] << State_vector[5] << " " << State_vector[4] << " " << J << " " << parameter << '\n';

			}

			break;

		case first:

			data[1] << -Image_coordiantes[0] << " " << Image_coordiantes[2] << " " << redshift << " " << std::fixed << std::setprecision(15) << Flux << " " << std::setprecision(6) << State_vector[2] << " " << State_vector[0] << '\n';

			if (lens_from_file) {

				momentum_data[1] << State_vector[5] << " " << State_vector[4] << " " << J << " " << parameter << '\n';

			}

			break;

		case second:

			data[2] << -Image_coordiantes[0] << " " << Image_coordiantes[2] << " " << redshift << " " << std::fixed << std::setprecision(15) << Flux << " " << std::setprecision(6) << State_vector[2] << " " << State_vector[0] << '\n';

			if (lens_from_file) {

				momentum_data[2] << State_vector[5] << " " << State_vector[4] << " " << J << " " << parameter << '\n';

			}

			break;

		case third:

			data[3] << -Image_coordiantes[0] << " " << Image_coordiantes[2] << " " << redshift << " " << std::fixed << std::setprecision(15) << Flux << " " << std::setprecision(6) << State_vector[2] << " " << State_vector[0] << '\n';

			if (lens_from_file) {

				momentum_data[3] << State_vector[5] << " " << State_vector[4] << " " << J << " " << parameter << '\n';

			}
			
			break;

	}

	return OK;

}
