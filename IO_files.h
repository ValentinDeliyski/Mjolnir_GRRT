#pragma once

#ifndef IO_FILES

	#define IO_FILES

	#include <filesystem>
	#include <string>
	#include <fstream>

	int open_output_files(e_Spacetimes e_metric, std::ofstream data[], std::ofstream momentum_data[], bool truncate);

	int close_output_files(std::ofstream data[], std::ofstream momentum_data[]);

	int get_geodesic_data(double J_data[], double p_theta_data[], int* Data_number);

	int write_to_file(double Image_coordiantes[], double redshift, double Flux, double State_vector[], double parameter, double J,
					  int Image_oder, bool lens_from_file, std::ofstream data[], std::ofstream momentum_data[]);

#endif
