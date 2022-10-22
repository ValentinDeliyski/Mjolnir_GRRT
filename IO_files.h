#pragma once

#ifndef IO_FILES

	#define IO_FILES

	#include <filesystem>
	#include <string>
	#include <fstream>
	#include "Lensing.h"

	int open_output_files(std::ofstream data[], std::ofstream momentum_data[]);

	int close_output_files(std::ofstream data[], std::ofstream momentum_data[]);

	int get_geodesic_data(double J_data[], double p_theta_data[], int* Data_number);

	int write_to_file(results Ray_results, std::ofstream data[], std::ofstream momentum_data[]);

#endif
