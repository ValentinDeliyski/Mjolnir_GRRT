#pragma once
#include <filesystem>
#include <string>
#include <fstream>
#include "Enumerations.h"
#include "Structs.h"

struct Results_type;

class File_manager_class {

private:

	std::string Base_File_Names[SPACETIME_NUMBER] = {

		"Kerr",
		"Wormhole",
		"Regular_Black_Hole",
		"Janis_Newman_Winicour",
		"Einstein_Gauss_Bonnet",
		"BH_w_Dark_Matter_Halo"

	};

	std::ofstream Image_Output_files[ORDER_NUM];
	std::ofstream Log_Output_File;
	Initial_conditions_type *p_Initial_Conditions;

	bool Truncate_files{};

	public:

		int sim_mode_2_ray_number;

		File_manager_class(Initial_conditions_type* p_Initial_Conditions);

		void open_image_output_files();

		void open_log_output_file();

		void write_image_data_to_file(Results_type* Ray_results);

		void write_simulation_metadata();

		void log_photon_path(Results_type* s_Ray_results);

		void close_image_output_files();

		void get_geodesic_data(double J_data[], double p_theta_data[]);

		void close_log_output_file();

};

