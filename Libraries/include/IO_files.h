#pragma once

#ifndef IO_FILES

	#define IO_FILES

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
			"RBH",
			"JNW",
			"Gauss_Bonnet",
			"BH_Dark_Matter_Halo"

		};

		std::ofstream Image_Output_files[ORDER_NUM];
		std::ofstream Log_Output_File;
		std::string Input_file_path_sim_mode_2;
		Initial_conditions_type *p_Initial_Conditions;

		bool Truncate_files{};

		public:

			File_manager_class(Initial_conditions_type* p_Initial_Conditions, std::string input_file_path, bool truncate);

			void open_image_output_files(int Sim_mode_2_number);

			void open_log_output_file();

			void write_image_data_to_file(Results_type* Ray_results);

			void write_simulation_metadata(int Sim_mode_2_number);

			void log_photon_path(Results_type* s_Ray_results);

			void close_image_output_files();

			void get_geodesic_data(double J_data[], double p_theta_data[], int* Data_number);

			void close_log_output_file();
	
	};

#endif
