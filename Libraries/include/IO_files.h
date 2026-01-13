#pragma once
#include <filesystem>
#include <string>
#include <fstream>
#include "Enumerations.h"
#include "Structs.h"

struct Results_type;

class File_manager_class {

private:

	std::string Base_File_Names[e_Spacetime_number] = {

		"Kerr",
		"Wormhole",
		"Regular_Black_Hole",
		"Janis_Newman_Winicour",
		"Einstein_Gauss_Bonnet",
		"BH_w_Dark_Matter_Halo",
		"Numerical",
		"Minkowski"

	};

	std::ofstream* Output_File;
	std::filesystem::path Output_File_Path;

	Initial_conditions_type *p_Initial_Conditions;

	bool Truncate_files{};

	void write_observer_metadata(std::ofstream* Output_file);
	void write_metric_metadata(std::ofstream* Output_file);
	void write_accretion_disk_metadata(std::ofstream* Output_file);
	void write_hotspot_metadata(std::ofstream* Output_file);
	void write_emission_models_metadata(std::ofstream* Output_file);
	void write_integrator_metadata(std::ofstream* Output_file);

	void write_simulation_metadata();
	void write_debug_metadata();

	public:

		int sim_mode_2_ray_number;

		File_manager_class(Initial_conditions_type* p_Initial_Conditions);

		void create_output_file();

		void write_image_data_to_file(Results_type* Ray_results);

		void write_debug_data_to_file(Debug_mode_struct* Debug_results);

		void log_photon_path(Results_type* s_Ray_results);

		void close_output_file();

		void open_output_file();

		void get_geodesic_data(double J_data[], double p_theta_data[]);

};

