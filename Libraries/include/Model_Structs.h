#pragma once

#include "Enumerations.h"
#include <string>

struct Numerical_disk_params_type {

	double Polytrope_coeff;
	double Polytrope_index;

	long long int Radial_grid_size;
	long long int Theta_grid_size;

	double* Radial_grid;
	double* Theta_grid;

	double* Raw_density_data;
	double* Raw_B_field_r_data;
	double* Raw_B_field_theta_data;
	double* Raw_B_field_phi_data;

	/* @brief String with the full path to the disk density XML. */
	std::string Density_file_path;

	/* @brief String with the full path to the disk magnetic field XML. */
	std::string Mag_field_file_path;

	Spline_selection_enums e_Spline_type;

};

struct Common_RIAF_params_type {

	double Disk_opening_angle;

	double Density_power_law_scale;

	double Density_power_law_power;

	double Temperature_power_law_scale;

	double Temperature_power_law_power;

	double Density_cutoff_radius;

	double Density_cutoff_scale;

	double Temperature_cutoff_radius;

	double Temperature_cutoff_scale;

};

struct Colab_test_1_params_type {

	double Radial_scale;

	double Vertical_scale;

};

struct Hotspot_model_params_type {

	double Density_power_law_scale;

	double Density_power_law_power;

	/*! Standard deviation of the Gaussian density profile. */
	double Density_gaussian_spread;

	double Temperature_power_law_scale;

	double Temperature_power_law_power;

	/*! Standard deviation of the Gaussian temperature profile. */
	double Temperature_gaussian_spread;

	/*! Standard deviation of the Gaussian temporal profile. */
	double Temporal_gaussian_spread;

	/*! Radius of the hotspot. Only affects the Spherical profile. */
	double Radius;

	/*! Coordinate time of maximum hotspot density */
	double Coord_time_offset;

};

struct Novikov_Thorne_params_type {

	double r_in;

	double r_out;

};

struct EOS_params_type {

	double Polytrope_Coeff;

	double Polytrope_Power;

};