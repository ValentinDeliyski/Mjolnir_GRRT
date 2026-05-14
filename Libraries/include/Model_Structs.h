#pragma once

#include "Enumerations.h"
#include <string>
#include <memory>

struct Numerical_disk_params_type {

	double Polytrope_coeff;
	double Polytrope_index;

	size_t Z_coord_grid_size;
	size_t R_coord_grid_size;

	std::shared_ptr<double[]> Z_coord_grid;
	std::shared_ptr<double[]> R_coord_grid;

	std::shared_ptr<double[]> Density_data;

	/* @brief String with the full path to the disk density XML. */
	std::string Density_file_path;

	Spline_selection_enums e_Spline_type;

};

struct Phenomenological_mag_field_params_type {

	/*! @brief Specifies the direction of the magnetic field. */
	Magnetic_field_geometry_enums e_Mag_field_geometry{};

	/*! @brief Specifies how the magnitude of the magnetic field is calculated. */
	Magnetic_field_magnitude_enums e_Mag_field_magnitude_profile{};

	/*! @brief The disk magnetization value [-]. */
	double Magnetization{};

	/*! @brief The constant magnetic field geometry in the Eularian frame.
	   The components are specified as [B_r, B_theta, B_phi].
	   This vector gets normalized with the metric before use. */
	double Mag_field_geometry[3]{};

	double Mag_field_B_0{};

	double Mag_field_power{};

	double Mag_field_r_0{};

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

	/*! @brief The peak density value in [1 / cm^3]. */
	double Electron_density_scale{};

	/*! @brief The peak temperature value in [K]. */
	double Electron_temperature_scale{};

	Phenomenological_mag_field_params_type Mag_field_params{};

};

struct Colab_test_1_params_type {

	double Density_scale;

	double Radial_scale;

	double Vertical_scale;

	Phenomenological_mag_field_params_type Mag_field_params{};

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

	double Mag_field_geometry[3];

};

struct EOS_params_type {

	double Polytrope_Coeff;

	double Polytrope_Power;

};