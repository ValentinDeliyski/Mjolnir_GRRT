#pragma once

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

struct Page_Thorne_params_type {

	double r_in;

	double r_out;

};