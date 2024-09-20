#pragma once

#ifndef GENERAL_MATH_FUNCTIONS

	#define GENERAL_MATH_FUNCTIONS

	int mat_vec_multiply_4D(double const Matrix[4][4], double const Vector[4], double result[4]);

	double vector_norm(double Vector[], int Vector_size);

	double my_max(double const vector[], int const element_number);

	class Novikov_Thorne_Model;

	bool interpolate_crossing(double State_Vector[], double Old_State_Vector[], double Crossing_coords[], double crossing_momenta[], Novikov_Thorne_Model* const NT_Model);

	double dot_product(double vector_1[3], double vector_2[3]);

	void convert_spherical_to_cartesian(double* Spherical_Coords, double* Cartesian_Coords);

	void convert_cartesian_to_spherical(double* Cartesian_Coords, double* Spherical_Coords);

#endif