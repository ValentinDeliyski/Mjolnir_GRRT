#pragma once

#ifndef GENERAL_MATH_FUNCTIONS

	#define GENERAL_MATH_FUNCTIONS

	int mat_vec_multiply_4D(double Matrix[4][4], double Vector[4], double result[4]);

	double vector_norm(double Vector[], int Vector_size);

	double my_max(double vector[]);

	bool interpolate_crossing(double State_Vector[], double Old_State_Vector[], double Crossing_coords[]);

	double dot_product(double vector_1[3], double vector_2[3]);

#endif