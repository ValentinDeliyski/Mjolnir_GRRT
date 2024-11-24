#include "General_math_functions.h"

//! Computes the Euclidian norm of a vector.
/*! Computes the Euclidian norm of the vector "Vector", with "Vector_size" number of elements
 *
 *   \param [in] Vector - Pointer to the vector.
 *   \param [in] Vector_size - The number of elements in the vector.
 *   \return The Euclidian norm of the vector "Vector".
 */
double vector_norm(const double* const Vector, const int Vector_size) {

	double norm{};

	for (int index = 0; index <= Vector_size - 1; index++) {

		norm += Vector[index] * Vector[index];

	}

	norm = sqrt(norm);

	return norm;

};

//! Main Multiplies a 4D vector by a 4x4 matrix
/*! Multiplies the 4D vector "Vector" by the 4x4 matrix "Matrix", and stores the result in the vector "Result"
 *
 *   \param [in] Matrix - The 4x4 matrix, represented as a 2D array.
 *   \param [in] Vector - The 4D vector to be multiplied by the matxi "Matrix".
 *   \param [out] Result - The result of the multiplication.
 *   \return Nothing.
 */
void mat_vec_multiply_4D(double const Matrix[4][4], const double* const Vector, double* const Result) {

	for (int row = 0; row <= 3; row += 1) {

		Result[row] = 0.0;

		for (int column = 0; column <= 3; column += 1) {

			Result[row] += Matrix[row][column] * Vector[column];

		}
	}

};

//! Returns the largest by absolute value element of a vector
/*! Returns the largest by absolute value element the vector "Vector". This is really intended to be used in 
 *	the step controller for the Dormond-Prince integrator, so the 0-th element of the vector is ignored 
 *	(it corresponds to the coordinate time of the ray, and becomes large in the emitting region - this throws off
 *   the adaptive step calculation).
 *
 *   \param [in] Vector - Pointer to the vector.
 *   \param [in] Element_number - The number of elements in the vector.
 *   \return The unsigned maximum by absolute value element.
 */
double get_max_element(const double* const Vector, int const Element_number) {

	double max = fabs(Vector[e_r]);
	double max_candidate = 0;

	for (int index = e_theta; index <= Element_number; index += 1) {

		max_candidate = fabs(Vector[index]);

		if (max_candidate > max) {
			
			max = max_candidate;
		
		}

	}

	return max;
}

//! Returns the largest by absolute value element of relative state error.
/*! Returns the largest by absolute value element of relative state error. This is intended to be used in
 *	the step controller for the Dormond-Prince integrator.
 *
 *   \param [in] error_state - Pointer to the state error vector.
 *   \param [in] current_state - Pointer to the state vector.
 *   \return The unsigned maximum by absolute value relative state error.
 */
double get_max_relative_error(const double* const Error_state, const double* const Current_state) {

	double max_rel_error = fabs(Error_state[e_t] / Current_state[e_t]);

	for (int index = e_r; index <= e_State_Number - 2; index += 1) {

		double temp_error = fabs(Error_state[index] / (Current_state[index]));

		if (temp_error > max_rel_error && !isinf(temp_error)) {

			max_rel_error = temp_error;

		}
	}

	return max_rel_error;
}

//! Interpolates the coordinates and three-momentum of the ray, at which it crosses the equator.
/*! Interpolates the coordinates and three-momentum of the ray, at which it crosses the equator.
 *
 *   \param [in] State_Vector - Pointer to the current state vector.
 *   \param [in] Old_State_Vector - Pointer to the previous state vector.
 *	 \param [out] Crossing_coords - Pointer to the three-vector that holds the coordinates of the crossing point.
 *	 \param [out] Crossing_coords - Pointer to the three-vector that holds the three-momentum of the ray at the crossing point.
 *   \return A boolian flag for weather the equator has been crossed or not.
 */
bool interpolate_crossing(const double* const State_Vector, 
						  const double* const Old_State_Vector, 
						  double* const Crossing_coords, 
						  double* const Crossing_momenta) {

	// Check weather the equator has been crossed
	if ((State_Vector[e_theta] - M_PI_2) * (Old_State_Vector[e_theta] - M_PI_2) > 0)
	{

		memset(Crossing_coords, 0, 3 * sizeof(double));
		memset(Crossing_momenta, 0, 3 * sizeof(double));

		return false;

	}

	/*
	
	Interpolate the equatorial crossing coorinates
	
	*/

	double x = State_Vector[e_r] * sin(State_Vector[e_theta]) * cos(State_Vector[e_phi]);
	double y = State_Vector[e_r] * sin(State_Vector[e_theta]) * sin(State_Vector[e_phi]);
	double z = State_Vector[e_r] * cos(State_Vector[e_theta]);

	double x_old = Old_State_Vector[e_r] * sin(Old_State_Vector[e_theta]) * cos(Old_State_Vector[e_phi]);
	double y_old = Old_State_Vector[e_r] * sin(Old_State_Vector[e_theta]) * sin(Old_State_Vector[e_phi]);
	double z_old = Old_State_Vector[e_r] * cos(Old_State_Vector[e_theta]);

	double gradient[3]   = { x - x_old, y - y_old, z - z_old };
	double const_term[3] = { x_old, y_old, z_old };

	double crossing_param = -const_term[2] / gradient[2];

	for (int index = 0; index <= 2; index++) {

		Crossing_coords[index] = gradient[index] * crossing_param + const_term[index];

	}

	/*
	
	Interpolate the covariant momenta at those coorinates
	
	*/

	double momentum_param = (M_PI_2 - Old_State_Vector[e_theta]) / (State_Vector[e_theta] - Old_State_Vector[e_theta]);

	for (int index = e_t; index <= e_phi; index++) {

		Crossing_momenta[index] = momentum_param * State_Vector[e_p_t + index] + (1 - momentum_param) * Old_State_Vector[e_p_t + index];
	}
	
	return true;
}

//! Computes a simple Eucliduan dot product between two vectors with element numbers "Vector_size".
/*! Computes a simple Eucliduan dot product between two vectors with element numbers "Vector_size".
 *
 *   \param [in] vector_1 - Pointer to the first vector.
 *   \param [in] vector_2 - Pointer to the second vector.
 *	 \param [in] Vector_size - The number of elements in each vector.
 *   \return The dot product between the two vectors.
 */
double dot_product(const double* const Vector_1, const double* const Vector_2, int Vector_size) {

	double result{};

	for (int index = 0; index <= Vector_size - 1; index++) {

		result += Vector_1[index] * Vector_2[index];

	}

	return  result;

}

//! Converts vector from spherical coordinates to cartesian.
/*! Converts vector, stored in an array pointed to by "Spherical_Coords", from spherical coordinates to cartesian, and stores 
 *  the results in the vector, pointed to by "Cartesian_Coords".
 *
 *   \param [in] Spherical_Coords - Pointer to the vector expressed in spherical coordinates.
 *   \param [in] Cartesian_Coords - Pointer to the vector expressed in carrtesian coordinates.
 *   \return Nothing.
 */
void convert_spherical_to_cartesian(double* Spherical_Coords, double* Cartesian_Coords) {

	double sin_theta = sin(Spherical_Coords[e_theta]);
	double cos_theta = cos(Spherical_Coords[e_theta]);

	double sin_phi = sin(Spherical_Coords[e_phi]);
	double cos_phi = cos(Spherical_Coords[e_phi]);

	Cartesian_Coords[x] = Spherical_Coords[e_r] * sin_theta * cos_phi;
	Cartesian_Coords[y] = Spherical_Coords[e_r] * sin_theta * sin_phi;
	Cartesian_Coords[z] = Spherical_Coords[e_r] * cos_theta;

}

//! Converts vector from spherical coordinates to cartesian.
/*! Converts vector, stored in an array pointed to by "Cartesian_Coords", from cartesian coordinates to spherical, and stores
 *  the results in the vector, pointed to by "Spherical_Coords".
 *
 *   \param [in] Cartesian_Coords - Pointer to the vector expressed in carrtesian coordinates.
 *   \param [in] Spherical_Coords - Pointer to the vector expressed in spherical coordinates.
 *   \return Nothing.
 */
void convert_cartesian_to_spherical(double* Cartesian_Coords, double* Spherical_Coords) {

	Spherical_Coords[e_r]  = Cartesian_Coords[x] * Cartesian_Coords[x];
	Spherical_Coords[e_r] += Cartesian_Coords[y] * Cartesian_Coords[y];
	Spherical_Coords[e_r] += Cartesian_Coords[z] * Cartesian_Coords[z];
	Spherical_Coords[e_r]  = sqrt(Spherical_Coords[e_r]);

	Spherical_Coords[e_theta] = acos(Cartesian_Coords[z] / Spherical_Coords[e_r]);

	Spherical_Coords[e_phi] = atan2(Cartesian_Coords[y], Cartesian_Coords[x]);

}

//! Adds two 4-vectors.
/*! Adds two 4-vectors.
 *
 *   \param [in] Vec_1 - Pointer to the first vector.
 *   \param [in] Vec_2 - Pointer to the second vector.
 *	 \param [out] Result - Pointer to the result vector
 *   \return Nothing.
 */
void add_4_vectors(const double* const Vec_1, const double* const Vec_2, double* const Result) {

	for (int idx = 0; idx <= 3; idx++) {

		Result[idx] = Vec_1[idx] + Vec_1[idx];

	}

}