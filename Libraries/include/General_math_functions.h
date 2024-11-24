#pragma once
#define _USE_MATH_DEFINES
#include "Enumerations.h"
#include <cmath>

//! Computes the Euclidian norm of a vector.
/*!  Computes the Euclidian norm of the vector "Vector", with "Vector_size" number of elements
 *
 *   \param [in] Vector - Pointer to the vector.
 *   \param [in] Vector_size - The number of elements in the vector.
 *   \return The Euclidian norm of the vector "Vector".
 */
double vector_norm(const double* const Vector, const int Vector_size);

//! Main Multiplies a 4D vector by a 4x4 matrix
/*! Multiplies the 4D vector "Vector" by the 4x4 matrix "Matrix", and stores the result in the vector "Result"
 *
 *   \param [in] Matrix - The 4x4 matrix, represented as a 2D array.
 *   \param [in] Vector - The 4D vector to be multiplied by the matxi "Matrix".
 *   \param [out] Result - The result of the multiplication.
 *   \return Nothing.
 */
void mat_vec_multiply_4D(double const Matrix[4][4], const double* const Vector, double* const Result);

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
double get_max_element(const double* const Vector, const int Element_number);

//! Returns the largest by absolute value element of relative state error.
/*! Returns the largest by absolute value element of relative state error. This is intended to be used in
 *	the step controller for the Dormond-Prince integrator.
 *
 *   \param [in] error_state - Pointer to the state error vector.
 *   \param [in] current_state - Pointer to the state vector.
 *   \return The unsigned maximum by absolute value relative state error.
 */
double get_max_relative_error(const double* const Error_state, const double* const Current_state);

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
						  double* const crossing_momenta);

//! Computes a simple Eucliduan dot product between two vectors with element numbers "Vector_size".
/*! Computes a simple Eucliduan dot product between two vectors with element numbers "Vector_size".
 *
 *   \param [in] vector_1 - Pointer to the first vector.
 *   \param [in] vector_2 - Pointer to the second vector.
 *	 \param [in] Vector_size - The number of elements in each vector.
 *   \return The dot product between the two vectors.
 */
double dot_product(const double* const Vector_1, const double* const Vector_2, int Vector_size);

//! Converts vector from spherical coordinates to cartesian.
/*! Converts vector, stored in an array pointed to by "Spherical_Coords", from spherical coordinates to cartesian, and stores
 *  the results in the vector, pointed to by "Cartesian_Coords".
 *
 *   \param [in] Spherical_Coords - Pointer to the vector expressed in spherical coordinates.
 *   \param [in] Cartesian_Coords - Pointer to the vector expressed in carrtesian coordinates.
 *   \return Nothing.
 */
void convert_spherical_to_cartesian(double* Spherical_Coords, double* Cartesian_Coords);

//! Converts vector from spherical coordinates to cartesian.
/*! Converts vector, stored in an array pointed to by "Cartesian_Coords", from cartesian coordinates to spherical, and stores
 *  the results in the vector, pointed to by "Spherical_Coords".
 *
 *   \param [in] Cartesian_Coords - Pointer to the vector expressed in carrtesian coordinates.
 *   \param [in] Spherical_Coords - Pointer to the vector expressed in spherical coordinates.
 *   \return Nothing.
 */
void convert_cartesian_to_spherical(double* Cartesian_Coords, double* Spherical_Coords);

//! Adds two 4-vectors.
/*! Adds two 4-vectors.
 *
 *   \param [in] Vec_1 - Pointer to the first vector.
 *   \param [in] Vec_2 - Pointer to the second vector.
 *	 \param [out] Result - Pointer to the result vector
 *   \return Nothing.
 */
void add_4_vectors(const double* const vec_1, const double* const vec_2, double* const result);
