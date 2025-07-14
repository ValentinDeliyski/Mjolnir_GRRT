/*! @file General_math_functions.h
 *	This file declares all the general (non-spacetime specific) supporting math functions, used throughout the code .
 */

#pragma once
#define _USE_MATH_DEFINES
#include "Enumerations.h"
#include <cmath>

 /*! @defgroup General_Math_Support_Functions General Math Support Functions
  *  This group contains all the general (non-spacetime specific) supporting math functions, used throughout the code.
  *  @{
  */

/*! @brief Computes the Euclidian norm of the vector "Vector", with "Vector_size" number of elements.
 *
 *  @param [in] Vector - Pointer to the vector.
 *  @param [in] Vector_size - The number of elements in the vector.
 *  @return The Euclidian norm of the vector "Vector".
 */
double vector_norm(const double* const Vector, const int Vector_size);

/*! @brief Multiplies two 4D matricies.
 *
 *  @param [in] Mat_A - The matrix that gets multiplied on the left.
 *  @param [in] Mat_B - The matrix that gets multiplied on the right.
 *	@param [out] Result - the matrix multiplication Mat_A * Mat_B.
 *  @return Nothing.
 */
void matrix_matrix_multiply(const double Mat_A[4][4], const double Mat_B[4][4], double Result[4][4]);

/*! @brief Main Multiplies a 4D vector by a 4x4 matrix.
 *  Multiplies the 4D vector "Vector" by the 4x4 matrix "Matrix", and stores the result in the vector "Result".
 *
 *  @param [in] Matrix - The 4x4 matrix, represented as a 2D array.
 *  @param [in] Vector - The 4D vector to be multiplied by the marix "Matrix".
 *  @param [out] Result - The result of the multiplication.
 *  @return Nothing.
 */
void mat_vec_multiply_4D(double const Matrix[4][4], const double* const Vector, double* const Result);

/*! @brief Returns the largest by absolute value element of a vector
 *  Returns the largest by absolute value element the vector "Vector". This is really intended to be used in
 *  the step controller for the Dormond-Prince integrator, so the 0-th element of the vector is ignored
 *  (it corresponds to the coordinate time of the ray, and becomes large in the emitting region - this throws off
 *  the adaptive step calculation).
 *
 *  @param [in] Vector - Pointer to the vector.
 *  @param [in] Element_number - The number of elements in the vector.
 *  @return The unsigned maximum by absolute value element.
 */
double get_max_element(const double* const Vector, const int Element_number);

/*! @brief Returns the largest by absolute value element of relative state error.
 *  Returns the largest by absolute value element of relative state error. This is intended to be used in
 *  the step controller for the Dormond-Prince integrator.
 *
 *  @param [in] error_state - Pointer to the state error vector.
 *  @param [in] current_state - Pointer to the state vector.
 *  @return The unsigned maximum by absolute value relative state error.
 */
double get_max_relative_error(const double* const Error_state, const double* const Current_state);

/*! @brief Interpolates the coordinates and three-momentum of the ray, at which it crosses the equator.
 *
 *  @param [in] State_Vector - Pointer to the current state vector.
 *  @param [in] Old_State_Vector - Pointer to the previous state vector.
 *  @param [out] Crossing_State - Pointer to interpolated state vector at the equator crossing point.
 *  @return A boolian flag for weather the equator has been crossed or not.
 */
bool interpolate_crossing(const double* const State_Vector, 
						  const double* const Old_State_Vector,  
						  double* const Crossing_State);

/*! @brief Computes a simple Eucliduan dot product between two vectors with element numbers "Vector_size".
 *
 *  @param [in] vector_1 - Pointer to the first vector.
 *  @param [in] vector_2 - Pointer to the second vector.
 *  @param [in] Vector_size - The number of elements in each vector.
 *  @return The dot product between the two vectors.
 */
double dot_product(const double* const Vector_1, const double* const Vector_2, int Vector_size);

/*! @brief Converts vector from spherical coordinates to cartesian.
 *	Converts vector, stored in an array pointed to by "Spherical_Coords", from spherical coordinates to cartesian, and stores
 *	the results in the array, pointed to by "Cartesian_Coords".
 *
 *  @param [in] Spherical_Coords - Pointer to the vector expressed in spherical coordinates.
 *  @param [out] Cartesian_Coords - Pointer to the vector expressed in carrtesian coordinates.
 *  @return Nothing.
 */
void convert_spherical_to_cartesian(const double* const Spherical_Coords, double* const Cartesian_Coords);
 
/*! @brief Converts vector from cartesian coordinates to spherical.
 *  Converts vector, stored in an array pointed to by "Cartesian_Coords", from cartesian coordinates to spherical, and stores
 *  the results in the array, pointed to by "Spherical_Coords".
 *
 *  @param [in] Cartesian_Coords - Pointer to the vector expressed in carrtesian coordinates.
 *  @param [out] Spherical_Coords - Pointer to the vector expressed in spherical coordinates.
 *  @return Nothing.
 */
void convert_cartesian_to_spherical(const double* const Cartesian_Coords, double* const Spherical_Coords);

/*! @brief Adds two vectors.
 *	Adds two vectors, stored in the arrays pointed to by "vec_1" and "vec_2", and store the result in the 
 *  array pointed to by "Result"
 *
 *  @param [in] Vec_1 - Pointer to the first vector.
 *  @param [in] Vec_2 - Pointer to the second vector.
 *  @param [in] size - the number of components in the vector.
 *  @param [out] Result - Pointer to the result vector
 *  @return Nothing.
 */
void add_vectors(const double* const vec_1, const double* const vec_2, const int size, double* const Result);

/*! @brief Raises a base to an integer power.
 *	Raaises the "base" variable to the integer power "exponent". 
 *  NOTE: The "exponent" variable is assumed positive.
 *
 *  @param [in] base - The number being exponentiated.
 *  @param [in] exponent - The power to which "base" is being raised.
 *  @return Result - base^exponent.
 */
double int_power(const double base, const int exponent);

/** @} */ // End of the General_Math_Support_Functions group