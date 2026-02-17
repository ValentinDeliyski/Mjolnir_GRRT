#pragma once
#include "General_math_functions.h"
#include "Enumerations.h"
#include "Structs.h"
#include "Spacetimes.h"
#include "Constants.h"
#include <iostream>
#include <complex>

void get_derivative_of_inverse_metric(const double Metric[4][4], const double deriv_Metric[4][4], double deriv_inv_Metric[4][4]);

/*! @brief Inverts the covariant metric tensor. The current implementation works only for axi-symmetric metrics
 *  in coordinates in which the cross terms appears in the [0][3] and [3][0] elements
 *
 *   \param [out] Inv_metric - Pointer to the 2D array that hold the inverse metric.
 *   \param [in] Metric - Pointer to the 2D array that holds the metric.
 *   \return Nothing.
 */
void invert_metric(double inv_metric[4][4], const double metric[4][4]);

/*! @brief Computes the determinant of the metric tensor. The current implementation works only for axi-symmetric metrics
 *  in coordinates in which the cross terms appears in the [0][3] and [3][0] elements
 *
 *   \param [in] Metric - Pointer to the 2D array that hold the metric.
 *   \return The determinant of the metric tensor.
 */
double get_metric_det(const double metric[4][4]);

/*! @brief Computes the determinant of the induced equatorial metric tensor. The current implementation works only for axi-symmetric metrics
 *  in coordinates in which the cross terms appears in the [0][3] and [3][0] elements
 *
 *   \param [in] Metric - Pointer to the 2D array that hold the metric.
 *   \return The determinant of the induced equatorial metric tensor.
 */
double get_eq_induced_metric_det(const double metric[4][4]);

/*! @brief Computes the norm of a 4-vector "Vector" with respect to the metric "Metric".
 *
 *   \param [in] Vector - Pointer to the vector.
 *   \param [in] Metric - Pointer to the 2D array that holds the metric.
 *   \param [in] Vector_type - Enum that determines the tensor type of the vector "Vector" (covariant or contravariant)
 *   \return The norm of the 4-vector.
 */
double get_complex_4vec_norm(const std::complex<double>* const Vector, const double Metric[4][4], Tensor_type_enums Vector_type);

void Manipulate_index(const Metric_type* const p_Metric, const double* const Vector, double* const Result, Index_manipulation_enums Index_operation);

void Normalize_complex_vector(std::complex<double>* const Vector, const double Metric[4][4], Tensor_type_enums Vector_type);

/*! @brief Converts a contravariant vector from the coordinate basis to the ZAMO basis.
 *
 *   \param [in] Metric - Pointer to the metric struct.
 *   \param [in] Contravariant_Vector - Pointer to the contravariant vector.
 *   \param [out] ZAMO_Vector - Pointer to the contravariant ZAMO vector.
 *   \return Nothing.
 */
template<typename Vec_type>
void Contravariant_coord_to_ZAMO(const Metric_type* const p_Metric, const Vec_type* const Contravariant_Vector, Vec_type* const ZAMO_Vector) {

    ZAMO_Vector[e_t] = p_Metric->Lapse_function * Contravariant_Vector[e_t];
    ZAMO_Vector[e_r] = sqrt(p_Metric->Metric[e_r][e_r]) * Contravariant_Vector[e_r];
    ZAMO_Vector[e_theta] = -sqrt(p_Metric->Metric[e_theta][e_theta]) * Contravariant_Vector[e_theta];
    ZAMO_Vector[e_phi] = sqrt(p_Metric->Metric[e_phi][e_phi]) * (Contravariant_Vector[e_phi] - p_Metric->Shift_function * Contravariant_Vector[e_t]);

};

/*! @brief Converts a contravariant vector from the ZAMO basis to the coordinate basis.
 *
 *   \param [in] Metric - Pointer to the metric struct.
 *   \param [in] ZAMO_Vector - Pointer to the contravariant ZAMO vector.
 *   \param [out] Contravariant_Vector - Pointer to the contravariant vector.
 *   \return Nothing.
 */
template<typename Vec_type>
void ZAMO_to_Contravariant_coord(const Metric_type* const p_Metric, const Vec_type* const ZAMO_Vector, Vec_type* const Contravariant_Vector){

    Contravariant_Vector[e_t] = ZAMO_Vector[e_t] / p_Metric->Lapse_function;
    Contravariant_Vector[e_r] = ZAMO_Vector[e_r] / sqrt(p_Metric->Metric[e_r][e_r]);
    Contravariant_Vector[e_theta] = -ZAMO_Vector[e_theta] / sqrt(p_Metric->Metric[e_theta][e_theta]);
    Contravariant_Vector[e_phi] = p_Metric->Shift_function / p_Metric->Lapse_function * ZAMO_Vector[e_t] + ZAMO_Vector[e_phi] / sqrt(p_Metric->Metric[e_phi][e_phi]);

};

/*! @brief Computes the initial 4-momentum of a ray, basaed on its direction angles at the observer, and stores it in the Initial Conditions struct.
 *
 *   \param [out] p_Initial_Conditions - Pointer to the Initial Conditions struct.
 *   \param [in] V_angle_cam - Vertical direction angle.
 *   \param [in] H_angle_cam - Horizontal direction angle.
 *   \return Nothing.
 */
void get_intitial_conditions_from_angles(Initial_conditions_type* p_Initial_Conditions, double V_angle, double H_angle);

/*! @brief Computes the coordinates of the image on the observation plane, and stores them in the Initial Conditions struct.
 *
 *   \param [in] p_Initial_Conditions - Pointer to the Initial Conditions struct.
 *   \param [out] Image_coords - Pointer to the array that hold the image coordinates on the observation plane.
 *   \return Nothing.
 */
void get_image_coordinates(Initial_conditions_type* p_Initial_Conditions, double* const Image_coords);

/*! @brief Computes redshift for for a ray, emitted by a source with 4-velocity "U_source", at a point specified by "State_Vector",
 *  for an observer, specified by the "Observer" class instance.
 *
 *   \param [in] State_Vector - Pointer to the ray state vector, containing covarian momentum components.
 *   \param [in] U_source - Pointer to the contravariant 4-velocity of the emmiting element.
 *   \param [in] p_Sim_Context - Pointer to the simulation context struct
 *   \return The redshift of the ray, as measured by the observer.
 */
double get_redshift(const double* const State_Vector, const double* const U_source, const Simulation_Context_type* const p_Sim_Context);

/*! @brief Computes a Lorentz boost matrix, given a contravariant 4-velocity (in an observers basis) U_cource, then stores it in Boost_matrix.
 *
 *   \param [out] Boost_matrix - Pointer to 2D array that holds the Lorentz matrix.
 *   \param [in] U_source - Pointer to the contravariant 4-velocity of the emmiting element.
 *	 \param [in] Inverse_boost - Boolean that determines weather the inverse boost matrix is calculated.
 *   \return Nothing.
 */
void get_Lorentz_boost_matrix(double Boost_matrix[4][4], const double* const U_source, bool Inverse_boost);

/*! @brief Checks if a turning point in the theta coordinate has been passed between the states "State_Vector" and "Old_State"
 *
 *   \param [in] State_Vector - Pointer to the current ray state vector.
 *   \param [in] Old_State - Pointer to the previous ray state vector.
 *   \return Boolian flag that states weather a turning point has been passed.
 */
bool Check_for_theta_turning_point(const double* const State_Vector, const double* const Old_State);

bool Check_for_equatorial_crossing(const double* const State_Vector, const double* const Old_State);

/*! @brief Computes the order of a ray, given its initial conditions and the current number of turning points in the theta coordinate the ray has passed.
 *
 *   \param [in] N_theta_turning_points - The number of turning points in the theta coordinate the ray has passed
 *   \param [in] p_Initial_Conditions - Pointer to the Initial Contitions struct.
 *   \return The order of the ray
 */
int compute_image_order(const int N_theta_turning_points, const int N_equatorial_crossings, Initial_conditions_type* const p_Initial_Conditions);

/*! @brief Computes the connection coefficients, given the metric and its radial and theta derivatives, then stores it in the 3D array "Connectrion_Coeffs"
 *
 *   \param [in] s_Metric - The struct that holds the metric.
 *   \param [in] s_dr_metric - The struct that holds the radial metric derivative.
 *   \param [in] s_dtheta_metric - The struct that holds the theta metric derivative.
 *   \param [out] Connectrion_Coeffs - The 3D array that holds the connection coefficients
 *   \return Nothing.
 */
void get_connection_coefficients(const Metric_type s_Metric, const Metric_type s_dr_metric, const Metric_type s_dtheta_metric, double Connectrion_Coeffs[4][4][4]);

void get_initial_conditions_from_image_coords(Initial_conditions_type* p_Initial_Conditions, double Image_X_coord, double Image_Y_coord);

std::complex<double> get_Penrose_Walker_constant(const double* const State_Vector, const Simulation_Context_type* const p_Sim_Context, const std::complex<double>* const Polarization_Vector);