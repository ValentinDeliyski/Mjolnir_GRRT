import sys
import os

""" Add the parent directory of this file to the search path, 
    so this script can be ran from the "Utilities" folder """
    
parent_directory = os.path.abspath('...')
sys.path.append(parent_directory)

from Support_functions.Parsers import Simulation_Parser
from Support_functions.Spacetimes_new import Kerr, Wormhole, Regular_Black_Hole, Gauss_Bonnet, Janis_Newman_Winicour, Spacetime
from numpy import array, sqrt, arctan2, dot, cross, append, zeros, sin, diag, pi, cos, flip
from numpy.linalg import inv, norm
from numpy.typing import NDArray

import matplotlib.pyplot as plt
from enum import Enum

class Coords(Enum):

    e_t = 0
    e_r = 1
    e_theta = 2
    e_phi = 3

    e_x = 0
    e_y = 1
    e_z = 2

class Polarization_class():
    
    def __init__(self, Spacetime_instance: Spacetime, Sim_metadata: dict):
        
        self.Spacetime_instance = Spacetime_instance
        self.Sim_metadata = Sim_metadata
        
    def get_Keplarian_ang_vel(self, dr_metric: NDArray):

        return (-dr_metric[Coords.e_t.value][Coords.e_phi.value] + sqrt((dr_metric[Coords.e_t.value][Coords.e_phi.value])**2 - dr_metric[Coords.e_t.value][Coords.e_t.value] * dr_metric[Coords.e_phi.value][Coords.e_phi.value])) / dr_metric[Coords.e_phi.value][Coords.e_phi.value]
        
    def get_disk_velocity(self, metric: NDArray, dr_metric: NDArray) -> NDArray:

        Angular_velocity = self.get_Keplarian_ang_vel(dr_metric)

        u_t = 1 / sqrt(-metric[Coords.e_t.value][Coords.e_t.value] - 2 * metric[Coords.e_t.value][Coords.e_phi.value] * Angular_velocity - metric[Coords.e_phi.value][Coords.e_phi.value] * Angular_velocity**2)

        return array([u_t, 0., 0., u_t * Angular_velocity])

    def Coord_to_ZAMO(self, Metric: NDArray, Coord_vec_contravariant: NDArray) -> NDArray:

        ZAMO_vec = array([0., 0., 0., 0.])

        Lapse = sqrt(Metric[Coords.e_t.value][Coords.e_phi.value]**2 / Metric[Coords.e_phi.value][Coords.e_phi.value] - Metric[Coords.e_t.value][Coords.e_t.value])
        Shift = -Metric[Coords.e_t.value][Coords.e_phi.value] / Metric[Coords.e_phi.value][Coords.e_phi.value]

        ZAMO_vec[Coords.e_t.value] = Lapse * Coord_vec_contravariant[Coords.e_t.value]
        ZAMO_vec[Coords.e_r.value] =  sqrt(Metric[Coords.e_r.value][Coords.e_r.value]) * Coord_vec_contravariant[Coords.e_r.value]
        ZAMO_vec[Coords.e_theta.value] = -sqrt(Metric[Coords.e_theta.value][Coords.e_theta.value]) * Coord_vec_contravariant[Coords.e_theta.value]
        ZAMO_vec[Coords.e_phi.value] = sqrt(Metric[Coords.e_phi.value][Coords.e_phi.value]) * (Coord_vec_contravariant[Coords.e_phi.value] - Shift * Coord_vec_contravariant[Coords.e_t.value])

        return ZAMO_vec

    def ZAMO_to_coord(self, Metric: NDArray, ZAMO_vec: NDArray) -> NDArray:

        Coord_vec = array([0., 0., 0., 0.])

        Lapse = sqrt(Metric[Coords.e_t.value][Coords.e_phi.value]**2 / Metric[Coords.e_phi.value][Coords.e_phi.value] - Metric[Coords.e_t.value][Coords.e_t.value])
        Shift = -Metric[Coords.e_t.value][Coords.e_phi.value] / Metric[Coords.e_phi.value][Coords.e_phi.value]

        Coord_vec[Coords.e_t.value] = ZAMO_vec[Coords.e_t.value] / Lapse
        Coord_vec[Coords.e_r.value] =  ZAMO_vec[Coords.e_r.value] / sqrt(Metric[Coords.e_r.value][Coords.e_r.value]) 
        Coord_vec[Coords.e_theta.value] = -ZAMO_vec[Coords.e_theta.value] / sqrt(Metric[Coords.e_theta.value][Coords.e_theta.value])
        Coord_vec[Coords.e_phi.value] = ZAMO_vec[Coords.e_phi.value] / sqrt(Metric[Coords.e_phi.value][Coords.e_phi.value]) + Shift / Lapse * ZAMO_vec[Coords.e_t.value]

        return Coord_vec
    
    def Manipulate_index(self, Metric: NDArray, Vector: NDArray, Operation: str):

        if Operation == "Raise":
            inv_Metric = inv(Metric)
            return inv_Metric.dot(Vector)
        
        elif Operation == "Lower":
            return Metric.dot(Vector)
        
        else:
            print("Unsupported index operation!")
            exit(1)   

    def get_boost_matrix(self, ZAMO_Disk_velocity: NDArray, inv: bool) -> NDArray:
        
        Normed_ZAMO_velocity = ZAMO_Disk_velocity / ZAMO_Disk_velocity[Coords.e_t.value]

        Beta_param = norm(Normed_ZAMO_velocity[Coords.e_r.value:])
        Gamma_param = 1 / sqrt(1 - Beta_param**2)

        if inv: 
            Beta_param = -Beta_param

        Boost_Matrix = zeros((4,4))

        Boost_Matrix[Coords.e_t.value][Coords.e_t.value] = Gamma_param

        Boost_Matrix[Coords.e_t.value][Coords.e_r.value]     = -Gamma_param * Normed_ZAMO_velocity[Coords.e_r.value]
        Boost_Matrix[Coords.e_t.value][Coords.e_theta.value] = -Gamma_param * Normed_ZAMO_velocity[Coords.e_theta.value]
        Boost_Matrix[Coords.e_t.value][Coords.e_phi.value]   = -Gamma_param * Normed_ZAMO_velocity[Coords.e_phi.value]

        Boost_Matrix[Coords.e_r.value][Coords.e_t.value]     = Boost_Matrix[Coords.e_t.value][Coords.e_r.value]
        Boost_Matrix[Coords.e_theta.value][Coords.e_t.value] = Boost_Matrix[Coords.e_t.value][Coords.e_theta.value]
        Boost_Matrix[Coords.e_phi.value][Coords.e_t.value]   = Boost_Matrix[Coords.e_t.value][Coords.e_phi.value]

        Boost_Matrix[Coords.e_r.value][Coords.e_r.value] = 1 + (Gamma_param - 1) * Normed_ZAMO_velocity[Coords.e_r.value]**2 / Beta_param**2
        Boost_Matrix[Coords.e_r.value][Coords.e_theta.value] = (Gamma_param - 1) * Normed_ZAMO_velocity[Coords.e_r.value] * Normed_ZAMO_velocity[Coords.e_theta.value] / Beta_param**2
        Boost_Matrix[Coords.e_r.value][Coords.e_phi.value]   = (Gamma_param - 1) * Normed_ZAMO_velocity[Coords.e_r.value] * Normed_ZAMO_velocity[Coords.e_phi.value] / Beta_param**2

        Boost_Matrix[Coords.e_theta.value][Coords.e_r.value] = Boost_Matrix[Coords.e_r.value][Coords.e_theta.value]
        Boost_Matrix[Coords.e_phi.value][Coords.e_r.value]   = Boost_Matrix[Coords.e_r.value][Coords.e_phi.value]

        Boost_Matrix[Coords.e_theta.value][Coords.e_theta.value] = 1 + (Gamma_param - 1) * Normed_ZAMO_velocity[Coords.e_theta.value]**2 / Beta_param**2
        Boost_Matrix[Coords.e_theta.value][Coords.e_phi.value] = (Gamma_param - 1) * Normed_ZAMO_velocity[Coords.e_theta.value] * Normed_ZAMO_velocity[Coords.e_phi.value] / Beta_param**2

        Boost_Matrix[Coords.e_phi.value][Coords.e_theta.value] = Boost_Matrix[Coords.e_theta.value][Coords.e_phi.value]

        Boost_Matrix[Coords.e_phi.value][Coords.e_phi.value] = 1 + (Gamma_param - 1) * Normed_ZAMO_velocity[Coords.e_phi.value]**2 / Beta_param**2

        return Boost_Matrix
    
    def get_PW_constant(self, Source_r_coord: float, P_ph_coord_con: NDArray, Pol_coord_con: NDArray) -> tuple[float, float]:
        
        if self.Spacetime_instance.Identify() == self.Spacetime_instance.Spacetime_enums.Kerr.value:
            
            kappa_1 = Source_r_coord * (P_ph_coord_con[Coords.e_t.value] * Pol_coord_con[Coords.e_r.value] - P_ph_coord_con[Coords.e_r.value] * Pol_coord_con[Coords.e_t.value])
            kappa_1 = kappa_1 + Source_r_coord * self.Spacetime_instance.a * (P_ph_coord_con[Coords.e_r.value] * Pol_coord_con[Coords.e_phi.value] - P_ph_coord_con[Coords.e_phi.value] * Pol_coord_con[Coords.e_r.value])
        
            kappa_2 = -Source_r_coord * ((Source_r_coord**2 + self.Spacetime_instance.a**2) * (P_ph_coord_con[Coords.e_phi.value] * Pol_coord_con[Coords.e_theta.value] - P_ph_coord_con[Coords.e_theta.value] * Pol_coord_con[Coords.e_phi.value]))
            kappa_2 = kappa_2 + Source_r_coord * self.Spacetime_instance.a * (P_ph_coord_con[Coords.e_t.value] * Pol_coord_con[Coords.e_theta.value] - P_ph_coord_con[Coords.e_theta.value] * Pol_coord_con[Coords.e_t.value])
            
        else:
            
            Metric = self.Spacetime_instance.get_metric(Source_coord, pi / 2)
            
            kappa_1 = sqrt(Metric[Coords.e_phi.value][Coords.e_phi.value] * (-Metric[Coords.e_r.value][Coords.e_r.value] * Metric[Coords.e_t.value][Coords.e_t.value]))
            kappa_1 = kappa_1 * (P_ph_coord_con[Coords.e_t.value] * Pol_coord_con[Coords.e_r.value] - P_ph_coord_con[Coords.e_r.value] * Pol_coord_con[Coords.e_t.value])
    
            kappa_2 = -sqrt(Metric[Coords.e_phi.value][Coords.e_phi.value])**3 * (P_ph_coord_con[Coords.e_phi.value] * Pol_coord_con[Coords.e_theta.value] - P_ph_coord_con[Coords.e_theta.value] * Pol_coord_con[Coords.e_phi.value])
    
        return kappa_1, kappa_2

    def get_polarization_vector(self,
                                Photon_Momentum_covariant: NDArray, 
                                Source_r_coord: float,
                                B_Field: list[float | int], 
                                Image_Coords: NDArray) -> NDArray:

        """ The math source for this implementation is from https://arxiv.org/pdf/2105.09440. """
        
        Metric_at_source = self.Spacetime_instance.get_metric(Source_r_coord, theta = pi / 2)
        dr_Metric_at_source = self.Spacetime_instance.get_dr_metric(Source_r_coord, theta = pi / 2)

        Photon_Momentum_covariant = append(-1, Photon_Momentum_covariant)

        P_ph_coord_con = self.Manipulate_index(Metric_at_source, Photon_Momentum_covariant, "Raise")
        P_ph_ZAMO = self.Coord_to_ZAMO(Metric_at_source, P_ph_coord_con)

        Disk_coord_velocity = self.get_disk_velocity(Metric_at_source, dr_Metric_at_source)
        Disk_ZAMO_velocity = self.Coord_to_ZAMO(Metric_at_source, Disk_coord_velocity)

        Boost_matrix = self.get_boost_matrix(Disk_ZAMO_velocity, inv = False)
        Photon_fluid_momentum = Boost_matrix.dot(P_ph_ZAMO)

        Pol_3_fluid = cross(Photon_fluid_momentum[Coords.e_r.value:], B_Field) / norm(Photon_fluid_momentum[Coords.e_r.value:])
        Pol_fluid = append(0, Pol_3_fluid)

        inv_Boost_matrix = inv(Boost_matrix)
        Pol_ZAMO = dot(inv_Boost_matrix, Pol_fluid)
        Pol_coord_con = self.ZAMO_to_coord(Metric_at_source, Pol_ZAMO)
        
        kappa_1, kappa_2 = self.get_PW_constant(Source_r_coord, P_ph_coord_con, Pol_coord_con)

        Redshift                 = 1 / Photon_fluid_momentum[Coords.e_t.value]
        Projected_Disk_Thickness = abs(Photon_fluid_momentum[Coords.e_t.value] / Photon_fluid_momentum[Coords.e_theta.value])

        if self.Spacetime_instance.Identify() == self.Spacetime_instance.Spacetime_enums.Kerr.value:

            alpha_coord = -(Image_Coords[Coords.e_x.value] + self.Spacetime_instance.a * sin(float(self.Sim_metadata["Observer Inclination [Deg]"]) / 180 * pi))
            beta_coord = Image_Coords[Coords.e_y.value]
            
        else:
            alpha_coord = -Image_Coords[Coords.e_x.value]
            beta_coord = Image_Coords[Coords.e_y.value]
            
        Transported_Polarization_Vector = zeros(2)
        Transported_Polarization_Vector[Coords.e_x.value] = (-alpha_coord * kappa_1 + beta_coord * kappa_2) / (alpha_coord**2 + beta_coord**2)
        Transported_Polarization_Vector[Coords.e_y.value] = ( alpha_coord * kappa_2 + beta_coord * kappa_1) / (alpha_coord**2 + beta_coord**2)

        return Transported_Polarization_Vector

if __name__ == "__main__":

    Sim_path = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\Utilities\\Reference_simulations\\Thin_Disk_Reference_Simulation\\Janis_Newman_Winicour"

    Sim_parser = Simulation_Parser(Sim_path)
    
    match Sim_parser.Simulation_metadata["Spacetime [-]"]:

        case "Kerr":
            Spacetime_instance = Kerr(float(Sim_parser.Simulation_metadata["Spin Parameter [M]"]))

        case "Wormhole":
            Spacetime_instance = Wormhole(r_throat = 1.0, parameter = float(Sim_parser.Simulation_metadata["Redshift Parameter [-]"]), use_global_coords = False)
        
        case "Einstein_Gauss_Bonnet":
            Spacetime_instance = Gauss_Bonnet(Param = float(Sim_parser.Simulation_metadata["Gamma [M^2]"]))
            
        case "Janis_Newman_Winicour":
            Spacetime_instance = Janis_Newman_Winicour(parameter = float(Sim_parser.Simulation_metadata["Gamma [-]"]))
            
        case _:
            exit(1)
            
    X_resolution: int = int(Sim_parser.Simulation_metadata["Simulation Resolution"].split(" ")[0])
    Y_resolution: int = int(Sim_parser.Simulation_metadata["Simulation Resolution"].split(" ")[2])

    Image_coords = array([Sim_parser.X_coords, Sim_parser.Y_coords]).T
    Photon_momentum = array([Sim_parser.Source_p_r, Sim_parser.Source_p_theta, Sim_parser.Source_p_phi]).T

    Analytical_Polarization_Vector_x = []
    Analytical_Polarization_Vector_y = []

    Numerical_Polarization_Vector_x = []
    Numerical_Polarization_Vector_y = []
    
    Polarization_instance = Polarization_class(Spacetime_instance, Sim_parser.Simulation_metadata)

    for P_photon, Source_coord, Image_coord, Flux, Numerical_pol_x, Numerical_pol_y in zip(Photon_momentum, Sim_parser.Source_r, Image_coords, Sim_parser.Disk_flux, Sim_parser.Polarization_vec_X, Sim_parser.Polarization_vec_Y):
    
        Numerical_Polarization_Vector_x.append(Numerical_pol_x)
        Numerical_Polarization_Vector_y.append(Numerical_pol_y)

        if (Flux != 0):
            
            Analytical_Polarization_Vector_x.append(Polarization_instance.get_polarization_vector(P_photon, Source_coord, [0.5, 0, 0.87], Image_coord)[Coords.e_x.value])
            Analytical_Polarization_Vector_y.append(Polarization_instance.get_polarization_vector(P_photon, Source_coord, [0.5, 0, 0.87], Image_coord)[Coords.e_y.value])

        else:
            Analytical_Polarization_Vector_x.append(0)
            Analytical_Polarization_Vector_y.append(0)

    Numerical_Polarization_Vector_x = flip(array(Numerical_Polarization_Vector_x).reshape(X_resolution, Y_resolution), axis = 0)
    Numerical_Polarization_Vector_y = flip(array(Numerical_Polarization_Vector_y).reshape(X_resolution, Y_resolution), axis = 0)

    Analytical_Polarization_Vector_x = flip(array(Analytical_Polarization_Vector_x).reshape(X_resolution, Y_resolution), axis = 0)
    Analytical_Polarization_Vector_y = flip(array(Analytical_Polarization_Vector_y).reshape(X_resolution, Y_resolution), axis = 0)

    # Set the X and Y axis limits, rescaling them for an observer, located at "Obs_effective_distance", rather than the simulation "Observer Distance [M]", and conver to to micro AS 
    axes_limits = Sim_parser.Simulation_metadata["Observation Window Dimentions (-X,+X,-Y,+Y) [M]"].split(",")
    axes_limits = array([float(Limit) for Limit in axes_limits])

    Fig = plt.figure()
    Comparison_plot_x = Fig.add_subplot(121)
    Comparison_plot_x.imshow(arctan2(Numerical_Polarization_Vector_x, Numerical_Polarization_Vector_y) - arctan2(Analytical_Polarization_Vector_x, Analytical_Polarization_Vector_y) , cmap = "seismic", vmin = -pi, vmax = pi, extent= tuple(axes_limits))

    Comparison_plot_y = Fig.add_subplot(122)
    Comparison_plot_y.imshow(Numerical_Polarization_Vector_y - Analytical_Polarization_Vector_y, cmap = "seismic", vmin = -1, vmax = 1)

    plt.show()

