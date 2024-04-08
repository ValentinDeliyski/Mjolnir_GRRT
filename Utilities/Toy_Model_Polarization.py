import matplotlib
import matplotlib.pyplot as plt

from numpy import zeros, array, dot, cross, column_stack, append
from numpy import sqrt, cos, sin, arctan2, arctan
from numpy import ndarray
from numpy import pi

from numpy.linalg import norm, inv

from Support_functions import Spacetimes
from Support_functions.Parsers import Simulation_Parser

def get_polarization_vector(Photon_Momentum: ndarray, 
                            Fluid_Velocity: ndarray, 
                            B_Field: ndarray, 
                            Image_Coords: ndarray, 
                            Metric: ndarray) -> ndarray:

    """
    The math source for this implementation is from https://arxiv.org/pdf/2206.09455.pdf.

    INPUT DESCRIPTION:

    -- Fluid_Velocity has components [x, y, z], such that x is orientated along the "-r" direction, y along "phi",
    and z along the "theta" direction. Units are in multiples of "c".
    
    """

    # =========== ENUMS =========== #

    e_x = 0
    e_y = 1
    e_z = 2

    e_t     = 0
    e_r     = 1
    e_theta = 2
    e_phi   = 3

    e_Momentum     = 0
    e_Image_Coords = 1
    e_Metric       = 2

    # ============================== #

    Polarization_Vectors = []
    Scaled_Polarization_Vectors = []

    for _, Phoron_Log_Entry in enumerate(zip(Photon_Momentum, Image_Coords, Metric)):

        Photon_ZAMO_Momentum = zeros(4)

        Photon_ZAMO_Momentum[e_t] = 1 / sqrt(-Phoron_Log_Entry[e_Metric][e_t])

        for index in range(1, 4):
            Photon_ZAMO_Momentum[index] = Phoron_Log_Entry[e_Momentum][index - 1] / sqrt(Phoron_Log_Entry[e_Metric][index])

        # ========================================== Construct the boost matrix ========================================== #

        beta_param  = norm(Fluid_Velocity)
        gamma_param = 1 / sqrt(1 - beta_param**2)

        Boost_Matrix = array([[          gamma_param,                                      -gamma_param * Fluid_Velocity[e_x],                    0,                      - gamma_param * Fluid_Velocity[e_y]                             ],
                              [-gamma_param * Fluid_Velocity[e_x], 1 + (gamma_param - 1) * Fluid_Velocity[e_x] * Fluid_Velocity[e_x] / beta_param**2, 0,     (gamma_param - 1) * Fluid_Velocity[e_x] * Fluid_Velocity[e_y] / beta_param**2],
                              [              0,                                                            0,                                      1,                                              0                                      ],
                              [-gamma_param * Fluid_Velocity[e_y],     (gamma_param - 1) * Fluid_Velocity[e_y] * Fluid_Velocity[e_x] / beta_param**2, 0, 1 + (gamma_param - 1) * Fluid_Velocity[e_y] * Fluid_Velocity[e_y] / beta_param**2]])

        # ============================ Boost the photon momentum from ZAMO to the fluid frame ============================ #

        Photon_4_Momentum_Fluid = dot(Boost_Matrix, Photon_ZAMO_Momentum)
        Photon_3_Momentum_Fluid = Photon_4_Momentum_Fluid[1:]

        # ============================= Evaluate the polarization vector in the fluid frame ============================== #

        Polarization_3_Vector_Fluid = cross(Photon_3_Momentum_Fluid, B_Field) / norm(Photon_3_Momentum_Fluid)
        Polarization_4_Vector_Fluid = append(0, Polarization_3_Vector_Fluid)

        # ============================= Boost the polarization vector back to the ZAMO frame ============================= #

        inv_Boost_Matrix = inv(Boost_Matrix)
        Polarization_4_Vector_ZAMO = dot(inv_Boost_Matrix, Polarization_4_Vector_Fluid)

        # ====================================== Evaluate Penrose-Walker constants ======================================= #

        kappa_1 = sqrt(Phoron_Log_Entry[e_Metric][e_phi]) * (Photon_ZAMO_Momentum[e_t]     * Polarization_4_Vector_ZAMO[e_r]   - Photon_ZAMO_Momentum[e_r]   * Polarization_4_Vector_ZAMO[e_t])
        kappa_2 = sqrt(Phoron_Log_Entry[e_Metric][e_phi]) * (Photon_ZAMO_Momentum[e_theta] * Polarization_4_Vector_ZAMO[e_phi] - Photon_ZAMO_Momentum[e_phi] * Polarization_4_Vector_ZAMO[e_theta])

        # ========================================= Evaluate the scaling factors ========================================= #

        Redshift                 = 1 / Photon_4_Momentum_Fluid[e_t]
        Projected_Disk_Thickness = abs(Photon_4_Momentum_Fluid[e_t] / Photon_4_Momentum_Fluid[e_theta])

        # ============================ Evaluate the parallel transported polarization vector ============================= #

        """ 
        Polarization_Vector_ZAMO already has a factor of sin(zeta) in it, from the cross product - equation (16), which translates to 
        Transported_Polarization_Vector trough the Penrose-Walker constant "kappa". So equation (33) is satisfied (even though it does not look like it).

        """
        Transported_Polarization_Vector = zeros(2)
        Transported_Polarization_Vector[0] = Redshift**2 * sqrt(Projected_Disk_Thickness) * ( Phoron_Log_Entry[e_Image_Coords][e_x] * kappa_1 + Phoron_Log_Entry[e_Image_Coords][e_y] * kappa_2) / norm(Phoron_Log_Entry[e_Image_Coords])**2
        Transported_Polarization_Vector[1] = Redshift**2 * sqrt(Projected_Disk_Thickness) * (-Phoron_Log_Entry[e_Image_Coords][e_x] * kappa_2 + Phoron_Log_Entry[e_Image_Coords][e_y] * kappa_1) / norm(Phoron_Log_Entry[e_Image_Coords])**2

        # ========================== Scale the polarization vector for visualization purposes =========================== #

        Scaled_Transported_Polarization_Vector = norm(Transported_Polarization_Vector) * Transported_Polarization_Vector

        # ============================ Append everything to a list that I return at the end ============================= #

        Polarization_Vectors.append(Transported_Polarization_Vector)
        Scaled_Polarization_Vectors.append(Scaled_Transported_Polarization_Vector)

    return array(Polarization_Vectors), array(Scaled_Polarization_Vectors)

def get_observational_quantities(Polarization_Vectors, Image_Coordinates):

    def split_EVPA(EVPA):

        DISCONTINUITY_TRESHOLD = pi / 2

        branch_index = []
        branch_index.append(0)

        for index, _ in enumerate(EVPA):

            if index > 1 and abs(EVPA[index] - EVPA[index - 1]) > DISCONTINUITY_TRESHOLD:

                branch_index.append(index)

        branch_index.append(len(EVPA) - 1)

        return branch_index

    EVPA = []

    Polarization_I = []
    Polarization_U = []
    Polarization_Q = []

    Image_phi_coord = []

    for f_vector, Coordinates in zip(Polarization_Vectors, Image_Coordinates):

        EVPA.append(arctan(-f_vector[0] / f_vector[1]))

        Polarization_I.append(norm(f_vector)**2)
        Polarization_U.append(-2 * f_vector[0] * f_vector[1])
        Polarization_Q.append(f_vector[1]**2 - f_vector[0]**2)

        Coord_angle = arctan2(Coordinates[1], Coordinates[0])
        if Coord_angle < 0:
            Coord_angle = Coord_angle + 2 * pi      

        Image_phi_coord.append(Coord_angle)

    EVPA_Branches = split_EVPA(EVPA = EVPA)

    return array(EVPA[:-1]), array(EVPA_Branches), array(Polarization_I[:-1]), array(Polarization_Q[:-1]), array(Polarization_U[:-1]), array(Image_phi_coord[:-1])

if __name__ == '__main__':
 
    Active_spacetime = "Schwarzshild"

    Schw_Sim_Path = "C:\\Users\\Valur\\Documents\\Repos\\Gravitational_Lenser\\Sim_Results\\Kerr_n0"
    Schw_Parser   = Simulation_Parser(Schw_Sim_Path)

    Other_Metric_Sim_Path = "C:\\Users\\Valur\\Documents\\Repos\\Gravitational_Lenser\\Sim_Results\\Wormhole_n0"
    Other_Metric_Parser   = Simulation_Parser(Other_Metric_Sim_Path)

    """  
    For metrics with more than one parameter these go as:
        * Wormhole: 1 = Spin, 2 = Redshift
        * Black Hole With Dark Matter Halo: 1 - Halo Mass, 2 = Halo Compactness
    """

    Other_Metric_Param_Number = 2
    Param_Sweep_Number = int(len(Other_Metric_Parser.X_coords) / Other_Metric_Parser.Photon_Number)

    # ======================================= Initial Conditions ======================================= #

    B_Field = array([0.5, 0, 0.87]) # This vector has components [r, theta, phi]

    beta_angle     = -120 / 180 * pi
    Fluid_Velocity = 0.3 * array([cos(beta_angle), sin(beta_angle), 0]) # This vector has components [x, y, z]

    # ====================================== Ray-Tracer Log Parser ====================================== #

    Schw_Photon_Momentum         = column_stack((Schw_Parser.Radial_Momentum, Schw_Parser.Theta_Momentum, Schw_Parser.Phi_Momentum))
    Schw_Image_Coordiantes       = column_stack((Schw_Parser.X_coords, Schw_Parser.Y_coords))
    Schw_Disk_Intersection_Point = column_stack((Schw_Parser.Source_R_Coord, Schw_Parser.Source_Phi_Coord))

    Other_Metric_Photon_Momentum         = column_stack((Other_Metric_Parser.Radial_Momentum, Other_Metric_Parser.Theta_Momentum, Other_Metric_Parser.Phi_Momentum))
    Other_Metric_Image_Coordiantes       = column_stack((Other_Metric_Parser.X_coords, Other_Metric_Parser.Y_coords))
    Other_Metric_Disk_Intersection_Point = column_stack((Other_Metric_Parser.Source_R_Coord, Other_Metric_Parser.Source_Phi_Coord))

    # ==================================== Polarization Calculations ==================================== #

    Schwarzschild_class  = Spacetimes.Schwarzschild()
    Schwarzschild_metric = Schwarzschild_class.metric(r = Schw_Parser.Source_R_Coord, theta = pi / 2)

    Figure = plt.figure()
    Tick_plot       = Figure.add_subplot(151)
    Intensity_Plot  = Figure.add_subplot(152)
    Delta_I_Plot    = Figure.add_subplot(153)
    EVPA_Plot       = Figure.add_subplot(154)
    Delta_EVPA_Plot = Figure.add_subplot(155)

    Schw_Polarization_Vectors, Schw_Scaled_Polzarization_Vectors = get_polarization_vector(Photon_Momentum = Schw_Photon_Momentum,
                                                                                           Fluid_Velocity = Fluid_Velocity,
                                                                                           B_Field = B_Field,
                                                                                           Image_Coords = Schw_Image_Coordiantes,
                                                                                           Metric = Schwarzschild_metric)

    Schw_EVPA, Schw_EVPA_Branches, Schw_Polarization_I, Schw_Polarization_Q, Schw_Polarization_U, Schw_Image_phi_coord = get_observational_quantities(Polarization_Vectors = Schw_Polarization_Vectors,
                                                                                                                                                      Image_Coordinates    = Schw_Image_Coordiantes)
    
    """ Plot the reference values for Schwarzschild """

    Tick_plot.quiver(Schw_Parser.X_coords[0::10] - Schw_Scaled_Polzarization_Vectors.T[0][0::10] / 2,
                     Schw_Parser.Y_coords[0::10] - Schw_Scaled_Polzarization_Vectors.T[1][0::10] / 2,
                     Schw_Scaled_Polzarization_Vectors.T[0][0::10],
                     Schw_Scaled_Polzarization_Vectors.T[1][0::10],
                     headwidth = 0,
                     headlength = 0,
                     headaxislength = 0,
                     angles = 'xy', 
                     scale_units = 'xy', 
                     scale = 1,
                     color = 'k')

    Intensity_Plot.plot(Schw_Image_phi_coord, Schw_Polarization_I, color = "k")

    for index in range(len(Schw_EVPA_Branches) - 1):

        EVPA_Plot.plot(Schw_Image_phi_coord[Schw_EVPA_Branches[index] : Schw_EVPA_Branches[index + 1] - 1], Schw_EVPA[Schw_EVPA_Branches[index] : Schw_EVPA_Branches[index + 1] - 1], color = 'k')
        

    """ Big fat loop over the parameter sweep of the other metric """

    for index in range(Param_Sweep_Number):

        """ Read the parameter value from the ray-tracer log """

        Param_Sweep_Value = Other_Metric_Parser.Param_1[index * Other_Metric_Parser.Photon_Number]
        Param_Sweep_Min_Value = Other_Metric_Parser.Param_1[0]
        Param_Sweep_Max_Value = Other_Metric_Parser.Param_1[-1]

        if Other_Metric_Param_Number == 2:
            Param_Sweep_Value = Other_Metric_Parser.Param_2[index * Other_Metric_Parser.Photon_Number]
            Param_Sweep_Min_Value = Other_Metric_Parser.Param_2[0]
            Param_Sweep_Max_Value = Other_Metric_Parser.Param_2[-1]

        cmap = matplotlib.colormaps['plasma']
        Color_range = cmap((Param_Sweep_Value - Param_Sweep_Min_Value) / (Param_Sweep_Max_Value - Param_Sweep_Min_Value))

        """ The classes need to be in here because they get init with thir parameter values """

        WH   = Spacetimes.Wormhole(r_throat = 1, parameter = Param_Sweep_Value)
        RBH  = Spacetimes.Regular_Black_Hole(parameter = Param_Sweep_Value)
        JNW  = Spacetimes.JNW_Naked_Singularity(parameter = Param_Sweep_Value)
        GBNS = Spacetimes.Gaus_Bonnet_Naked_Singularity(parameter = Param_Sweep_Value)

        Spacetime_dict ={"Wormhole":     WH,
                         "RBH":          RBH,
                         "JNW":          JNW,
                         "Gauss_Bonnet": GBNS}
        
        Other_Metric = Spacetime_dict[Other_Metric_Parser.metric].metric(r = Other_Metric_Parser.Source_R_Coord, theta = pi / 2)

        Polarization_Vectors, Scaled_Polzarization_Vectors = get_polarization_vector(Photon_Momentum = Other_Metric_Photon_Momentum[index * Other_Metric_Parser.Photon_Number : (index + 1) * Other_Metric_Parser.Photon_Number],
                                                                                     Fluid_Velocity = Fluid_Velocity,
                                                                                     B_Field = B_Field,
                                                                                     Image_Coords = Other_Metric_Image_Coordiantes[index * Other_Metric_Parser.Photon_Number : (index + 1) * Other_Metric_Parser.Photon_Number],
                                                                                     Metric = Other_Metric[index * Other_Metric_Parser.Photon_Number : (index + 1) * Other_Metric_Parser.Photon_Number])
       

        EVPA, EVPA_Branches, Polarization_I, Polarization_Q, Polarization_U, Image_phi_coord = get_observational_quantities(Polarization_Vectors = Polarization_Vectors,
                                                                                                                            Image_Coordinates = Other_Metric_Image_Coordiantes[index * Other_Metric_Parser.Photon_Number : (index + 1) * Other_Metric_Parser.Photon_Number])

        Tick_plot.quiver(Other_Metric_Parser.X_coords[index * Other_Metric_Parser.Photon_Number : (index + 1) * Other_Metric_Parser.Photon_Number: 10] - Scaled_Polzarization_Vectors.T[0][0::10] / 2,
                         Other_Metric_Parser.Y_coords[index * Other_Metric_Parser.Photon_Number : (index + 1) * Other_Metric_Parser.Photon_Number: 10] - Scaled_Polzarization_Vectors.T[1][0::10] / 2,
                         Scaled_Polzarization_Vectors.T[0][0::10],
                         Scaled_Polzarization_Vectors.T[1][0::10],
                         headwidth = 0,
                         headlength = 0,
                         headaxislength = 0,
                         angles = 'xy', 
                         scale_units = 'xy', 
                         scale = 1,
                         color = Color_range)
        
        Intensity_Plot.plot(Image_phi_coord, Polarization_I, color = Color_range)

        for index in range(len(EVPA_Branches) - 1):

            EVPA_Plot.plot(Image_phi_coord[EVPA_Branches[index] : EVPA_Branches[index + 1] - 1], EVPA[EVPA_Branches[index] : EVPA_Branches[index + 1] - 1], color = Color_range)

        Delta_I_Plot.plot(Image_phi_coord, Polarization_I - Schw_Polarization_I, color = Color_range)

        Delta_EVPA = EVPA - Schw_EVPA

        for index, _ in enumerate(Delta_EVPA):

            if Delta_EVPA[index] > pi / 2:
                Delta_EVPA[index] = Delta_EVPA[index] - pi

            elif Delta_EVPA[index] < -pi / 2:
                Delta_EVPA[index] = Delta_EVPA[index] + pi

        Delta_EVPA_Plot.plot(Image_phi_coord, Delta_EVPA, color = Color_range)

    plt.show()