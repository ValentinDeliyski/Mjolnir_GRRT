import matplotlib
import matplotlib.pyplot as plt

from numpy import zeros, array, dot, cross, column_stack, append, argmax, arange
from numpy import sqrt, cos, sin, arctan2, arctan
from numpy import ndarray
from numpy import pi

from numpy.linalg import norm, inv

from Support_functions import Spacetimes
from Support_functions.Parsers import Simulation_Parser

def split_Deltas(Delta):

    branch_index = []
    branch_index.append(0)

    for index, _ in enumerate(Delta):

        if index > 1 and Delta[index] * Delta[index - 1] < 0:

            branch_index.append(index)

    branch_index.append(len(Delta) - 1)

    return branch_index

def split_EVPA(EVPA):

    DISCONTINUITY_TRESHOLD = pi / 2

    branch_index = []
    branch_index.append(0)

    for index, _ in enumerate(EVPA):

        if index > 1 and abs(EVPA[index] - EVPA[index - 1]) > DISCONTINUITY_TRESHOLD:

            branch_index.append(index)

    branch_index.append(len(EVPA) - 1)

    return branch_index

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

    return array(Polarization_Vectors[:-1]), array(Scaled_Polarization_Vectors[:-1])

def get_observational_quantities(Polarization_Vectors, Image_Coordinates):

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

        Image_phi_coord.append(Coord_angle / pi)

    EVPA_Branches = split_EVPA(EVPA = EVPA)

    return array(EVPA), array(EVPA_Branches), array(Polarization_I), array(Polarization_Q), array(Polarization_U), array(Image_phi_coord)

def plot_polarization_ticks(Schw_Parser, Other_Metric_Parser, B_Fields, beta_angles, Other_Metric_Param_Number = 1, QUIVER_SAMPLE_SKIP = 20, Scalar = 1, Tick_plot = None, colormap = "plasma", Fontsize = 26):

    # ====================================== Ray-Tracer Log Parser ====================================== #

    Schw_Photon_Momentum   = column_stack((-Schw_Parser.Radial_Momentum, -Schw_Parser.Theta_Momentum, Schw_Parser.Phi_Momentum))
    Schw_Image_Coordiantes = column_stack((Schw_Parser.X_coords, Schw_Parser.Y_coords))

    Other_Metric_Photon_Momentum   = column_stack((-Other_Metric_Parser.Radial_Momentum, -Other_Metric_Parser.Theta_Momentum, Other_Metric_Parser.Phi_Momentum))
    Other_Metric_Image_Coordiantes = column_stack((Other_Metric_Parser.X_coords, Other_Metric_Parser.Y_coords))

    for _, (B_Field, beta_angle) in enumerate(zip(B_Fields, beta_angles)):

        Fluid_Velocity = 0.3 * array([cos(beta_angle), sin(beta_angle), 0]) # This vector has components [x, y, z]

        if Tick_plot == None:

            Figure_Pattern, (Tick_plot) = plt.subplots(1, 1, gridspec_kw = {'width_ratios': [2.4]}, constrained_layout = True)
            Figure_Pattern.set_figwidth(12)
            Figure_Pattern.set_figheight(25)

            colorbar_map = matplotlib.cm.ScalarMappable(cmap = cmap)
            colorbar_map.set_clim([0,3])
            
            Colorbar = Figure_Pattern.colorbar(colorbar_map, ax = Tick_plot)
            Colorbar.set_label(r"$\gamma$", fontsize = Fontsize)
            Colorbar.ax.tick_params(labelsize = Fontsize)     

        Tick_plot.set_aspect(1)
        Tick_plot.set_ylabel(r"$y\,[M]$", fontsize = Fontsize)
        Tick_plot.set_xlabel(r"$x\,[M]$", fontsize = Fontsize)
        Tick_plot.tick_params(axis='x', labelsize=Fontsize)
        Tick_plot.tick_params(axis='y', labelsize=Fontsize)
        Tick_plot.set_xlim([-8, 8])

        Schwarzschild_class  = Spacetimes.Schwarzschild()
        Schwarzschild_metric = Schwarzschild_class.metric(r = Schw_Parser.Source_R_Coord, theta = pi / 2)

        _, Schw_Scaled_Polzarization_Vectors = get_polarization_vector(Photon_Momentum = Schw_Photon_Momentum,
                                                                                               Fluid_Velocity = Fluid_Velocity,
                                                                                               B_Field = B_Field,
                                                                                               Image_Coords = Schw_Image_Coordiantes,
                                                                                               Metric = Schwarzschild_metric)
        
        Schw_Scaled_Polzarization_Vectors = Scalar * Schw_Scaled_Polzarization_Vectors

        Param_Sweep_Number = int(len(Other_Metric_Parser.X_coords) / Other_Metric_Parser.Photon_Number)

        for Param_sweep_index in range(Param_Sweep_Number):

            """ Read the parameter value from the ray-tracer log """

            Param_Sweep_Value = Other_Metric_Parser.Param_1[Param_sweep_index * Other_Metric_Parser.Photon_Number]
            Param_Sweep_Min_Value = Other_Metric_Parser.Param_1[0]
            Param_Sweep_Max_Value = Other_Metric_Parser.Param_1[-1]

            if Other_Metric_Param_Number == 2:
                Param_Sweep_Value = Other_Metric_Parser.Param_2[Param_sweep_index * Other_Metric_Parser.Photon_Number]
                Param_Sweep_Min_Value = Other_Metric_Parser.Param_2[0]
                Param_Sweep_Max_Value = Other_Metric_Parser.Param_2[-1]

            cmap = matplotlib.colormaps[colormap]

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

            Polarization_Vectors, Scaled_Polzarization_Vectors = get_polarization_vector(Photon_Momentum = Other_Metric_Photon_Momentum[Param_sweep_index * Other_Metric_Parser.Photon_Number : (Param_sweep_index + 1) * Other_Metric_Parser.Photon_Number],
                                                                                          Fluid_Velocity = Fluid_Velocity,
                                                                                          B_Field = B_Field,
                                                                                          Image_Coords = Other_Metric_Image_Coordiantes[Param_sweep_index * Other_Metric_Parser.Photon_Number : (Param_sweep_index + 1) * Other_Metric_Parser.Photon_Number],
                                                                                          Metric = Other_Metric[Param_sweep_index * Other_Metric_Parser.Photon_Number : (Param_sweep_index + 1) * Other_Metric_Parser.Photon_Number])
            
            Scaled_Polzarization_Vectors = Scalar * Scaled_Polzarization_Vectors
        
            Quiver_plot_WH = Tick_plot.quiver(Other_Metric_Parser.X_coords[Param_sweep_index * Other_Metric_Parser.Photon_Number : (Param_sweep_index + 1) * Other_Metric_Parser.Photon_Number : QUIVER_SAMPLE_SKIP] - Scaled_Polzarization_Vectors.T[0][0::QUIVER_SAMPLE_SKIP] / 2,
                                Other_Metric_Parser.Y_coords[Param_sweep_index * Other_Metric_Parser.Photon_Number : (Param_sweep_index + 1) * Other_Metric_Parser.Photon_Number : QUIVER_SAMPLE_SKIP] - Scaled_Polzarization_Vectors.T[1][0::QUIVER_SAMPLE_SKIP] / 2,
                                Scaled_Polzarization_Vectors.T[0][0::QUIVER_SAMPLE_SKIP],
                                Scaled_Polzarization_Vectors.T[1][0::QUIVER_SAMPLE_SKIP],
                                headwidth = 0,
                                headlength = 0,
                                headaxislength = 0,
                                angles = 'xy', 
                                scale_units = 'xy', 
                                scale = 1,
                                color = Color_range,
                                width = 0.01)

        """ Plot the reference values for Schwarzschild """

        Tick_plot.quiver(Schw_Parser.X_coords[0::QUIVER_SAMPLE_SKIP] - Schw_Scaled_Polzarization_Vectors.T[0][0::QUIVER_SAMPLE_SKIP] / 2,
                        Schw_Parser.Y_coords[0::QUIVER_SAMPLE_SKIP] - Schw_Scaled_Polzarization_Vectors.T[1][0::QUIVER_SAMPLE_SKIP] / 2,
                        Schw_Scaled_Polzarization_Vectors.T[0][0::QUIVER_SAMPLE_SKIP],
                        Schw_Scaled_Polzarization_Vectors.T[1][0::QUIVER_SAMPLE_SKIP],
                        headwidth = 0,
                        headlength = 0,
                        headaxislength = 0,
                        angles = 'xy', 
                        scale_units = 'xy',
                        scale=1,)
             
def plot_delta_figures(Schw_Parser, Other_Metric_Parser, B_Fields, beta_angles, Other_Metric_Param_Number = 1, QUIVER_SAMPLE_SKIP = 20, PARAWM_SWEEP_FIGURE_SKIP = 50, Fontsize = 26, Scalar = 3):

    # ====================================== Ray-Tracer Log Parser ====================================== #

    Schw_Photon_Momentum         = column_stack((Schw_Parser.Radial_Momentum, Schw_Parser.Theta_Momentum, Schw_Parser.Phi_Momentum))
    Schw_Image_Coordiantes       = column_stack((Schw_Parser.X_coords, Schw_Parser.Y_coords))

    Other_Metric_Photon_Momentum         = column_stack((Other_Metric_Parser.Radial_Momentum, Other_Metric_Parser.Theta_Momentum, Other_Metric_Parser.Phi_Momentum))
    Other_Metric_Image_Coordiantes       = column_stack((Other_Metric_Parser.X_coords, Other_Metric_Parser.Y_coords))

    # ==================================== Polarization Calculations ==================================== #

    Schwarzschild_class  = Spacetimes.Schwarzschild()
    Schwarzschild_metric = Schwarzschild_class.metric(r = Schw_Parser.Source_R_Coord, theta = pi / 2)

    Figure_Param_Sweep, (Delta_I_Plot_Sweep, Delta_EVPA_Plot_Sweep) = plt.subplots(1, 2, gridspec_kw = {'width_ratios': [1, 1]}, constrained_layout = True)
    Figure_Param_Sweep.set_figwidth(25)
    Figure_Param_Sweep.set_figheight(12)

    Delta_I_Plot_Sweep.set_ylabel(r"$\max \Delta I\,[-]$", fontsize = Fontsize)
    Delta_I_Plot_Sweep.set_xlabel(r"$\gamma\,[-]$", fontsize = Fontsize)
    Delta_I_Plot_Sweep.tick_params(axis='x', labelsize=Fontsize)
    Delta_I_Plot_Sweep.tick_params(axis='y', labelsize=Fontsize)

    Delta_EVPA_Plot_Sweep.set_ylabel(r"$\max \Delta EVPA$ [rad]", fontsize = Fontsize)
    Delta_EVPA_Plot_Sweep.set_xlabel(r"$\gamma\,[-]$", fontsize = Fontsize)
    Delta_EVPA_Plot_Sweep.tick_params(axis='x', labelsize=Fontsize)
    Delta_EVPA_Plot_Sweep.tick_params(axis='y', labelsize=Fontsize)

    colorbar_map = matplotlib.cm.ScalarMappable(cmap = "plasma")
    
    Tick_plot_list       = []
    Intensity_Plot_list  = []
    Delta_I_Plot_list    = []
    EVPA_Plot_list       = []
    Delta_EVPA_Plot_list = []

    for B_Field_Idx, (B_Field, beta_angle) in enumerate(zip(B_Fields, beta_angles)):

        Fluid_Velocity = 0.3 * array([cos(beta_angle), sin(beta_angle), 0]) # This vector has components [x, y, z]

        Figure_Pattern, (Tick_plot, Intensity_Plot, Delta_I_Plot, EVPA_Plot, Delta_EVPA_Plot) = plt.subplots(1, 5, gridspec_kw = {'width_ratios': [2.4, 1, 1, 1, 1]}, constrained_layout = True)
        Figure_Pattern.set_figwidth(40)
        Figure_Pattern.set_figheight(8)

        Tick_plot.set_aspect(1)
        Tick_plot.set_ylabel(r"$y\,[M]$", fontsize = Fontsize)
        Tick_plot.set_xlabel(r"$x\,[M]$", fontsize = Fontsize)
        Tick_plot.set_xlim([-8, 8])
        Tick_plot.set_xticks(arange(-8, 9, 2.0))
        Tick_plot.set_ylim([-8, 8])
        Tick_plot.set_yticks(arange(-8, 9, 2.0))
        Tick_plot.tick_params(axis='x', labelsize=Fontsize)
        Tick_plot.tick_params(axis='y', labelsize=Fontsize)

        Intensity_Plot.set_ylabel(r"$I$ [-]", fontsize = Fontsize)
        Intensity_Plot.set_xlabel(r"$\phi$ [$\pi$]", fontsize = Fontsize)
        Intensity_Plot.set_xlim([0, 2])
        Intensity_Plot.tick_params(axis='x', labelsize=Fontsize)
        Intensity_Plot.tick_params(axis='y', labelsize=Fontsize)

        Delta_I_Plot.set_ylabel(r"$\Delta I$ [-]", fontsize = Fontsize)
        Delta_I_Plot.set_xlabel(r"$\phi$ [$\pi$]", fontsize = Fontsize)
        Delta_I_Plot.set_xlim([0, 2])
        Delta_I_Plot.tick_params(axis='x', labelsize=Fontsize)
        Delta_I_Plot.tick_params(axis='y', labelsize=Fontsize)

        EVPA_Plot.set_ylabel(r"$EVPA$ [rad]", fontsize = Fontsize)
        EVPA_Plot.set_xlabel(r"$\phi$ [$\pi$]", fontsize = Fontsize)
        EVPA_Plot.set_ylim([-pi / 2, pi / 2])
        EVPA_Plot.set_xlim([0, 2])
        EVPA_Plot.tick_params(axis='x', labelsize=Fontsize)
        EVPA_Plot.tick_params(axis='y', labelsize=Fontsize)

        Delta_EVPA_Plot.set_ylabel(r"$\Delta EVPA$ [rad]", fontsize = Fontsize)
        Delta_EVPA_Plot.set_xlabel(r"$\phi$ [$\pi$]", fontsize = Fontsize)
        Delta_EVPA_Plot.set_xlim([0, 2])
        Delta_EVPA_Plot.tick_params(axis='x', labelsize=Fontsize)
        Delta_EVPA_Plot.tick_params(axis='y', labelsize=Fontsize)

        Colorbar = Figure_Param_Sweep.colorbar(colorbar_map, ax = Delta_EVPA_Plot)
        Colorbar.set_label(r"$\gamma$", fontsize = Fontsize)
        Colorbar.ax.tick_params(labelsize = Fontsize)     

        Schw_Polarization_Vectors, Schw_Scaled_Polzarization_Vectors = get_polarization_vector(Photon_Momentum = Schw_Photon_Momentum,
                                                                                              Fluid_Velocity = Fluid_Velocity,
                                                                                              B_Field = B_Field,
                                                                                              Image_Coords = Schw_Image_Coordiantes,
                                                                                              Metric = Schwarzschild_metric)

        Schw_Scaled_Polzarization_Vectors = Scalar * Schw_Scaled_Polzarization_Vectors

        Schw_EVPA, Schw_EVPA_Branches, Schw_Polarization_I, Schw_Polarization_Q, Schw_Polarization_U, Schw_Image_phi_coord = get_observational_quantities(Polarization_Vectors = Schw_Polarization_Vectors,
                                                                                                                                                        Image_Coordinates    = Schw_Image_Coordiantes)
        
        """ Plot the reference values for Schwarzschild """

        Tick_plot.quiver(Schw_Parser.X_coords[0::QUIVER_SAMPLE_SKIP] - Schw_Scaled_Polzarization_Vectors.T[0][0::QUIVER_SAMPLE_SKIP] / 2,
                        Schw_Parser.Y_coords[0::QUIVER_SAMPLE_SKIP] - Schw_Scaled_Polzarization_Vectors.T[1][0::QUIVER_SAMPLE_SKIP] / 2,
                        Schw_Scaled_Polzarization_Vectors.T[0][0::QUIVER_SAMPLE_SKIP],
                        Schw_Scaled_Polzarization_Vectors.T[1][0::QUIVER_SAMPLE_SKIP],
                        headwidth = 0,
                        headlength = 0,
                        headaxislength = 0,
                        angles = 'xy', 
                        scale_units = 'xy', 
                        scale = 1,
                        color = 'k')

        Intensity_Plot.plot(Schw_Image_phi_coord, Schw_Polarization_I, "-", color = "k")

        for index in range(len(Schw_EVPA_Branches) - 1):

            EVPA_Plot.plot(Schw_Image_phi_coord[Schw_EVPA_Branches[index] : Schw_EVPA_Branches[index + 1] - 1], Schw_EVPA[Schw_EVPA_Branches[index] : Schw_EVPA_Branches[index + 1] - 1], "--", color = 'k')
            
        """ Big fat loop over the parameter sweep of the other metric """

        Max_Delta_I    = []
        Max_Delta_EVPA = []
        Param_Sweep_Values = []

        Param_Sweep_Number = int(len(Other_Metric_Parser.X_coords) / Other_Metric_Parser.Photon_Number)
        Param_Sweep_Color_Cycle = ["r", "b", "k"]

        for Param_sweep_index in range(Param_Sweep_Number):

            """ Read the parameter value from the ray-tracer log """

            Param_Sweep_Value = Other_Metric_Parser.Param_1[Param_sweep_index * Other_Metric_Parser.Photon_Number]
            Param_Sweep_Min_Value = Other_Metric_Parser.Param_1[0]
            Param_Sweep_Max_Value = Other_Metric_Parser.Param_1[-1]

            if Other_Metric_Param_Number == 2:
                Param_Sweep_Value = Other_Metric_Parser.Param_2[Param_sweep_index * Other_Metric_Parser.Photon_Number]
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

            Polarization_Vectors, Scaled_Polzarization_Vectors = get_polarization_vector(Photon_Momentum = Other_Metric_Photon_Momentum[Param_sweep_index * Other_Metric_Parser.Photon_Number : (Param_sweep_index + 1) * Other_Metric_Parser.Photon_Number],
                                                                                        Fluid_Velocity = Fluid_Velocity,
                                                                                        B_Field = B_Field,
                                                                                        Image_Coords = Other_Metric_Image_Coordiantes[Param_sweep_index * Other_Metric_Parser.Photon_Number : (Param_sweep_index + 1) * Other_Metric_Parser.Photon_Number],
                                                                                        Metric = Other_Metric[Param_sweep_index * Other_Metric_Parser.Photon_Number : (Param_sweep_index + 1) * Other_Metric_Parser.Photon_Number])
        
            Scaled_Polzarization_Vectors = Scalar * Scaled_Polzarization_Vectors

            EVPA, EVPA_Branches, Polarization_I, Polarization_Q, Polarization_U, Image_phi_coord = get_observational_quantities(Polarization_Vectors = Polarization_Vectors,
                                                                                                                                Image_Coordinates = Other_Metric_Image_Coordiantes[Param_sweep_index * Other_Metric_Parser.Photon_Number : (Param_sweep_index + 1) * Other_Metric_Parser.Photon_Number])

            Delta_I = Polarization_I - Schw_Polarization_I
            Delta_EVPA = EVPA - Schw_EVPA

            for index, _ in enumerate(Delta_EVPA):

                if Delta_EVPA[index] > pi / 2:
                    Delta_EVPA[index] = Delta_EVPA[index] - pi

                elif Delta_EVPA[index] < -pi / 2:
                    Delta_EVPA[index] = Delta_EVPA[index] + pi

            if Param_sweep_index % PARAWM_SWEEP_FIGURE_SKIP == 0:

                Tick_plot.quiver(Other_Metric_Parser.X_coords[Param_sweep_index * Other_Metric_Parser.Photon_Number : (Param_sweep_index + 1) * Other_Metric_Parser.Photon_Number : QUIVER_SAMPLE_SKIP] - Scaled_Polzarization_Vectors.T[0][0::QUIVER_SAMPLE_SKIP] / 2,
                                Other_Metric_Parser.Y_coords[Param_sweep_index * Other_Metric_Parser.Photon_Number : (Param_sweep_index + 1) * Other_Metric_Parser.Photon_Number : QUIVER_SAMPLE_SKIP] - Scaled_Polzarization_Vectors.T[1][0::QUIVER_SAMPLE_SKIP] / 2,
                                Scaled_Polzarization_Vectors.T[0][0::QUIVER_SAMPLE_SKIP],
                                Scaled_Polzarization_Vectors.T[1][0::QUIVER_SAMPLE_SKIP],
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

                Delta_I_Plot.plot(Image_phi_coord, Delta_I, color = Color_range)

                Delta_EVPA_Plot.plot(Image_phi_coord, Delta_EVPA, color = Color_range)

            Max_Delta_I_idx    = argmax(abs(Delta_I))  
            Max_Delta_EVPA_idx = argmax(abs(Delta_EVPA))      

            Max_Delta_I.append(Delta_I[Max_Delta_I_idx])
            Max_Delta_EVPA.append(Delta_EVPA[Max_Delta_EVPA_idx])
            Param_Sweep_Values.append(Param_Sweep_Value)

        """ ================ Append the plots to lists, so I can scale all their axis later, after I know all the limits ================ """

        Tick_plot_list.append(Tick_plot)
        Intensity_Plot_list.append(Intensity_Plot)
        Delta_I_Plot_list.append(Delta_I_Plot)
        EVPA_Plot_list.append(EVPA_Plot)
        Delta_EVPA_Plot_list.append(Delta_EVPA_Plot)

        """ ==================================== """

        Max_Delta_EVPA_Branches = split_Deltas(Delta = Max_Delta_EVPA)
        Max_Delta_I_Branches    = split_Deltas(Delta = Max_Delta_I) 

        for index in range(len(Max_Delta_EVPA_Branches) - 1):

            if index == 0:
                Delta_EVPA_Plot_Sweep.plot(Param_Sweep_Values[Max_Delta_EVPA_Branches[index] : Max_Delta_EVPA_Branches[index + 1] - 1], 
                                        Max_Delta_EVPA[Max_Delta_EVPA_Branches[index] : Max_Delta_EVPA_Branches[index + 1] - 1], 
                                        color = Param_Sweep_Color_Cycle[B_Field_Idx], 
                                        label = "B = " + 
                                                "[{}, {}, {}], ".format(B_Field[0], B_Field[1], B_Field[2]) + 
                                                r"$\chi$" + 
                                                "= {}".format(round(beta_angle * 180 / pi)) + 
                                                r"$^\circ$")
            else:
                Delta_EVPA_Plot_Sweep.plot(Param_Sweep_Values[Max_Delta_EVPA_Branches[index] : Max_Delta_EVPA_Branches[index + 1] - 1], 
                                        Max_Delta_EVPA[Max_Delta_EVPA_Branches[index] : Max_Delta_EVPA_Branches[index + 1] - 1], 
                                        color = Param_Sweep_Color_Cycle[B_Field_Idx])
                
        for index in range(len(Max_Delta_I_Branches) - 1):

            if index == 0:
                Delta_I_Plot_Sweep.plot(Param_Sweep_Values[Max_Delta_I_Branches[index] : Max_Delta_I_Branches[index + 1] - 1], 
                                        Max_Delta_I[Max_Delta_I_Branches[index] : Max_Delta_I_Branches[index + 1] - 1], 
                                        color = Param_Sweep_Color_Cycle[B_Field_Idx], 
                                        label = "B = " + 
                                                "[{}, {}, {}], ".format(B_Field[0], B_Field[1], B_Field[2]) + 
                                                r"$\chi$" + 
                                                "= {}".format(round(beta_angle * 180 / pi)) + 
                                                r"$^\circ$")
                                                    
            else:
                Delta_I_Plot_Sweep.plot(Param_Sweep_Values[Max_Delta_I_Branches[index] : Max_Delta_I_Branches[index + 1] - 1], 
                                        Max_Delta_I[Max_Delta_I_Branches[index] : Max_Delta_I_Branches[index + 1] - 1], 
                                        color = Param_Sweep_Color_Cycle[B_Field_Idx])
            
        Delta_EVPA_Plot_Sweep.set_xlim([Param_Sweep_Min_Value, Param_Sweep_Max_Value])
        Delta_I_Plot_Sweep.set_xlim([Param_Sweep_Min_Value, Param_Sweep_Max_Value])

        colorbar_map.set_clim([Param_Sweep_Min_Value, Param_Sweep_Max_Value])

    Delta_I_Plot_Sweep.set_xlim([Param_Sweep_Min_Value, Param_Sweep_Max_Value])
    Delta_EVPA_Plot_Sweep.set_xlim([Param_Sweep_Min_Value, Param_Sweep_Max_Value])

    Delta_EVPA_Plot_Sweep.legend(loc = "upper right", fontsize = 26)
    Delta_I_Plot_Sweep.legend(loc = "lower left", fontsize = 26)

    Delta_I_Plot_Sweep.plot([min(Param_Sweep_Values), max(Param_Sweep_Values)], [0, 0], "--", color = "k")
    Delta_EVPA_Plot_Sweep.plot([min(Param_Sweep_Values), max(Param_Sweep_Values)], [0, 0], "--", color = "k")

    ratio = 1.5

    """ ================================= Scale all Y limits equally ================================= """

    I_y_limits = [0, 0]
    Delta_I_y_limits = [0, 0]
    Delta_EVPA_y_limits = [0, 0]

    for index, _ in enumerate(Intensity_Plot_list):

        I_y_min, I_y_max = Intensity_Plot_list[index].get_ylim()
        Delta_I_y_min, Delta_I_y_max = Delta_I_Plot_list[index].get_ylim()
        Delta_EVPA_min, Delta_EVPA_max = Delta_EVPA_Plot_list[index].get_ylim()

        EVPA_Plot_list[index].set_ylim([-pi / 2, pi / 2])        
        EVPA_Plot_list[index].set_aspect(2 / pi * ratio)

        if index > 0:

            if I_y_min < I_y_limits[0]:
                I_y_limits[0] = I_y_min
            
            if I_y_max > I_y_limits[1]:
                I_y_limits[1] = I_y_max

            if Delta_I_y_min < Delta_I_y_limits[0]:
                Delta_I_y_limits[0] = Delta_I_y_min
            
            if Delta_I_y_max > Delta_I_y_limits[1]:
                Delta_I_y_limits[1] = Delta_I_y_max

            if Delta_EVPA_min < Delta_EVPA_y_limits[0]:
                Delta_EVPA_y_limits[0] = Delta_EVPA_min
            
            if Delta_EVPA_max > Delta_EVPA_y_limits[1]:
                Delta_EVPA_y_limits[1] = Delta_EVPA_max
        
        else:

            I_y_limits = [I_y_min, I_y_max] 
            Delta_I_y_limits = [Delta_I_y_min, Delta_I_y_max]
            Delta_EVPA_y_limits = [Delta_EVPA_min, Delta_EVPA_max]

    I_y_limits[0] = (I_y_limits[0] * 0.9) * (I_y_limits[0] > 0) + (I_y_limits[0] * 1.1) * (I_y_limits[0] < 0)
    I_y_limits[1] = (I_y_limits[1] * 0.9) * (I_y_limits[1] < 0) + (I_y_limits[1] * 1.1) * (I_y_limits[1] > 0)

    Delta_I_y_limits[0] = (Delta_I_y_limits[0] * 0.9) * (Delta_I_y_limits[0] > 0) + (Delta_I_y_limits[0] * 1.1) * (Delta_I_y_limits[0] < 0)
    Delta_I_y_limits[1] = (Delta_I_y_limits[1] * 0.9) * (Delta_I_y_limits[1] < 0) + (Delta_I_y_limits[1] * 1.1) * (Delta_I_y_limits[1] > 0)

    Delta_EVPA_y_limits[0] = (Delta_EVPA_y_limits[0] * 0.9) * (Delta_EVPA_y_limits[0] > 0) + (Delta_EVPA_y_limits[0] * 1.1) * (Delta_EVPA_y_limits[0] < 0)
    Delta_EVPA_y_limits[1] = (Delta_EVPA_y_limits[1] * 0.9) * (Delta_EVPA_y_limits[1] < 0) + (Delta_EVPA_y_limits[1] * 1.1) * (Delta_EVPA_y_limits[1] > 0)

    for index, _ in enumerate(Intensity_Plot_list):

        Intensity_Plot_list[index].set_ylim([I_y_limits[0], I_y_limits[1]])
        Intensity_Plot_list[index].set_aspect(abs(2 / (I_y_limits[0] - I_y_limits[1])) * ratio)

        Delta_I_Plot_list[index].set_ylim([Delta_I_y_limits[0], Delta_I_y_limits[1]])
        Delta_I_Plot_list[index].set_aspect(abs(2 / (Delta_I_y_limits[0] - Delta_I_y_limits[1])) * ratio)

        Delta_EVPA_Plot_list[index].set_ylim([Delta_EVPA_y_limits[0], Delta_EVPA_y_limits[1]])
        Delta_EVPA_Plot_list[index].set_aspect(abs(2 / (Delta_EVPA_y_limits[0] - Delta_EVPA_y_limits[1])) * ratio)

if __name__ == '__main__':
 
    params = {"ytick.color" : "black",
              "xtick.color" : "black",
              "axes.labelcolor" : "black",
              "axes.edgecolor" : "black",
              "text.usetex" : True,
              "font.family" : "serif",
              "font.serif" : ["Computer Modern Serif"]}
    plt.rcParams.update(params)

    # ======================================================= Setup the Figure ======================================================= #

    # Figure_Pattern, (Tick_plot_1, Tick_plot_2, Tick_plot_3) = plt.subplots(1, 3, gridspec_kw = {'width_ratios': [1, 1, 1]}, constrained_layout = True)
    # Figure_Pattern.set_figwidth(12.5)
    # Figure_Pattern.set_figheight(25)

    # colorbar_map = matplotlib.cm.ScalarMappable(cmap = matplotlib.colormaps['plasma'])
    # colorbar_map.set_clim([0,3])
            
    # Colorbar = Figure_Pattern.colorbar(colorbar_map, ax = Tick_plot_3)
    # Colorbar.set_label(r"$\gamma$", fontsize = 32)
    # Colorbar.ax.tick_params(labelsize = 26)    

    # Tick_plot_2.axes.get_yaxis().set_visible(False)
    # Tick_plot_3.axes.get_yaxis().set_visible(False)

    # Tick_plot_1.axes.get_xaxis().set_ticks(arange(-8, 10, 2.0))
    # Tick_plot_1.axes.get_yaxis().set_ticks(arange(-8, 10, 2.0))
    # Tick_plot_2.axes.get_xaxis().set_ticks(arange(-8, 10, 2.0))
    # Tick_plot_3.axes.get_xaxis().set_ticks(arange(-8, 10, 2.0))

    # ================================================================================================================================= #

    # Schw_Sim_Path = "C:\\Users\\Valur\\Documents\\Repos\\Gravitational_Lenser\\Sim_Results\\Schwarzschild_n0_r4.5_20_deg"
    # Schw_Parser_4_5   = Simulation_Parser(Schw_Sim_Path)

    # Other_Metric_Sim_Path = "C:\\Users\\Valur\\Documents\\Repos\\Gravitational_Lenser\\Sim_Results\\Wormhole_n0_r4.5_gamma_scan_20_deg"
    # Other_Metric_Parser_4_5   = Simulation_Parser(Other_Metric_Sim_Path)

    Schw_Sim_Path = "C:\\Users\\Valur\\Documents\\Repos\\Gravitational_Lenser\\Sim_Results\\Kerr_n0"
    Schw_Parser_6 = Simulation_Parser(Schw_Sim_Path)

    Other_Metric_Sim_Path = "C:\\Users\\Valur\\Documents\\Repos\\Gravitational_Lenser\\Sim_Results\\JNW_n0"
    Other_Metric_Parser_6   = Simulation_Parser(Other_Metric_Sim_Path)

    QUIVER_SAMPLE_SKIP = 8
    PARAWM_SWEEP_FIGURE_SKIP = 20

    """  
    For metrics with more than one parameter these go as:
        * Wormhole: 1 = Spin, 2 = Redshift
        * Black Hole With Dark Matter Halo: 1 - Halo Mass, 2 = Halo Compactness
    """

    Other_Metric_Param_Number = 1

    B_Fields = [array([0.87, 0, 0.5]),
                array([0.71, 0, 0.71]),
                array([0.5, 0, 0.87])] # This vector has components [r, theta, phi]

    beta_angles = [-150. / 180 * pi,
                   -135. / 180 * pi,
                   -120. / 180 * pi]

    plot_delta_figures(Schw_Parser = Schw_Parser_6, 
                       Other_Metric_Parser = Other_Metric_Parser_6,
                       B_Fields = B_Fields,
                       beta_angles = beta_angles,
                       Other_Metric_Param_Number = Other_Metric_Param_Number,
                       Fontsize = 32,
                       Scalar = 1.5,
                       PARAWM_SWEEP_FIGURE_SKIP = PARAWM_SWEEP_FIGURE_SKIP, 
                       QUIVER_SAMPLE_SKIP = QUIVER_SAMPLE_SKIP)

    # ======================================= Initial Conditions ======================================= #

    # B_Fields = [array([0.5, 0, 0.87])] # This vector has components [r, theta, phi]

    # beta_angles = [-120 / 180 * pi]

    # Tick_plot_1.set_title(r'B = [0.5, 0.87, 0] $\beta = 0.3,\,\chi = -120^\circ$', fontsize = 36)

    # plot_polarization_ticks(Schw_Parser = Schw_Parser_6, 
    #                         Other_Metric_Parser = Other_Metric_Parser_6, 
    #                         B_Fields = B_Fields, 
    #                         beta_angles = beta_angles, 
    #                         Other_Metric_Param_Number = 2, 
    #                         QUIVER_SAMPLE_SKIP = 18, 
    #                         Scalar = 1, 
    #                         Tick_plot = Tick_plot_1,
    #                         Fontsize = 36)
    # plot_polarization_ticks(Schw_Parser = Schw_Parser_4_5, 
    #                         Other_Metric_Parser = Other_Metric_Parser_4_5, 
    #                         B_Fields = B_Fields, beta_angles = beta_angles, 
    #                         Other_Metric_Param_Number = 2, 
    #                         QUIVER_SAMPLE_SKIP = 18, 
    #                         Scalar = 3, 
    #                         Tick_plot = Tick_plot_1,
    #                         Fontsize = 36)
    
    #  # ======================================= Initial Conditions ======================================= #

    # B_Fields = [array([0.71, 0, 0.71])] # This vector has components [r, theta, phi]

    # beta_angles = [-135 / 180 * pi]

    # Tick_plot_2.set_title(r'B = [0.71, 0.71, 0] $\beta = 0.3,\,\chi = -135^\circ$', fontsize = 36)

    # plot_polarization_ticks(Schw_Parser = Schw_Parser_6, 
    #                         Other_Metric_Parser = Other_Metric_Parser_6, 
    #                         B_Fields = B_Fields, 
    #                         beta_angles = beta_angles, 
    #                         Other_Metric_Param_Number = 2, 
    #                         QUIVER_SAMPLE_SKIP = 18, 
    #                         Scalar = 1, 
    #                         Tick_plot = Tick_plot_2,
    #                         Fontsize = 36)
    # plot_polarization_ticks(Schw_Parser = Schw_Parser_4_5, 
    #                         Other_Metric_Parser = Other_Metric_Parser_4_5, 
    #                         B_Fields = B_Fields, beta_angles = beta_angles, 
    #                         Other_Metric_Param_Number = 2, 
    #                         QUIVER_SAMPLE_SKIP = 18, 
    #                         Scalar = 3, 
    #                         Tick_plot = Tick_plot_2,
    #                         Fontsize = 36)
    

    #  # ======================================= Initial Conditions ======================================= #

    # B_Fields = [array([0.87, 0, 0.5])] # This vector has components [r, theta, phi]

    # beta_angles = [-150 / 180 * pi]

    # Tick_plot_3.set_title(r'B = [0.87, 0.5, 0] $\beta = 0.3,\,\chi = -150^\circ$', fontsize = 36)

    # plot_polarization_ticks(Schw_Parser = Schw_Parser_6, 
    #                         Other_Metric_Parser = Other_Metric_Parser_6, 
    #                         B_Fields = B_Fields, 
    #                         beta_angles = beta_angles, 
    #                         Other_Metric_Param_Number = 2, 
    #                         QUIVER_SAMPLE_SKIP = 18, 
    #                         Scalar = 1, 
    #                         Tick_plot = Tick_plot_3,
    #                         Fontsize = 36)
    # plot_polarization_ticks(Schw_Parser = Schw_Parser_4_5, 
    #                         Other_Metric_Parser = Other_Metric_Parser_4_5, 
    #                         B_Fields = B_Fields, beta_angles = beta_angles, 
    #                         Other_Metric_Param_Number = 2, 
    #                         QUIVER_SAMPLE_SKIP = 18, 
    #                         Scalar = 3, 
    #                         Tick_plot = Tick_plot_3,
    #                         Fontsize = 36)
    
    plt.show()
