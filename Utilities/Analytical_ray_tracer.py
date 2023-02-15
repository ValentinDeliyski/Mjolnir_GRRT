import warnings
import numpy as np

from scipy.stats import beta
from scipy.integrate import quad
from scipy.interpolate import CubicSpline

import matplotlib.pyplot as plt

from Support_functions.Spacetimes import *

def Integrand(r, spacetime, impact_param):

    metric = spacetime.metric(r, np.pi / 2)
    V_eff = 1 + metric[0] / metric[2] * impact_param**2 

    if impact_param == 0:
        return 0

    function = -1 / impact_param**2 * metric[3]**2 / metric[0] / metric[1] * V_eff

    if function < 0 or function != function or function == 0:
        return 0

    return 1 / np.sqrt(function)
    
def propagate_rays(Spacetime, GRANULARITY, r_source, r_obs):

    warnings.filterwarnings('ignore')

    metric_source       = Spacetime.metric(r_source, np.pi / 2)
    impact_param_source = np.sqrt(-metric_source[3] / metric_source[0])

    if Spacetime.HAS_PHOTON_SPHERE:

        metric_photon_sphere       = Spacetime.metric(Spacetime.photon_sphere(), np.pi / 2)
        impact_param_photon_sphere = np.sqrt(-metric_photon_sphere[3] / metric_photon_sphere[0])

    print("Starting to compute direct photon paths...")

    # Calculate the impact parameters for photons that do not encounter a radian turning point

    # NOTE: To get a smoother image at the end it is better to distribute the impact paramters with a beta distrib.
    # A density parameter of 2 for both arguments seems like a good distribution - weighed more at the ends of the interval

    density_parameter  = 2
    distribution_range = np.linspace(0, 1, GRANULARITY)

    if r_source < 0:
        direct_impact_params = (beta.cdf(distribution_range, density_parameter, density_parameter) * impact_param_photon_sphere).tolist()

    else:
        direct_impact_params = (beta.cdf(distribution_range, density_parameter, density_parameter) * impact_param_source).tolist()

    results = []

    for impact_param in direct_impact_params:

        azimuthal_distance = quad(Integrand, 
                                  r_source, # Lower integration limit
                                  r_obs,    # Upper integration limit
                                  args  = (Spacetime, impact_param), 
                                  limit = 10000,
                                  epsabs = 1e-15,
                                  epsrel = 1e-15)[0]

        results.append(azimuthal_distance)

    print("Finished computing direct photon paths!")

    if r_source > 0:

        print("Starting to compute higher order photon paths...")

        if Spacetime.HAS_PHOTON_SPHERE:

            # In the presence of a photon sphere, the turning points will be between it and the source

            # higher_order_turning_points = np.linspace(r_source, Spacetime.photon_sphere(), GRANULARITY)
            density_parameter  = 5.5
            higher_order_turning_points = r_source + (beta.cdf(distribution_range, density_parameter, density_parameter)) * (Spacetime.photon_sphere() - r_source)

            metrics_at_turning_points   = Spacetime.metric(higher_order_turning_points, np.pi / 2)

            # Calculate the impact parameters for the correspoinding turning points

            higher_order_impact_params  = np.sqrt(-metrics_at_turning_points[3] / metrics_at_turning_points[0])

        else:

            # NOTE: Currently only applies to the Janis-Newman-Winicour Naked Singularities
            # In the absence of a photon sphere the turning points are located between the singularity and the source

            density_parameter  = 105.5
            distribution_range = np.linspace(0, 1, GRANULARITY)
            r_singularity      = 2 * Spacetime.MASS / Spacetime.PARAMETER

            # NOTE: To properly construct the images we need a very uneven distribution of photon turning points
            # The majority of the points must be distributed at the ends of the interval - this is neatly done with a beta distribution

            higher_order_turning_points = r_source + (beta.cdf(distribution_range, density_parameter, density_parameter)) * (r_singularity - r_source)

            # higher_order_turning_points = np.linspace(r_source, r_singularity)

            # Calculate the impact parameters for the correspoinding turning points

            higher_order_impact_params = (higher_order_turning_points * pow(1 - 2 * Spacetime.MASS / Spacetime.PARAMETER / higher_order_turning_points, 1/2 - Spacetime.PARAMETER)).tolist()

        for impact_param, turning_point in zip(higher_order_impact_params, higher_order_turning_points):

            branch_1_integral = quad(Integrand, 
                                     turning_point, # Lower integration limit
                                     r_obs,         # Upper integration limit
                                     args  = (Spacetime, impact_param), 
                                     limit = 50000,
                                     epsabs = 1e-20,
                                     epsrel = 1e-20)[0]

            branch_2_integral = quad(Integrand, 
                                     turning_point, # Lower integration limit
                                     r_source,      # Upper integration limit
                                     args  = (Spacetime, impact_param), 
                                     limit = 50000,
                                     epsabs = 1e-20,
                                     epsrel = 1e-20)[0]

            results.append(branch_1_integral + branch_2_integral)
        
        impact_parameters = direct_impact_params
        impact_parameters.extend(higher_order_impact_params)

    else:
        impact_parameters = direct_impact_params

    return np.array(impact_parameters), np.array(results)

def construct_images(Impact_parameters, Azimuthal_distance, Inclination):

    np.seterr(all = "ignore")

    Image_order = 0
    Image_angles    = []
    Image_distances = []
  
    Solution_condition = np.arcsin(1 / np.tan(Azimuthal_distance) / np.tan(Inclination)).tolist()

    start_index = 0
    end_index   = 0

    solution_starts = 0
    solution_ends   = 0

    for index, _ in enumerate(Solution_condition):

        if index != 0:
            solution_starts =     np.isnan(Solution_condition[index - 1]) and not np.isnan(Solution_condition[index])
            solution_ends   = not np.isnan(Solution_condition[index - 1]) and     np.isnan(Solution_condition[index])

            if index == len(Solution_condition) - 1 and not np.isnan(Solution_condition[index]):
                solution_ends = True

        if solution_starts:
            start_index = index

        if solution_ends:
            end_index = index

            Image_angles.append(Solution_condition[ start_index:end_index ])
            Image_distances.append(Impact_parameters[ start_index:end_index ])

            if len(Image_angles[Image_order]) > 20:

                # Image_angles[Image_order], Image_distances[Image_order] = zip(* sorted(zip(Image_angles[Image_order], Image_distances[Image_order])) )

                # Image_angles[Image_order] = Image_angles[Image_order] + np.linspace(0, 1e-10, len(Image_angles[Image_order]))

                # Spline = CubicSpline(Image_angles[Image_order], Image_distances[Image_order])

                # Angles_spline = np.linspace(-np.pi / 2, np.pi / 2, 1000)
                # Distances_spline = Spline(Angles_spline)

                x_coords =  np.array(Image_distances[Image_order]) * np.cos(Image_angles[Image_order])
                y_coords = -np.array(Image_distances[Image_order]) * np.sin(Image_angles[Image_order])

                # x_coords =  np.array(Distances_spline) * np.cos(Angles_spline)
                # y_coords = -np.array(Distances_spline) * np.sin(Angles_spline)

                x_coords = np.append( np.append(-x_coords, np.flip(x_coords)), -x_coords[0] ) 
                y_coords = np.append( np.append( y_coords, np.flip(y_coords)),  y_coords[0] )

                plt.plot(x_coords, y_coords)

                Image_order = Image_order + 1

    print("Number of images constructed: = ", Image_order)

    plt.show()

#------    Constants      -------#

DEG_TO_RAD = np.pi / 180

#------ Metric Parameters -------#

WH_THROAT = 1    # [ M ]
WH_ALPHA  = 2    # [ - ]
     
RBH_PARAM = 1    # [ - ]

JNW_PARAM = 0.48  # [ - ]

#------      Metrics      -------#

SCH = Schwarzschild()
WH  = Wormhole(r_throat = WH_THROAT, parameter = WH_ALPHA)
RBH = Regular_Black_Hole(parameter = RBH_PARAM)
JNW = JNW_Naked_Singularity(parameter = JNW_PARAM)

Spacetime_dict ={"Schwarzshild":       SCH, 
                 "Wormhole":            WH,
                 "Regular Black Hole": RBH,
                 "Naked Singularity":  JNW}

#----- Observer / Source  -------#

r_obs = 1e3                         # [ M ]
inclination_obs = 80 * DEG_TO_RAD   # [ rad ]

Impact_Parameters, Azimuthal_distance = propagate_rays(Spacetime = Spacetime_dict["Naked Singularity"], 
                                                       GRANULARITY = 5000, 
                                                       r_source = 25, 
                                                       r_obs = r_obs)

construct_images(Impact_parameters = Impact_Parameters, Azimuthal_distance = Azimuthal_distance, Inclination = inclination_obs)