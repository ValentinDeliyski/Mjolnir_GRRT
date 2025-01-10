from Toy_Model_Polarization import plot_delta_figures, plot_polarization_ticks, plot_UV_diagram
from Support_functions.Parsers import Simulation_Parser

from numpy import array, sqrt
from numpy import pi

import matplotlib.pyplot as plt

if __name__ == '__main__':
 
    params = {"ytick.color" : "black",
              "xtick.color" : "black",
              "axes.labelcolor" : "black",
              "axes.edgecolor" : "black",
              "text.usetex" : True,
              "font.family" : "serif",
              "font.serif" : ["Computer Modern Serif"]}
    
    plt.rcParams.update(params)

    """  
    For metrics with more than one parameter these go as:
        * Wormhole: 1 = Spin, 2 = Redshift
        * Black Hole With Dark Matter Halo: 1 - Halo Mass, 2 = Halo Compactness
    """

    Schw_Sim_Path = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\Sim_Results\\Test_Simulation\\Kerr_n0"
    Schw_Parser_6 = Simulation_Parser(Schw_Sim_Path)

    B_Fields = [array([0, 1, 0])] # This vector has components [r, theta, phi]

    beta_angles = [-pi / 2]
    
    beta_norms = [0.5]

    plot_UV_diagram(Sim_Parser = Schw_Parser_6,
                    B_Fields = B_Fields,
                    beta_angles = beta_angles,
                    beta_norms = beta_norms)

    
    plt.show()
