from Support_functions.Parsers import Simulation_Parser
from Support_functions.Spacetimes_new import Kerr, Coords
import matplotlib.pyplot as plt
from matplotlib import colors

from numpy import pi, array, cos, sin, arccos, arctan, sqrt, linspace, zeros

if __name__ == "__main__":
    
    plt.rcParams['axes.titlepad'] = 15
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = 'serif' # or 'sans-serif' or 'monospace'
    plt.rcParams['font.serif'] = 'cmr10'
    plt.rcParams['font.sans-serif'] = 'cmss10'
    plt.rcParams['font.monospace'] = 'cmtt10'
    plt.rcParams["axes.formatter.use_mathtext"] = True
    plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}' + '\n' + r'\usepackage{xcolor}'

    Sim_parser_Kerr = Simulation_Parser("C:/Users/Valur/Documents/Repos/Mjolnir_GRRT/Utilities/Numerical_metric_runs/Zero_curvature/Config_V_70_deg_B_0.0_1.0_0.0_min_order_0_max_order_0_zoom/Kerr_analog_single_geodesic/Kerr_photon_log")
    Position_tuple_Kerr, Momentum_tuple, _, _, Affine_param_Kerr, _, Polarization_tuple_Kerr, PW_Constant_tuple = Sim_parser_Kerr.get_photon_log()
  
    # ==================================================================  Polarization plotting ===================================================================
    
    f_x, f_y = Polarization_tuple_Kerr
    p_t, p_r, p_theta, p_phi = Momentum_tuple
    t, r, theta, phi = Position_tuple_Kerr
    
    M = float(Sim_parser_Kerr.Simulation_metadata["Mass [M]"])
    a = float(Sim_parser_Kerr.Simulation_metadata["Spin Parameter [M]"])
    
    Kerr_class_instance = Kerr(mass = M, spin_param = a)
    
    PW_Const_1, PW_Const_2 = PW_Constant_tuple
    PW_Const_1 = array(PW_Const_1)
    PW_Const_2 = array(PW_Const_2)
    
    plt.plot(PW_Const_1[PW_Const_1 != 0])
    
    plt.show()