import csv
import numpy as np
from itertools import islice
from matplotlib import pyplot as plt

class Simulation_Parser():

    def __init__(self, File_name: str) -> None:

        with open("Sim_Results\\" + File_name + ".txt", 'r') as file:

            self.HEADER_ROW_COUNT = 5

            csvreader = csv.reader(file, delimiter = ":")
            row_count = sum(1 for row in file) - self.HEADER_ROW_COUNT
            file.seek(0)

            self.OBS_DISTANCE    = float(csvreader.__next__()[1])
            self.OBS_INCLICATION = float(csvreader.__next__()[1])

            self.WINDOW_LIMITS = [float(limit) for limit in csvreader.__next__()[1].split(',')]
            
            Resolution_list = csvreader.__next__()[1].split(' ')

            self.X_PIXEL_COUNT   = int(Resolution_list[1])
            self.Y_PIXEL_COUNT   = int(Resolution_list[3])

            self.Legend = csvreader.__next__()

            csvreader = csv.reader(file, delimiter = " ")

            self.X_coords        = np.zeros(row_count)
            self.Y_coords        = np.zeros(row_count)
            self.NT_Flux         = np.zeros(row_count)
            self.Intensity       = np.zeros(row_count)
            self.NT_Redshift     = np.zeros(row_count)
            self.NT_Flux_Shifted = np.zeros(row_count)

            index = 0

            for row in csvreader:

                self.X_coords[index] = row[0]
                self.Y_coords[index] = row[1]

                self.NT_Redshift[index]  = row[2]
                self.NT_Flux[index]      = row[3]
                self.Intensity[index]    = row[4]

                self.NT_Flux_Shifted[index] = self.NT_Redshift[index]**4*self.NT_Flux[index]

                index += 1

    def get_plottable_sim_data(self) -> tuple:

        # Arrays need to be flipped, because mpl treats y = 0 as the top, 
        # and the simulator (aka openGL) treats it as the bottom

        Intensity = self.Intensity.reshape(self.X_PIXEL_COUNT, self.Y_PIXEL_COUNT)
        Intensity = np.flip(Intensity, 0)

        NT_Flux         = self.NT_Flux.reshape(self.X_PIXEL_COUNT,self.Y_PIXEL_COUNT)
        NT_Flux         = np.flip(NT_Flux, 0)

        NT_Redshift     = self.NT_Redshift.reshape(self.X_PIXEL_COUNT,self.Y_PIXEL_COUNT)
        NT_Redshift     = np.flip(NT_Redshift, 0)

        NT_Flux_Shifted = self.NT_Flux_Shifted.reshape(self.X_PIXEL_COUNT,self.Y_PIXEL_COUNT)
        NT_Flux_Shifted = np.flip(NT_Flux_Shifted, 0)

        Metadata = self.OBS_DISTANCE, self.OBS_INCLICATION, self.WINDOW_LIMITS, self.Legend
        
        return Intensity, NT_Flux, NT_Redshift, NT_Flux_Shifted, Metadata
    
    def export_ehtim_data(self, data: np.array):

        ehtim_x_fov = 2 * 7.500000e-05
        ehtim_y_fov = 2 * 7.500000e-05

        R_M87_LY = 53490000
        LY_TO_M  = 9.461e+15
        R_M87_M  = R_M87_LY * LY_TO_M

        M_SUN_KG = 2e30
        M_M87_KG = 6.5e9 * M_SUN_KG

        G = 6.67e-11
        c = 3e8

        R_SCH_M87 = M_M87_KG * G / c**2

        R_OBJ_DIMENTIONLESS = R_M87_M / R_SCH_M87

        RAD_TO_ARCS     = 206265
        CARTESIAN_TO_AS = RAD_TO_ARCS / R_OBJ_DIMENTIONLESS
        
        dx = (max(self.X_coords) - min(self.X_coords)) * CARTESIAN_TO_AS
        dy = (max(self.Y_coords) - min(self.Y_coords)) * CARTESIAN_TO_AS

        stupid_eht_scaling = 1e9

        formatted_sim_data = data.reshape(1, self.X_PIXEL_COUNT * self.Y_PIXEL_COUNT).flatten()

        array_to_export = np.array([self.X_coords / max(self.X_coords) * ehtim_x_fov / 2, 
                                    self.Y_coords / max(self.Y_coords) * ehtim_y_fov / 2, 
                                    formatted_sim_data * dx * dy / R_OBJ_DIMENTIONLESS**2 * stupid_eht_scaling]).T

        header = ("SRC: M87 \n"                   + 
                  "RA: 12 h 30 m 49.3920 s \n"    +
                  "DEC: 12 deg 23 m 27.9600 s \n" +
                  "MJD: 58211.000000 \n"          + 
                  "RF: 230.0000 GHz \n"           +
                  "FOVX: {} pix 0.000150 as \n".format(self.X_PIXEL_COUNT) +
                  "FOVY: {} pix 0.000150 as \n".format(self.Y_PIXEL_COUNT) +
                  "------------------------------------ \n" +
                  "x (as)     y (as)       I (Jy/pixel)")

        with open('data_for_ehtim.csv', 'w') as my_file:
            np.savetxt(my_file, array_to_export, fmt = '%0.4e', header = header)

        print('Array exported to file')

class ehtim_Parser():

    def __init__(self, File_name: str) -> None:

        def get_csv_line(path, line_number) -> str:
            with open(path) as file:
                return next(islice(csv.reader(file), line_number, None))[0]

        with open("Sim_Results\\" + File_name + ".txt", 'r') as file:

            self.HEADER_ROW_COUNT = 9

            csvreader = csv.reader(file, delimiter = " ")
            row_count = sum(1 for row in file) - self.HEADER_ROW_COUNT
            file.seek(0)
    
            self.X_PIXEL_COUNT = int(get_csv_line("Sim_Results\\" + File_name + ".txt", 5).split(" ")[2])
            self.Y_PIXEL_COUNT = int(get_csv_line("Sim_Results\\" + File_name + ".txt", 6).split(" ")[2])

            X_range = float(get_csv_line("Sim_Results\\" + File_name + ".txt", 5).split(" ")[4])
            Y_range = float(get_csv_line("Sim_Results\\" + File_name + ".txt", 6).split(" ")[4])

            self.WINDOW_LIMITS = [-X_range, X_range, -Y_range, Y_range]

            self.Legend = get_csv_line("Sim_Results\\" + File_name + ".txt", 8).split("  ")
            self.Legend = [self.Legend[0][2:8], self.Legend[2][1:7], self.Legend[5][1:13]]
            
            next(islice(csvreader, 10, None))

            self.X_coords        = np.zeros(row_count)
            self.Y_coords        = np.zeros(row_count)
            self.Intensity       = np.zeros(row_count)

            index = 0

            for row in csvreader:

                self.X_coords[index]  = row[0]
                self.Y_coords[index]  = row[1]
                self.Intensity[index] = row[2]

                index += 1

    def get_plottable_ehtim_data(self) -> tuple:

        # Arrays need to be flipped, because mpl treats y = 0 as the top, 
        # and the simulator (aka openGL) treats it as the bottom

        Intensity = self.Intensity.reshape(self.X_PIXEL_COUNT, self.Y_PIXEL_COUNT)
        # Intensity = np.flip(Intensity, 0)

        Metadata = self.WINDOW_LIMITS, self.Legend
        
        return Intensity, Metadata
    
def add_Wormhole_Shadow(spin: float, alpha: float, obs_distance: float, obs_inclanation: float, figure: plt) -> None:

    from scipy.optimize import fsolve

    def get_integrals_of_motion(r_ph: np.ndarray):

        omega = 2*spin/r_ph**3
        dr_omega = -3*omega/r_ph

        N    = np.exp(-1/r_ph - alpha/r_ph**2)
        dr_N = (1/r_ph**2 + 2*alpha/r_ph**3)*N
        
        Sigma = dr_N/N - 1/r_ph

        ksi = Sigma/(Sigma*omega - dr_omega)
        eta = r_ph**2/N**2*(1 - omega*ksi)**2 

        return ksi.flatten(), eta.flatten()

    def get_retrograde_photon_orbit() -> float:

        def func(r_ph: float):

            #--- Shadow Rim Parametric Equations ---#

            ksi, eta = get_integrals_of_motion(r_ph)

            #--- Square of the y coordinate of the shadow, viewed from the equator    ---#
            #--- The root of this functions gives the outter photon orbit radius r_ph ---#

            return eta - ksi**2
        
        return fsolve(func, 10)

    def get_shadow_branches_intersection() -> tuple:

        def func(r_ph: float):

            #--- Shadow Rim Parametric Equations ---#

            ksi, eta = get_integrals_of_motion(r_ph)

            N_th = np.exp(-1/r_throat - alpha/r_throat**2)
            omega_th = 2*spin/r_throat**3

            #--- The above calculates the integrals of motion, using the parametric equations for photon orbits outside the throat,
            #--- then imposes the condition they must satisfy ON the throat. This finds the intersection points (parametrized via r_ph) ---#

            return eta*N_th**2/r_throat**2 - (1-omega_th*ksi)**2
        
        return fsolve(func, 10)    

    def generate_shadow_rim(ksi: np.ndarray, eta: np.ndarray) -> tuple:

        #--- Metric Functions at the Observer ---#

        N_obs = np.exp(-1/obs_distance - alpha/obs_distance**2)
        omega_obs = 2*spin/obs_distance**3

        sin_0 = np.sin(obs_inclanation)

        g_tt   = -N_obs**2 + omega_obs**2*obs_distance**2*sin_0**2
        g_thth = obs_distance**2
        g_phph = obs_distance**2*sin_0**2
        g_tph  = omega_obs*obs_distance**2*sin_0**2

        g2    = g_tph**2 - g_tt*g_phph
        gamma = -g_tph/g_phph*np.sqrt(g_phph/g2)
        ksi_metric   = np.sqrt(g_phph/g2)

        Rad_potential   = (-eta*N_obs**2/obs_distance**2 + (1 - omega_obs * ksi)**2)
        Rad_potential = Rad_potential + (Rad_potential < 0)*1e10

        Theta_potential = (eta - ksi**2/sin_0**2)
        Theta_potential = Theta_potential * (Theta_potential > 0)

        #--- Local Contravariant Momenta ---#

        p_r = np.sqrt(Rad_potential)/N_obs
        p_theta = np.sqrt(Theta_potential)/np.sqrt(g_thth)
        p_phi = ksi/np.sqrt(g_phph)
        p_t = ksi_metric - gamma*ksi

        #--- Celestial Angular Coordinates ---#

        alpha_args = p_phi[(p_theta/p_t)**2 < 1] / p_r[(p_theta/p_t)**2 < 1]
        beta_args  = p_theta[(p_theta/p_t)**2 < 1] / p_t[(p_theta/p_t)**2 < 1]

        alpha_coord = np.arctan(alpha_args)
        beta_coord  = np.arcsin(beta_args)

        return alpha_coord.flatten(), beta_coord.flatten()

    def generate_shadow_due_to_photons_outside_throat(r_ph: np.ndarray) -> tuple:

        #--- Shadow Rim Parametric Equations due to Photons Orbiting Outside the Throat ---#

        ksi, eta = get_integrals_of_motion(r_ph)

        alpha_coord, beta_coord = generate_shadow_rim(ksi,eta)

        #--- Get rid of the weird arc that goes "backwards" and doesnt contribute to the shadow ---#
    
        if spin > 0:
            index = np.argmax(alpha_coord)
        else:
            index = np.argmin(alpha_coord)

        alpha_coord = alpha_coord[index:]
        beta_coord  = beta_coord[index:]
        ksi = ksi[index:]

        return alpha_coord, beta_coord, ksi

    def generate_shadow_due_to_photons_on_throat() -> tuple:

        ksi_max = ksi_out[np.argmax(beta_out)]
    
        N_th = np.exp(-1/r_throat - alpha/r_throat**2)
        omega_th = 2*spin/r_throat**3

        A = N_th**2/r_throat**2/np.sin(obs_inclanation)**2-omega_th**2
        B = 2*omega_th
        C = -1

        roots = np.roots([A,B,C])

        if (roots[0]*roots[0] > roots[1]*roots[1]):
            ksi_0 = roots[0]
        else:
            ksi_0 = roots[1]

        ksi = np.linspace(ksi_0, ksi_max, 10000)

        eta = r_throat**2/N_th**2*(1-omega_th*ksi)**2

        alpha_th, beta_th = generate_shadow_rim(ksi = ksi, eta = eta)

        return alpha_th, beta_th

    r_throat = 1
    r_ph_max = get_retrograde_photon_orbit()
    r_ph = np.linspace(r_throat, r_ph_max, 10000)

    alpha_out, beta_out, ksi_out = generate_shadow_due_to_photons_outside_throat(r_ph)

    alpha_th, beta_th = generate_shadow_due_to_photons_on_throat()

    #--- Get the itnersection points between the two shadow branches ---#

    alpha_int, beta_int, ksi_int = generate_shadow_due_to_photons_outside_throat(get_shadow_branches_intersection())

    #--- Delete the branches of the curves that do not contribute to the shadow ---#

    if spin > 0:

        throat_bitmask  = beta_th**2  < beta_int**2
        outside_bitmask = alpha_out < alpha_int

    else:

        throat_bitmask  = beta_th**2  < beta_int**2
        outside_bitmask = alpha_out > alpha_int

    alpha_out = alpha_out[outside_bitmask]
    beta_out = beta_out[outside_bitmask]

    alpha_th = alpha_th[throat_bitmask]
    beta_th = beta_th[throat_bitmask]

    alpha_shadow = np.concatenate([alpha_th, alpha_out,  np.flip(alpha_out), np.flip(alpha_th)])
    beta_shadow  = np.concatenate([beta_th , beta_out , -np.flip(beta_out), -np.flip(beta_th) ])

    alpha_shadow = alpha_shadow[beta_shadow != 0]
    beta_shadow  = beta_shadow[beta_shadow != 0]

    figure.plot(alpha_shadow, beta_shadow, "b")

def add_Kerr_Shadow(spin: float, obs_distance: float, obs_inclanation: float, figure: plt) -> None:

    #--- Equatorial Photon Orbits ---#

    r_ph_max = 2*(1 + np.cos(2/3 * np.arccos(-spin)))
    r_ph_min = 2*(1 + np.cos(2/3 * np.arccos( spin)))

    r_ph = np.linspace(r_ph_min, r_ph_max, 10000)

    #--- Critical Integrals of Motion ---#

    J = (spin**2 * (r_ph + 1) + r_ph**3 - 3*r_ph**2)/(spin*(r_ph - 1))
    K = (-9*r_ph**4 + 4*spin**2*r_ph**3 + 6*r_ph**5 - r_ph**6)/(r_ph - 1)**2/spin**2

    cos_0 = np.cos(obs_inclanation)
    sin_0 = np.sin(obs_inclanation)
    tan_0 = np.tan(obs_inclanation)

    #--- Metric at the Observer ---#
    g_phph = (obs_distance**2 + spin**2 + 2*obs_distance*spin**2*sin_0**2/(obs_distance**2 + spin**2*cos_0**2))*sin_0**2
    g_tph  = -2*spin*obs_distance*sin_0**2/(obs_distance**2 + spin**2*cos_0**2)
    g_tt   = -(1 - 2*obs_distance/(obs_distance**2 + spin**2*cos_0**2))

    g2    = g_tph**2 - g_tt*g_phph
    gamma = -g_tph/g_phph*np.sqrt(g_phph/g2)
    ksi   = np.sqrt(g_phph/g2)

    rho2  = obs_distance**2 + spin**2*cos_0**2
    delta = obs_distance**2 - 2*obs_distance + spin**2 

    Rad_potential = obs_distance**4 + (spin**2 - K - J**2)*obs_distance**2 + 2*obs_distance*(K + (spin - J)**2) - K*spin**2
    
    #--- Celestial Coordinates of the Shadow Rim ---#

    alpha = np.arctan(J * np.sqrt(rho2 * delta/Rad_potential/g_phph))

    #--- This should be non-negative on the shadow rim ---#
    problem_term = (K + spin**2*cos_0**2 - J**2/tan_0**2) * (K + spin**2*cos_0**2 - J**2/tan_0**2 > 0)
    beta = np.arcsin(np.sqrt(problem_term)/np.sqrt(obs_distance**2 + spin**2*cos_0**2)/(ksi - J*gamma)) 

    alpha = alpha[(beta != 0)]
    beta  = beta [(beta != 0)]

    alpha = np.concatenate([alpha, np.flip(alpha), [alpha[0]]])
    beta  = np.concatenate([beta ,-np.flip(beta),  [beta[0]] ])

    figure.plot(alpha,beta,"k")

def add_JNW_Shadow() -> None:

    #TODO

    pass

def add_celestial_sphere_pattern() -> None:

    #TODO

    pass

def overlay_thin_disk_images() -> None:

    #TODO

    pass 

# SIM_CASE = 5

Sim_Parser_0   = Simulation_Parser("Wormhole_Data0")
Sim_Parser_1   = Simulation_Parser("Wormhole_Data1")
Sim_Parser_2   = Simulation_Parser("Wormhole_Data2")
# ehtim_parser = ehtim_Parser("Case {}/Kerr_Case_{}_ehtim_blur".format(SIM_CASE,SIM_CASE))
# ehtim_parser_wh = ehtim_Parser("Wormhole_20_deg_tight_blur")

Intensity_0, NT_Flux_0, NT_Redshift_0, NT_Flux_Shifted_0, Metadata_0 = Sim_Parser_0.get_plottable_sim_data()
Intensity_1, NT_Flux_1, NT_Redshift_1, NT_Flux_Shifted_1, Metadata_1 = Sim_Parser_1.get_plottable_sim_data()
Intensity_2, NT_Flux_2, NT_Redshift_2, NT_Flux_Shifted_2, Metadata_2 = Sim_Parser_2.get_plottable_sim_data()

# Sim_Parser.export_ehtim_data(Intensity)

# ehtim_Intensity, ehtim_Metadata                            = ehtim_parser.get_plottable_ehtim_data()
# ehtim_Intensity_wh, ehtim_Metadata_wh                      = ehtim_parser_wh.get_plottable_ehtim_data()

axes_limits       = np.array([(limit) for limit in Metadata_0[2]])

# ehtim_axes_limits = np.array([(limit) for limit in ehtim_Metadata[0]])
# ehtim_axes_limits_wh = np.array([(limit) for limit in ehtim_Metadata_wh[0]])

Data_to_plot = Intensity_1
sim_figure   = plt.imshow(Data_to_plot, interpolation = 'bilinear', cmap = 'hot', extent = axes_limits, vmin = 0, vmax = np.max(Data_to_plot))
# ehtim_figure = plt.imshow(ehtim_Intensity, interpolation = 'bilinear', cmap = 'hot', extent = ehtim_axes_limits)

# figure = plt.figure()

# subfig_sim = figure.add_subplot(1, 2, 1)
# imgplot = plt.imshow(Intensity, interpolation = 'bilinear', cmap = 'hot', extent = axes_limits)
# subfig_sim.set_title('Ray Tracer')
# subfig_sim.set_xlabel(r'$\alpha$ [rad]')
# subfig_sim.set_ylabel(r'$\beta$ [rad]')

# subfig_eht = figure.add_subplot(1, 2, 2)
# imgplot = plt.imshow(ehtim_Intensity, interpolation = 'bilinear', cmap = 'hot', extent = ehtim_axes_limits*1e6)
# subfig_eht.set_title('Eht Imager')
# subfig_eht.set_xlabel(r'$\alpha$ [$\mu$as]')
# subfig_eht.set_ylabel(r'$\beta$  [$\mu$as]')

# # add_Kerr_Shadow(0.998, obs_distance = Metadata[0], obs_inclanation = np.deg2rad(Metadata[1]))
# # add_Wormhole_Shadow(spin = -0.98, alpha = 2, obs_distance = Metadata[0], obs_inclanation = np.deg2rad(Metadata[1]))

# #---- Figure Labels -----#

plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
plt.show()

"""
    Image cases:

        Case 1: 
            @ Thin disk - DISK_HEIGHT_SCALE = 100. / 3
            @ Radial Scale = 10 [M]
            @ Inclination = 80 [deg]
            @ Obs Window - X: [-30 30] [M], Y: [-30 30] [M]
            @ Absorbtion Coeff = 0
            @ Spin = 0.9
            @ Wormhole Alpha = 2

        Case 2: 
            @ Thick disk - DISK_HEIGHT_SCALE = 10. / 3
            @ Radial Scale = 10 [M]
            @ Inclination = 80 [deg]
            @ Obs Window - X: [-30 30] [M], Y: [-30 30] [M]
            @ Absorbtion Coeff = 0
            @ Spin = 0.9
            @ Wormhole Alpha = 2

        Case 3: 
            @ Thick disk - DISK_HEIGHT_SCALE = 10. / 3
            @ Radial Scale = 10 [M]
            @ Inclination = 20 [deg]
            @ Obs Window - X: [-30 30] [M], Y: [-30 30] [M]
            @ Absorbtion Coeff = 0
            @ Spin = 0.9
            @ Wormhole Alpha = 2

        Case 4: 
            @ Thick disk - DISK_HEIGHT_SCALE = 10. / 3
            @ Radial Scale = 5 [M]
            @ Inclination = 20 [deg]
            @ Obs Window - X: [-30 30] [M], Y: [-30 30] [M]
            @ Absorbtion Coeff = 10e5
            @ Spin = 0.9
            @ Wormhole Alpha = 2

        Case 5: 
            @ Thick disk - DISK_HEIGHT_SCALE = 10. / 3
            @ Radial Scale = 5 [M]
            @ Inclination = 20 [deg]
            @ Obs Window - X: [-30 30] [M], Y: [-30 30] [M]
            @ Absorbtion Coeff = 1e3
            @ Spin = 0.9
            @ Wormhole Alpha = 2
            
"""
