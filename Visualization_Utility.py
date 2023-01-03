import csv
import numpy as np
from matplotlib import pyplot as plt

def Parse_Sim_Results(File_name: str):

    with open(File_name + ".txt", 'r') as file:

        HEADER_ROW_COUNT = 4

        csvreader = csv.reader(file, delimiter = ":")
        row_count = sum(1 for row in file) - HEADER_ROW_COUNT
        file.seek(0)

        OBS_DISTANCE    = float(csvreader.__next__()[1])
        OBS_INCLICATION = float(csvreader.__next__()[1])
        X_PIXEL_COUNT   = 401
        Y_PIXEL_COUNT   = 401

        WINDOW_LIMITS = [float(limit) for limit in csvreader.__next__()[1].split(',')]

        Legend = csvreader.__next__()

        csvreader = csv.reader(file, delimiter = " ")

        NT_Flux         = np.zeros(row_count)
        Intensity       = np.zeros(row_count)
        NT_Redshift     = np.zeros(row_count)
        NT_Flux_Shifted = np.zeros(row_count)

        index = 0

        for row in csvreader:

            NT_Redshift[index]  = row[2]
            NT_Flux[index]      = row[3]
            Intensity[index]    = row[4]

            NT_Flux_Shifted[index] = NT_Redshift[index]**4*NT_Flux[index]

            index += 1

    # Arrays need to be flipped, because mpl treats y = 0 as the top, 
    # and the simulator (aka openGL) treats it as the bottom

    Intensity = Intensity.reshape(X_PIXEL_COUNT, Y_PIXEL_COUNT)
    Intensity = np.flip(Intensity, 0)

    NT_Flux         = NT_Flux.reshape(X_PIXEL_COUNT,Y_PIXEL_COUNT)
    NT_Flux         = np.flip(NT_Flux, 0)
    NT_Redshift     = NT_Redshift.reshape(X_PIXEL_COUNT,Y_PIXEL_COUNT)
    NT_Redshift     = np.flip(NT_Redshift, 0)
    NT_Flux_Shifted = NT_Flux_Shifted.reshape(X_PIXEL_COUNT,Y_PIXEL_COUNT)
    NT_Flux_Shifted = np.flip(NT_Flux_Shifted, 0)

    Metadata = (OBS_DISTANCE, OBS_INCLICATION, WINDOW_LIMITS, Legend)

    return Intensity, NT_Flux, NT_Redshift, NT_Flux_Shifted, Metadata

def add_Wormhole_Shadow(spin: float, alpha: float, obs_distance: float, obs_inclanation: float):

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

    def get_retrograde_photon_orbit(spin: float, alpha: float):

        def func(r_ph: float):

            #--- Shadow Rim Parametric Equations ---#

            ksi, eta = get_integrals_of_motion(r_ph)

            #--- Square of the y coordinate of the shadow, viewed from the equator    ---#
            #--- The root of this functions gives the outter photon orbit radius r_ph ---#

            return eta - ksi**2
        
        return fsolve(func, 10)

    def get_shadow_branches_intersection(spin: float, alpha: float):

        def func(r_ph: float):

            #--- Shadow Rim Parametric Equations ---#

            ksi, eta = get_integrals_of_motion(r_ph)

            N_th = np.exp(-1/r_throat - alpha/r_throat**2)
            omega_th = 2*spin/r_throat**3

            #--- The above calculates the integrals of motion, using the parametric equations for photon orbits outside the throat,
            #--- then imposes the condition they must satisfy ON the throat. This finds the intersection points (parametrized via r_ph) ---#

            return eta*N_th**2/r_throat**2 - (1-omega_th*ksi)**2
        
        return fsolve(func, 10)    

    def generate_shadow_rim(ksi: np.ndarray, eta: np.ndarray):

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

        Rad_potential   = (-eta*N_obs**2/obs_distance**2 + (1 - omega_obs * ksi)**2) * (-eta*N_obs**2/obs_distance**2 + (1 - omega_obs * ksi)**2 > 0)
        Theta_potential = (eta - ksi**2/sin_0**2) * (eta - ksi**2/sin_0**2 > 0)

        #--- Local Contravariant Momenta ---#

        p_r = np.sqrt(Rad_potential)/N_obs
        p_theta = np.sqrt(Theta_potential)/np.sqrt(g_thth)
        p_phi = ksi/np.sqrt(g_phph)
        p_t = ksi_metric - gamma*ksi

        #--- Celestial Angular Coordinates ---#

        alpha_coord = np.arctan(p_phi/p_r)
        beta_coord  = np.arcsin(p_theta/p_t)

        return alpha_coord.flatten(), beta_coord.flatten()

    def generate_shadow_due_to_photons_outside_throat(r_ph: np.ndarray):

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
        ksi  = ksi[index:]

        return alpha_coord, beta_coord, ksi

    def generate_shadow_due_to_photons_on_throat():

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
    r_ph_max = get_retrograde_photon_orbit(spin = spin, alpha = alpha)
    r_ph = np.linspace(r_throat, r_ph_max, 10000)

    alpha_out, beta_out, ksi_out = generate_shadow_due_to_photons_outside_throat(r_ph)

    alpha_th, beta_th = generate_shadow_due_to_photons_on_throat()

    #--- Get the itnersection points between the two shadow branches ---#

    alpha_int, beta_int, ksi_int = generate_shadow_due_to_photons_outside_throat(get_shadow_branches_intersection(spin = spin, alpha = alpha))

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

    plt.plot(alpha_shadow, beta_shadow, "k")

def add_Kerr_Shadow(spin: float, obs_distance: float, obs_inclanation: float):

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

    plt.plot(alpha,beta,"k")

def add_JNW_Shadow():

    #TODO

    pass

def add_celestial_sphere_pattern():

    #TODO

    pass

Intensity, NT_Flux, NT_Redshift, NT_Flux_Shifted, Metadata = Parse_Sim_Results("Wormhole_Data0")

axes_limits = [np.round(limit) for limit in Metadata[2]]
plt.imshow(Intensity, interpolation = 'nearest', cmap = 'jet', extent = axes_limits)

# add_Kerr_Shadow(0.998, Metadata[0], np.deg2rad(Metadata[1]))
add_Wormhole_Shadow(spin = -0.5, alpha = 4, obs_distance = 15, obs_inclanation = np.deg2rad(89.9))

#---- Figure Labels -----#

plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

plt.xlabel(r'$\alpha$ [rad]')
plt.ylabel(r'$\beta$ [rad]')
plt.show()
