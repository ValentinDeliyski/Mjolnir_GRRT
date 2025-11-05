from _csv import reader, Reader, Writer
from csv import DictReader

from dataclasses import dataclass

from numpy import array, zeros, sum, flip, linspace, vstack, repeat, savetxt, log, pi, float64, sqrt
from numpy.typing import NDArray

class Simulation_Parser():

    def __init__(self, File_name: str) -> None:

        self.Simulation_metadata: dict = {}
        self.Raw_simulation_header: str = ""
        
        with open(File_name + ".txt", 'r') as file:

            Header_parser: Reader = reader(file, delimiter = ":")
            
            for line in Header_parser:
                
                for String in line: 
                    self.Raw_simulation_header = self.Raw_simulation_header + String + ":"
                    
                self.Raw_simulation_header = self.Raw_simulation_header[:-1] + "\n"
                
                if len(line) == 2:
                    self.Simulation_metadata.update({str(line[0]).strip(): str(line[1]).strip()})
                
                if -1 != line[0].find("Simulation Results"):
                    break

            if (int(self.Simulation_metadata["Active Simulation Mode"]) == 0):
     
                X_resolution: int = int(self.Simulation_metadata["Simulation Resolution"].split(" ")[0])
                Y_resolution: int = int(self.Simulation_metadata["Simulation Resolution"].split(" ")[2])
                
                Array_size = X_resolution * Y_resolution    
            
                self.X_coords: NDArray[float64] = zeros(Array_size)
                self.Y_coords: NDArray[float64] = zeros(Array_size)
                self.I_Intensity: NDArray[float64] = zeros(Array_size)
                self.Q_Intensity: NDArray[float64] = zeros(Array_size)
                self.U_Intensity: NDArray[float64] = zeros(Array_size)
                self.V_Intensity: NDArray[float64] = zeros(Array_size)

                self.Final_t_coord: NDArray[float64] = zeros(Array_size)

                self.Disk_redshift: NDArray[float64] = zeros(Array_size)
                self.Disk_flux: NDArray[float64] = zeros(Array_size)
                self.Source_t: NDArray[float64] = zeros(Array_size)
                self.Source_r: NDArray[float64] = zeros(Array_size)
                self.Source_phi: NDArray[float64] = zeros(Array_size)
                self.Source_p_r: NDArray[float64] = zeros(Array_size)
                self.Source_p_theta: NDArray[float64] = zeros(Array_size)
                self.Source_p_phi: NDArray[float64] = zeros(Array_size)

                self.Polarization_vec_X: NDArray[float64] = zeros(Array_size)
                self.Polarization_vec_Y: NDArray[float64] = zeros(Array_size)

                self.Celestial_theta: NDArray[float64] = zeros(Array_size)
                self.Celestial_phi: NDArray[float64] = zeros(Array_size)

                Data_parser = DictReader(file, delimiter = ",")
                index = 0
                
                for row in Data_parser:

                    try:

                        self.X_coords[index] = float(row["Image X Coord [M]"])
                        self.Y_coords[index] = float(row["Image Y Coord [M]"])

                        if ("Novikov-Thorne" == self.Simulation_metadata["Active disk model"]):
                            
                            self.Disk_redshift[index] = float(row["Disk Redshift [-]"])
                            self.Disk_flux[index] = float(row["Disk Flux [M_dot/M^2]"])

                            self.Source_t[index] = float(row["Source t Coord [M]"])
                            self.Source_r[index] = float(row["Source r Coord [M]"])
                            self.Source_phi[index] = float(row["Source phi Coord [Rad]"])
                            self.Source_p_r[index] = float(row["Radial Momentum (covariant)"])
                            self.Source_p_theta[index] = float(row["Theta Momentum (covariant)"])
                            self.Source_p_phi[index]  = float(row["Phi Momentum (covariant)"])

                            self.Polarization_vec_X[index] = float(row["Polarization vector X [-]"])
                            self.Polarization_vec_Y[index] = float(row["Polarization vector Y [-]"])
                            
                        else:

                            self.I_Intensity[index]  = float(row["Synchotron Intensity I [Jy/sRad]"])
                            self.Q_Intensity[index]  = float(row["Synchotron Intensity Q [Jy/sRad]"])
                            self.U_Intensity[index]  = float(row["Synchotron Intensity U [Jy/sRad]"])
                            self.V_Intensity[index]  = float(row["Synchotron Intensity V [Jy/sRad]"])       
                            
                            self.Final_t_coord[index] = float(row["Final t Coordinate [M]"])  
                            
                            self.Celestial_theta[index] = float(row["Celestial Sphere Crossing Theta [Rad]"])       
                            self.Celestial_phi[index]  = float(row["Celestial Sphere Crossing Phi [Rad]"])
                            
                        index += 1

                    except:
                        break
        
            if (int(self.Simulation_metadata["Active Simulation Mode"]) == 2):

                self.t_coord: list[float]     = []
                self.r_coord: list[float]     = []
                self.theta_coord: list[float] = []
                self.phi_coord: list[float]   = []
                
                self.p_t: list[float]     = []
                self.p_r: list[float]     = []
                self.p_theta: list[float] = []
                self.p_phi: list[float]   = []
                
                self.integration_step: list[float] = []
                self.affine_param: list[float] = []
                
                self.I_Intensity_log: list[float] = []
                self.Q_Intensity_log: list[float] = []
                self.U_Intensity_log: list[float] = []
                self.V_Intensity_log: list[float] = []
                self.State_error_log: list[float] = []
                self.Rejected_steps_log: list[int] = []
                
                Data_parser = DictReader(file, delimiter = ",")
                
                for row in Data_parser:
                    
                    self.t_coord.append(float(row["t_coord [M]"]))
                    self.r_coord.append(float(row["r_coord [M]"]))
                    self.theta_coord.append(float(row["theta_coord [rad]"]))
                    self.phi_coord.append(float(row["phi_coord [rad]"]))

                    self.p_t.append(float(row["p_t [-]"]))
                    self.p_r.append(float(row["p_r [-]"]))
                    self.p_theta.append(float(row["p_theta [rad/M]"]))
                    self.p_phi.append(float(row["p_phi [rad/M]"]))

                    self.integration_step.append(float(row["Integration Step [M]"]))
                    self.affine_param.append(float(row["Affine Parameter [M]"]))

                    self.I_Intensity_log.append(float(row["Synchotron Intensity I [Jy/sRad]"]))
                    self.Q_Intensity_log.append(float(row["Synchotron Intensity Q [Jy/sRad]"]))
                    self.U_Intensity_log.append(float(row["Synchotron Intensity U [Jy/sRad]"]))
                    self.V_Intensity_log.append(float(row["Synchotron Intensity V [Jy/sRad]"]))
                    self.State_error_log.append(float(row["State Error [-]"]))
                    
                    try:
                        self.Rejected_steps_log.append(int(row["Number of rejected steps [-]"]))    
                    except:
                        self.Rejected_steps_log.append(int(float(row["Number of rejected steps [-]"])))   
                
    def get_total_flux(self, obs_pos: float, unit: str = "Jy") -> float:
        
        """ The observation window limits are given in geometric length units, 
            so one divides by the effective observer distance to get the angular size. """
            
        Window_limits = self.Simulation_metadata["Observation Window Dimentions (-X,+X,-Y,+Y) [M]"].split(",")
        Window_limits = [float(Limit) for Limit in Window_limits]
        
        X_resolution: int = int(self.Simulation_metadata["Simulation Resolution"].split(" ")[0])
        Y_resolution: int = int(self.Simulation_metadata["Simulation Resolution"].split(" ")[2])
        
        Pixel_area: float = (abs(Window_limits[1] - Window_limits[0]) * 
                             abs(Window_limits[3] - Window_limits[2]) / X_resolution / Y_resolution / obs_pos**2)

        """ The base flux unit, returned by the ray-tracer is Jy. """
        Total_Intensity_Jy: float = float(sum(self.I_Intensity) * Pixel_area)

        match unit:
            
            case "Jy":
                return Total_Intensity_Jy
            
            case "mJy":
              return 1e3 * Total_Intensity_Jy
        
            case _:
                print("Unsupported flux unit!")    
                return 0
            
    def get_polarized_flux(self, obs_pos: float, unit: str = "Jy") -> float:
        
        """ The observation window limits are given in geometric length units, 
            so one divides by the effective observer distance to get the angular size. """
            
        Window_limits = self.Simulation_metadata["Observation Window Dimentions (-X,+X,-Y,+Y) [M]"].split(",")
        Window_limits = [float(Limit) for Limit in Window_limits]
        
        X_resolution: int = int(self.Simulation_metadata["Simulation Resolution"].split(" ")[0])
        Y_resolution: int = int(self.Simulation_metadata["Simulation Resolution"].split(" ")[2])
        
        Pixel_area: float = (abs(Window_limits[1] - Window_limits[0]) * 
                             abs(Window_limits[3] - Window_limits[2]) / X_resolution / Y_resolution / obs_pos**2)

        """ The base flux unit, returned by the ray-tracer is Jy. """
        Total_Intensity_Jy: float = float(sum(sqrt(self.Q_Intensity**2 + self.U_Intensity**2 + self.V_Intensity**2)) * Pixel_area)

        match unit:
            
            case "Jy":
                return Total_Intensity_Jy
            
            case "mJy":
              return 1e3 * Total_Intensity_Jy
        
            case _:
                print("Unsupported flux unit!")    
                return 0

    def get_plottable_sim_data(self) -> tuple[NDArray[float64], ...]:

        """ The arrays first need to be reshaped into 2D ones, then flipped along the x axis, 
            because mpl treats y = 0 as the top, and the ray-tracer (openGL) treats it as the bottom. """
        
        X_resolution: int = int(self.Simulation_metadata["Simulation Resolution"].split(" ")[0])
        Y_resolution: int = int(self.Simulation_metadata["Simulation Resolution"].split(" ")[2])
            
        I_Intensity = self.I_Intensity.reshape(Y_resolution, X_resolution)
        I_Intensity = flip(I_Intensity, axis = 0)

        Q_Intensity = self.Q_Intensity.reshape(Y_resolution, X_resolution)
        Q_Intensity = flip(Q_Intensity, axis =  0)

        U_Intensity = self.U_Intensity.reshape(Y_resolution, X_resolution)
        U_Intensity = flip(U_Intensity, axis =  0)

        V_Intensity = self.V_Intensity.reshape(Y_resolution, X_resolution)
        V_Intensity = flip(V_Intensity, axis =  0)

        Disk_flux   = self.Disk_flux.reshape(Y_resolution, X_resolution)
        Disk_flux     = flip(Disk_flux, axis =  0)

        Disk_redshift = self.Disk_redshift.reshape(Y_resolution, X_resolution)
        Disk_redshift   = flip(Disk_redshift, axis =  0)
        
        Celestial_theta = self.Celestial_theta.reshape(Y_resolution, X_resolution)
        Celestial_theta = flip(Celestial_theta, axis =  0)
        
        Celestial_phi = self.Celestial_phi.reshape(Y_resolution, X_resolution)
        Celestial_phi = flip(Celestial_phi, axis =  0)

        return I_Intensity, Q_Intensity, U_Intensity, V_Intensity, Disk_redshift, Disk_flux, Celestial_theta, Celestial_phi
    
    def get_photon_log(self) -> tuple[tuple, tuple, tuple, list, list, tuple]:
        
        Position_tuple = self.t_coord, self.r_coord, self.theta_coord, self.phi_coord
        Momentum_tuple = self.p_t, self.p_r, self.p_theta, self.p_phi
        Emission_tuple = self.I_Intensity_log, self.Q_Intensity_log, self.U_Intensity_log, self.V_Intensity_log
        Debug_tuple = self.State_error_log, self.Rejected_steps_log
        
        return Position_tuple, Momentum_tuple, Emission_tuple, self.integration_step, self.affine_param, Debug_tuple
        
    def export_ehtim_data(self, Spacetime: str, data: NDArray, path: str) -> None:

        ehtim_x_fov = 2 * 5.000000e-05
        ehtim_y_fov = 2 * 5.000000e-05
        
        Window_limits = self.Simulation_metadata["Observation Window Dimentions (-X,+X,-Y,+Y) [M]"].split(",")
        Window_limits = [float(Limit) for Limit in Window_limits]
        
        X_resolution: int = int(self.Simulation_metadata["Simulation Resolutoin"].split(" ")[0])
        Y_resolution: int = int(self.Simulation_metadata["Simulation Resolutoin"].split(" ")[2])
        
        Units = Units_class()
        
        Pixel_area = abs(Window_limits[1] - Window_limits[0]) * abs(Window_limits[3] - Window_limits[2]) / X_resolution / Y_resolution

        formatted_sim_data = data.reshape(1, X_resolution * Y_resolution).flatten()

        X_coords = linspace(-1, 1, X_resolution) * ehtim_x_fov / 2
        X_coords = vstack([X_coords] * Y_resolution).flatten()

        Y_coords = linspace(-1, 1, Y_resolution) * ehtim_y_fov / 2
        Y_coords = repeat(Y_coords, X_resolution, axis = 0)

        array_to_export = array([X_coords, 
                                 Y_coords, 
                                 formatted_sim_data * Pixel_area / Units.M87_DISTANCE_GEOMETRICAL**2]).T

        Obs_frequency: float = float(self.Simulation_metadata["Observation Frequency [Hz]"])

        header = ("SRC: M87 \n"                   + 
                  "RA: 12 h 30 m 49.3920 s \n"    +
                  "DEC: 12 deg 23 m 27.9600 s \n" +
                  "MJD: 58211.000000 \n"          + 
                  "RF: {} GHz \n".format(Obs_frequency / 1e9)    +
                  "FOVX: {} pix 0.000100 as \n".format(X_resolution) +
                  "FOVY: {} pix 0.000100 as \n".format(Y_resolution) +
                  "------------------------------------ \n" +
                  "x (as)     y (as)       I (Jy/pixel)")

        with open(path + '{}_data_for_ehtim_{}.csv'.format(Spacetime, int(Obs_frequency / 1e9)), 'w') as my_file:
                  savetxt(my_file, array_to_export, fmt = '%0.4e', header = header)

        print('Array exported to file!')

class ehtim_Parser():

    def __init__(self, File_name: str) -> None:

        with open(File_name + ".txt", 'r') as file:

            csvreader = reader(file, delimiter = " ")

            for _ in range(4):
                    _ = csvreader.__next__()

            self.OBS_FREQUENCY = float(csvreader.__next__()[2])

            X_data_line = csvreader.__next__()
    
            self.X_PIXEL_COUNT = int(X_data_line[2])
            X_range            = float(X_data_line[4]) / 2

            Y_data_line = csvreader.__next__()

            self.Y_PIXEL_COUNT = int(Y_data_line[2])
            Y_range            = float(Y_data_line[4]) / 2

            self.WINDOW_LIMITS = [-X_range, X_range, -Y_range, Y_range]

            for _ in range(2):
                    _ = csvreader.__next__()

            self.X_coords  = zeros(self.X_PIXEL_COUNT * self.Y_PIXEL_COUNT)
            self.Y_coords  = zeros(self.X_PIXEL_COUNT * self.Y_PIXEL_COUNT)
            self.Intensity = zeros(self.X_PIXEL_COUNT * self.Y_PIXEL_COUNT)

            index = 0

            for row in csvreader:

                self.X_coords[index]  = float(row[0])
                self.Y_coords[index]  = float(row[1])
                self.Intensity[index] = float(row[2])

                index += 1

    def get_plottable_ehtim_data(self) -> tuple[NDArray, list]:

        Intensity = self.Intensity.reshape(self.X_PIXEL_COUNT, self.Y_PIXEL_COUNT)

        Metadata = self.WINDOW_LIMITS
        
        return Intensity, Metadata
    
    def get_total_flux(self) -> float64:

        return sum(self.Intensity)

@dataclass
class VIDA_params_Parser():

    d0: float
    Sigma: float
    Tau: float
    rot_angle: float
    slash: float
    slash_angle: float
    x0: float
    y0: float
    div: float

    def __init__(self, File_name: str) -> None:

        with open(File_name + ".csv", 'r') as file:

            csvreader = reader(file, delimiter = " ")
            
            self.d0          = 2 * float(csvreader.__next__()[0])
            self.Sigma       = float(csvreader.__next__()[0]) 
            self.Tau         = float(csvreader.__next__()[0])
            self.rot_angle   = float(csvreader.__next__()[0])
            self.slash       = float(csvreader.__next__()[0])
            self.slash_angle = float(csvreader.__next__()[0])
            self.x0          = float(csvreader.__next__()[0])
            self.y0          = float(csvreader.__next__()[0])
            self.div         = float(csvreader.__next__()[0]) 

@dataclass
class Units_class():

    """ ============== Useful scaling constants ============== """

    KILO: float = 1e3
    MEGA: float = 1e6
    GIGA: float = 1e9
    
    """ =================== Physical constants ================ """
    
    C_LIGHT_SI: float   = 299792458
    G_NEWTON_SI: float  = 6.6743e-11
    BOLTZMANN_SI: float = 1.380649e-23
    PLANCK_SI: float    = 6.62607015e-34

    M_SUN_SI: float     = 1.988475e30
    M_M87_BH_SI: float  = 6.2e9 * M_SUN_SI
    M_SGRA_BH_SI: float = 4.297e6 * M_SUN_SI
    
    """ ================== Time conversions ================== """

    YEAR_TO_SEC: float = 31556952
    
    """ ================= Angular conversions ================ """

    RAD_TO_DEG: float      = 180 / pi
    ARCSEC_TO_RAD: float   = pi / 180 / 3600
    DEG_TO_AS: float       = 3600
    RAD_TO_MICRO_AS: float = RAD_TO_DEG * DEG_TO_AS * 1e6

    """ ================ Distance conversions ================ """

    LY_TO_METER: float      = C_LIGHT_SI * YEAR_TO_SEC
    PC_TO_METER: float      = 3.26156 * LY_TO_METER
    GR_MASS_TO_METER: float = G_NEWTON_SI / C_LIGHT_SI**2

    M87_DISTANCE_LY: float = 53.49e6
    M87_DISTANCE_PC: float = 16.9e6
    M87_DISTANCE_GEOMETRICAL: float = M87_DISTANCE_PC * PC_TO_METER / GR_MASS_TO_METER / M_M87_BH_SI

    SGRA_DISTANCE_LY: float = 26673
    SGRA_DISTANCE_PC: float = 8.277e3
    SGRA_DISTANCE_GEOMETRICAL = SGRA_DISTANCE_PC * PC_TO_METER / GR_MASS_TO_METER / M_SGRA_BH_SI

    """ ================== Flux conversions ================== """

    W_M2_TO_JY: float = 1e26
    J_TO_ERG: float   = 1e7
    
    """ ================== Misc =================="""

    M_ELECTRON_SI = 9.1093837e-31 
    
    def Spectral_density_to_T(self, I_nu: NDArray, frequency: float) -> NDArray:

        I_nu = abs(I_nu) + 1e-40 # To avoid division by 0 errors

        return self.PLANCK_SI * frequency / self.BOLTZMANN_SI / log(1 + 2 * self.PLANCK_SI * frequency**3 / self.C_LIGHT_SI**2 / (I_nu))
    