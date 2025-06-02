from csv import reader
from dataclasses import dataclass
from numpy import array, zeros, sum, flip, linspace, vstack, repeat, savetxt, log, pi, float64
from numpy.typing import NDArray

class Simulation_Parser():

    def __init__(self, File_name: str) -> None:

        with open(File_name + ".txt", 'r') as file:

            csvreader = reader(file, delimiter = ":")

            _ = csvreader.__next__()

            self.metric = str(csvreader.__next__()[1])[1:]

            match self.metric:

                case "Kerr":
                    self.Spin = float(csvreader.__next__()[1][1:])
                case "Wormhole":
                    self.Spin = float(csvreader.__next__()[1][1:])
                    self.Redshift_parameter = float(csvreader.__next__()[1][1:])
                case "Janis_Newman_Winicour":
                    self.Gamma = float(csvreader.__next__()[1][1:])
                case "Einstein_Gauss_Bonnet":
                    self.Gamma = float(csvreader.__next__()[1][1:])
                case "Regular_Black_Hole":
                    self.Parameter = float(csvreader.__next__()[1][1:])
                case "BH_w_Dark_Matter_Halo":
                    self.Halo_mass = float(csvreader.__next__()[1][1:])
                    self.Halo_Compactness = float(csvreader.__next__()[1][1:])
     
            self.Active_Sim_Mode = int(csvreader.__next__()[1])

            _ = csvreader.__next__() # Image Order
            _ = csvreader.__next__() # Observer Params 

            self.OBS_TIME = float(csvreader.__next__()[1])
            self.OBS_DISTANCE    = float(csvreader.__next__()[1])
            self.OBS_INCLICATION = float(csvreader.__next__()[1])

            _ = csvreader.__next__() # Observer Azimuth

            self.OBS_FREQUENCY   = float(csvreader.__next__()[1])
            
            if(self.Active_Sim_Mode == 2):
                self.Photon_Number = int(csvreader.__next__()[1])
                self.Param_Sweep_Number = int(csvreader.__next__()[1])

            if (self.Active_Sim_Mode == 1):

                self.WINDOW_LIMITS = [float(limit) for limit in csvreader.__next__()[1].split(',')]
                
                Resolution_list = csvreader.__next__()[1].split(' ')

                self.X_PIXEL_COUNT   = int(Resolution_list[1])
                self.Y_PIXEL_COUNT   = int(Resolution_list[3])

            _ = csvreader.__next__() # Accretion Disk Parameters Header
            self.ACTIVE_DISK_MODEL = csvreader.__next__()[1] # Active model string
            
            _ = csvreader.__next__() # Model parameters header

            test = self.ACTIVE_DISK_MODEL[1:-2]

            if self.ACTIVE_DISK_MODEL[1:-2] == "Phenomenological_RIAF":
                self.disk_opening_angle = float(csvreader.__next__()[1])
                self.disk_density_power_law_scale = float(csvreader.__next__()[1])
                self.disk_density_power_law_power = float(csvreader.__next__()[1])
                self.disk_density_cutoff_radius = float(csvreader.__next__()[1])
                self.disk_density_cutoff_scale = float(csvreader.__next__()[1])
                
                self.disk_temperature_power_law_scale = float(csvreader.__next__()[1])
                self.disk_temperature_power_law_power = float(csvreader.__next__()[1])
                self.disk_temperature_cutoff_radius = float(csvreader.__next__()[1])
                self.disk_temperature_cutoff_scale = float(csvreader.__next__()[1])
                
                self.disk_max_density = float(csvreader.__next__()[1])
                self.disk_max_temperature = float(csvreader.__next__()[1])
                
                self.disk_ensamble = csvreader.__next__()[1]

            else:
                self.disk_density_exp_height_scale = float(csvreader.__next__()[1])
                self.disk_density_exp_radial_scale = float(csvreader.__next__()[1])

            _ = csvreader.__next__() # Magnetic field parameters header

            self.disk_magnetization = float(csvreader.__next__()[1])
            self.disk_magnetic_field = csvreader.__next__()[1]
            self.disk_magnetic_field_magnitude_profile = csvreader.__next__()[1]
            
            if("Power law based" == self.disk_magnetic_field_magnitude_profile[1:]):
                self.disk_magnetic_field_scale = float(csvreader.__next__()[1])
                self.magnetic_field_radial_scale = float(csvreader.__next__()[1])
                self.magnetic_field_power_law_power = float(csvreader.__next__()[1])

            _ = csvreader.__next__() # Hotspot Parameters Header
            _ = csvreader.__next__() # Density Model Parameters Header

            self.hotspot_density_profile = csvreader.__next__()[1][1:]
            if self.hotspot_density_profile == "Gaussian":
                self.hotspot_density_spread = float(csvreader.__next__()[1])
            elif self.hotspot_density_profile == "Spherical":
                self.hotspot_dentiy_radius = float(csvreader.__next__()[1])
            elif self.hotspot_temperature_profile == "Hybrid power law gaussian":
                self.hotspot_density_spread = float(csvreader.__next__()[1])  
                self.hotspot_density_power_law_power = float(csvreader.__next__()[1])  
                self.hotspot_density_power_law_scale = float(csvreader.__next__()[1])   
            
            self.hotspot_max_density = float(csvreader.__next__()[1])

            _ = csvreader.__next__() # Temperature Model Parameters Header

            self.hotspot_temperature_profile = csvreader.__next__()[1][1:]
            if self.hotspot_temperature_profile == "Gaussian":
                self.hotspot_temperature_spread = float(csvreader.__next__()[1])       
            elif self.hotspot_temperature_profile == "Spherical":
                self.hotspot_temperature_radius = float(csvreader.__next__()[1])
            elif self.hotspot_temperature_profile == "Hybrid power law gaussian":
                self.hotspot_temperature_spread = float(csvreader.__next__()[1])  
                self.hotspot_temperature_power_law_power = float(csvreader.__next__()[1])  
                self.hotspot_temperature_power_law_scale = float(csvreader.__next__()[1])   

            self.hotspot_max_temperature = float(csvreader.__next__()[1])

            _ = csvreader.__next__() # Hotspot Synchrotron Emission Model Parameters Header

            self.hotspot_ensamble = csvreader.__next__()[1][1:]
            if self.hotspot_ensamble == "Kappa":
                self.kappa = csvreader.__next__()[1][1:]
            elif self.hotspot_ensamble == "Phenomenological":
                for _ in range(4):
                    _ = csvreader.__next__()
            elif self.hotspot_ensamble == "Thermal":
                pass

            self.hotspot_magnetization = float(csvreader.__next__()[1])
            self.hotspot_magnetic_field = csvreader.__next__()[1][1:]

            _ = csvreader.__next__() # Hotspot Position Header

            self.hotspot_distance    = float(csvreader.__next__()[1])
            self.hotspot_inclination = float(csvreader.__next__()[1])
            self.hotspot_azimuth     = float(csvreader.__next__()[1])
            self.coord_time_offset   = float(csvreader.__next__()[1])

            _ = csvreader.__next__() # Novikov - Thorner Model Parameters Header

            try:
                self.NT_r_in  = float(csvreader.__next__()[1])
                self.NT_r_out = float(csvreader.__next__()[1])
            except:
                pass
            
            _ = csvreader.__next__() # Simulation Results Header
            self.Legend = csvreader.__next__()

            csvreader = reader(file, delimiter = " ")

            if (self.Active_Sim_Mode != 2):
                Array_size = self.X_PIXEL_COUNT * self.Y_PIXEL_COUNT
            else:
                Array_size = self.Photon_Number * self.Param_Sweep_Number

            self.X_coords        = zeros(Array_size)
            self.Y_coords        = zeros(Array_size)
            self.NT_Flux         = zeros(Array_size)
            self.I_Intensity     = zeros(Array_size)
            self.Q_Intensity     = zeros(Array_size)
            self.U_Intensity     = zeros(Array_size)
            self.V_Intensity     = zeros(Array_size)
            self.NT_Redshift     = zeros(Array_size)
            self.NT_Flux_Shifted = zeros(Array_size)

            self.Source_R_Coord   = zeros(Array_size)
            self.Source_Phi_Coord = zeros(Array_size)
            self.Radial_Momentum  = zeros(Array_size)
            self.Theta_Momentum   = zeros(Array_size)
            self.Phi_Momentum     = zeros(Array_size)
            self.Param_1          = zeros(Array_size)
            self.Param_2          = zeros(Array_size)

            index = 0

            for row in csvreader:

                try:

                    self.X_coords[index] = row[0]
                    self.Y_coords[index] = row[1]

                    self.NT_Redshift[index]  = row[2]
                    self.NT_Flux[index]      = row[3]
                    self.I_Intensity[index]  = row[4]
                    self.Q_Intensity[index]  = row[5]
                    self.U_Intensity[index]  = row[6]
                    self.V_Intensity[index]  = row[7]

                    if self.Active_Sim_Mode == 2:

                        self.Source_R_Coord[index]   = row[8]
                        self.Source_Phi_Coord[index] = row[9]
                        self.Radial_Momentum[index]  = row[10]
                        self.Theta_Momentum[index]   = row[11]
                        self.Phi_Momentum[index]     = row[12]
                        self.Param_1[index]          = row[13]

                        try:
                            self.Param_2[index] = row[14]
                        except:
                            self.Param_2[index] = 0

                    self.NT_Flux_Shifted[index] = self.NT_Redshift[index]**4 * self.NT_Flux[index]

                    index += 1

                except:
                    break
                
    def get_total_flux(self, obs_pos: float, unit: str = "Jy") -> float:
        
        """ The observation window limits are given in geometric length units, 
            so one divides by the effective observer distance to get the angular size. """
        Pixel_area: float = ((self.WINDOW_LIMITS[1] - self.WINDOW_LIMITS[0]) * 
                             (self.WINDOW_LIMITS[3] - self.WINDOW_LIMITS[2]) / self.X_PIXEL_COUNT / self.Y_PIXEL_COUNT / obs_pos**2)

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

    def get_plottable_sim_data(self) -> tuple[NDArray[float64], ...]:

        """ The arrays first need to be reshaped into 2D ones, then flipped along the x axis, 
            because mpl treats y = 0 as the top, and the ray-tracer (openGL) treats it as the bottom. """
            
        I_Intensity = self.I_Intensity.reshape(self.Y_PIXEL_COUNT, self.X_PIXEL_COUNT)
        I_Intensity = flip(I_Intensity, axis = 0)

        Q_Intensity = self.Q_Intensity.reshape(self.Y_PIXEL_COUNT, self.X_PIXEL_COUNT)
        Q_Intensity = flip(Q_Intensity, axis =  0)

        U_Intensity = self.U_Intensity.reshape(self.Y_PIXEL_COUNT, self.X_PIXEL_COUNT)
        U_Intensity = flip(U_Intensity, axis =  0)

        V_Intensity = self.V_Intensity.reshape(self.Y_PIXEL_COUNT, self.X_PIXEL_COUNT)
        V_Intensity = flip(V_Intensity, axis =  0)

        NT_Flux         = self.NT_Flux.reshape(self.Y_PIXEL_COUNT,self.X_PIXEL_COUNT)
        NT_Flux         = flip(NT_Flux, axis =  0)

        NT_Redshift     = self.NT_Redshift.reshape(self.Y_PIXEL_COUNT,self.X_PIXEL_COUNT)
        NT_Redshift     = flip(NT_Redshift, axis =  0)

        NT_Flux_Shifted = self.NT_Flux_Shifted.reshape(self.Y_PIXEL_COUNT,self.X_PIXEL_COUNT)
        NT_Flux_Shifted = flip(NT_Flux_Shifted, axis =  0)

        return I_Intensity, Q_Intensity, U_Intensity, V_Intensity, NT_Redshift, NT_Flux, NT_Flux_Shifted
    
    def export_ehtim_data(self, Spacetime: str, data: NDArray, path: str) -> None:

        ehtim_x_fov = 2 * 5.000000e-05
        ehtim_y_fov = 2 * 5.000000e-05

        Units = Units_class()
        
        Pixel_area = (self.WINDOW_LIMITS[1] - self.WINDOW_LIMITS[0]) * (self.WINDOW_LIMITS[3] - self.WINDOW_LIMITS[2]) / self.X_PIXEL_COUNT / self.Y_PIXEL_COUNT

        formatted_sim_data = data.reshape(1, self.X_PIXEL_COUNT * self.Y_PIXEL_COUNT).flatten()

        X_coords = linspace(-1, 1, self.X_PIXEL_COUNT) * ehtim_x_fov / 2
        X_coords = vstack([X_coords] * self.Y_PIXEL_COUNT).flatten()

        Y_coords = linspace(-1, 1, self.Y_PIXEL_COUNT) * ehtim_y_fov / 2
        Y_coords = repeat(Y_coords, self.X_PIXEL_COUNT, axis = 0)

        array_to_export = array([X_coords, 
                                 Y_coords, 
                                 formatted_sim_data * Pixel_area / Units.SGRA_DISTANCE_GEOMETRICAL**2]).T

        header = ("SRC: M87 \n"                   + 
                  "RA: 12 h 30 m 49.3920 s \n"    +
                  "DEC: 12 deg 23 m 27.9600 s \n" +
                  "MJD: 58211.000000 \n"          + 
                  "RF: {} GHz \n".format(self.OBS_FREQUENCY / 1e9)    +
                  "FOVX: {} pix 0.000100 as \n".format(self.X_PIXEL_COUNT) +
                  "FOVY: {} pix 0.000100 as \n".format(self.Y_PIXEL_COUNT) +
                  "------------------------------------ \n" +
                  "x (as)     y (as)       I (Jy/pixel)")

        with open(path + '{}_data_for_ehtim_{}.csv'.format(Spacetime, int(self.OBS_FREQUENCY / 1e9)), 'w') as my_file:
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

        I_nu += 1e-10 # To avoid division by 0 errors

        return self.PLANCK_SI * frequency / self.BOLTZMANN_SI / log(1 + 2 * self.PLANCK_SI * frequency**3 / self.C_LIGHT_SI**2 / I_nu)
    