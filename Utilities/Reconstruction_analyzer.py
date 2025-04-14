from csv import reader
from copy import deepcopy
from numpy import average, std, pi, argsort, array
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from sys import path as sys_path
from dataclasses import dataclass
from os import walk, path as os_path

from Support_functions.Parsers import VIDA_params_Parser, ehtim_Parser, Units_class
from Support_functions.Image_processing import get_template_pixel_mask, get_brigness_depression_ratio

""" Add the parent directory of this file to the search path, 
    so this script can be ran from the "Utilities" folder. """
parent_directory = os_path.abspath('...')
sys_path.append(parent_directory)

@dataclass
class Template_parameters_class:

    Simulation_name: list
    Radius: list
    Sigma: list
    Tau: list
    rot_angle: list
    slash: list
    slash_angle: list
    x0: list
    y0: list
    div: list
    
    def __init__(self):
        
        self.Simulation_name = []
        self.Radius = []
        self.Sigma = []
        self.Tau = []
        self.rot_angle = []
        self.slash = []
        self.slash_angle = []
        self.x0 = []
        self.y0 = []
        self.div = []
    
@dataclass        
class Reconstruction_parameters_calss:
    
    Simulation_name: list
    f_measure: list
    flux: list
    chi2_amp: list
    chi2_phase: list
    
    def __init__(self):
        
        self.Simulation_name = []
        self.f_measure = []
        self.flux = []
        self.chi2_amp = []
        self.chi2_phase = []
    
def Scan_trough_VIDA_results(Simulation_case: str, EHT_array: str, Frequency: int = 230) -> list:
    
    Fit_param_file_list = []
    
    """ This scans trough all the VIDA results for a certain metric. """
    for VIDA_results_root_folder, _, VIDA_file_list in walk(parent_directory + Simulation_case):
  
        """ Check that this results folder is for the desired EHT array """
        if EHT_array in VIDA_results_root_folder:
        
            for file in VIDA_file_list:
                
                filename_key = "fit_params"
                
                if EHT_array == "ngEHT":
                
                    filename_key = filename_key + "_{}".format(Frequency)
                        
                if filename_key in file:
                        
                    Fit_param_file_list.append(VIDA_results_root_folder + "\\{}".format(filename_key))    

    return Fit_param_file_list

def Scan_trough_Ehtim_results(Simulation_case: str, EHT_array: str, Frequency: int = 230) -> tuple[list, list]:
    
    Reconstruction_file_list = []
    Chi2_file_list = []
    
    """ This scans trough all the ehtim results for a certain metric. """
    for Ehtim_results_root_folder, _, Ehtim_file_list in walk(parent_directory + Simulation_case):
                
        """ Check that this results folder is for the desired EHT array """
        if EHT_array in Ehtim_results_root_folder:
        
            for file in Ehtim_file_list:
                
                reconstruction_filename_key = "Results_blur"
                chi2_filename_key = "Chi2"
                
                if EHT_array == "ngEHT":
                
                    reconstruction_filename_key = "Results_blur_{}".format(Frequency)
                    chi2_filename_key = "Chi2_{}".format(Frequency)
                        
                if reconstruction_filename_key in file and ".fits" not in file:
                        
                    Reconstruction_file_list.append(Ehtim_results_root_folder + "\\{}".format(reconstruction_filename_key))   
                
                if chi2_filename_key in file:
                        
                    Chi2_file_list.append(Ehtim_results_root_folder + "\\{}".format(chi2_filename_key))   
         
    return Reconstruction_file_list, Chi2_file_list

def Get_reconstruction_parameters(Reconstruction_file_list: list, Chi2_file_list: list, Fit_param_file_list: list) -> tuple[Reconstruction_parameters_calss, Template_parameters_class]:
    
    Reconstruction_parameters = Reconstruction_parameters_calss()
    Template_parameters = Template_parameters_class()
    
    for idx, (Reconstruction_file, Chi2_file, Fit_param_file) in enumerate(zip(Reconstruction_file_list, Chi2_file_list, Fit_param_file_list)):
        
        """ ====================== Parse the ehtim logs ====================== """
        
        Ehtim_Parser = ehtim_Parser(File_name = Reconstruction_file)
        
        Reconstruction_parameters.flux.append(Ehtim_Parser.get_total_flux())
        Reconstruction_parameters.Simulation_name.append(Reconstruction_file.split("\\")[-2])
        
        """ These are used for the template pixel masks, which in turn are used to compute the f measure. """
        
        Ehtim_intensity, _ = Ehtim_Parser.get_plottable_ehtim_data()
        Ehtim_image_FOV = abs(Ehtim_Parser.WINDOW_LIMITS[0] - Ehtim_Parser.WINDOW_LIMITS[1]) * Units_class.MEGA 
        
        """ ====================== Parse the VIDA logs ====================== """
        
        VIDA_parser = VIDA_params_Parser(File_name = Fit_param_file)
        Template_parameters.Radius.append(VIDA_parser.d0 / 2)
        Template_parameters.Sigma.append(VIDA_parser.Sigma)
        Template_parameters.Tau.append(VIDA_parser.Tau)
        
        if VIDA_parser.rot_angle < 0: 
            Template_parameters.rot_angle.append(VIDA_parser.rot_angle + pi)
        else:
            Template_parameters.rot_angle.append(VIDA_parser.rot_angle)
            
        Template_parameters.slash.append(VIDA_parser.slash)
        Template_parameters.slash_angle.append(VIDA_parser.slash_angle)
        Template_parameters.x0.append(VIDA_parser.x0)
        Template_parameters.y0.append(VIDA_parser.y0)
        Template_parameters.div.append(VIDA_parser.div)
        Template_parameters.Simulation_name.append(Fit_param_file.split("\\")[-2])
        
        """ ====================== Comput the f measure ====================== """
        
        ring_mask, dark_spot_mask = get_template_pixel_mask(VIDA_parser, Ehtim_image_FOV, Ehtim_Parser.X_PIXEL_COUNT)
        Reconstruction_parameters.f_measure.append(get_brigness_depression_ratio(ring_mask = ring_mask, dark_spot_mask = dark_spot_mask, Ehtim_intensity = Ehtim_intensity))
        
        """ ====================== Parse the Chi2 logs ====================== """
        
        with open(Chi2_file + ".csv", 'r') as file:
            csvreader = reader(file, delimiter = " ")
            _ = csvreader.__next__()
                
            Reconstruction_parameters.chi2_amp.append(float(csvreader.__next__()[3].split("|")[0]))
            Reconstruction_parameters.chi2_phase.append(float(csvreader.__next__()[3].split("|")[0]))
            
        """ This asserts that the Ehtim and VIDA folders were traversed in the same order and the results match up. """
        try: 
            assert(Fit_param_file.split("\\")[-2] == Reconstruction_file.split("\\")[-2])
        except:
            print("Ehtim and VIDA traversion order mismatch!")
            exit(-1)
        
    return Reconstruction_parameters, Template_parameters

def Process_parameters(Reconstruction_parameters: Reconstruction_parameters_calss, Template_parameters: Template_parameters_class, func: callable) -> tuple[Reconstruction_parameters_calss, Template_parameters_class]:
    
    Processed_reconstruction_parameters = Reconstruction_parameters_calss()
    Processed_reconstruction_parameters.chi2_amp = func(Reconstruction_parameters.chi2_amp)
    Processed_reconstruction_parameters.chi2_phase = func(Reconstruction_parameters.chi2_phase)
    Processed_reconstruction_parameters.flux = func(Reconstruction_parameters.flux)
    Processed_reconstruction_parameters.f_measure = func(Reconstruction_parameters.f_measure)
    
    Processed_template_parameters = Template_parameters_class()
    Processed_template_parameters.Radius = func(Template_parameters.Radius)
    Processed_template_parameters.Sigma = func(Template_parameters.Sigma)
    Processed_template_parameters.Tau = func(Template_parameters.Tau)
    Processed_template_parameters.rot_angle = func(Template_parameters.rot_angle)
    Processed_template_parameters.slash = func(Template_parameters.slash)
    Processed_template_parameters.slash_angle = func(Template_parameters.slash_angle)
    Processed_template_parameters.x0 = func(Template_parameters.x0)
    Processed_template_parameters.y0 = func(Template_parameters.y0)
    Processed_template_parameters.div = func(Template_parameters.div)
    
    return Processed_reconstruction_parameters, Processed_template_parameters

def Compute_quality_metric(Reconstruction_parameters: Reconstruction_parameters_calss, 
                           Avg_reconstruction_parameters: Reconstruction_parameters_calss,
                           Std_reconstruction_parameters: Reconstruction_parameters_calss,
                           Template_parameters: Template_parameters_class,
                           Avg_template_parameters: Template_parameters_class,
                           Std_template_parameters: Template_parameters_class) -> list:
    
    flux_term = []
    for flux in Reconstruction_parameters.flux:
        flux_term.append(((flux - Avg_reconstruction_parameters.flux) / Std_reconstruction_parameters.flux)**2)
    
    f_measure_term = []
    for f_measure in Reconstruction_parameters.f_measure:
        f_measure_term.append(((f_measure - Avg_reconstruction_parameters.f_measure) / Std_reconstruction_parameters.f_measure)**2)
    
    radius_term = []
    for radius in Template_parameters.Radius:
        radius_term.append(((radius - Avg_template_parameters.Radius) / Std_template_parameters.Radius)**2)
    
    sigma_term = []
    for sigma in Template_parameters.Sigma:
        sigma_term.append(((sigma - Avg_template_parameters.Sigma) / Std_template_parameters.Sigma)**2)
    
    slash_term = []
    for slash in Template_parameters.slash:
        slash_term.append(((slash - Avg_template_parameters.slash) / Std_template_parameters.slash)**2)
    
    slash_angle_term = []
    for slash_angle in Template_parameters.slash_angle:
        slash_angle_term.append(((slash_angle - Avg_template_parameters.slash_angle) / Std_template_parameters.slash_angle)**2)
        
    tau_term = []
    for tau in Template_parameters.Tau:
        tau_term.append(((tau - Avg_template_parameters.Tau) / Std_template_parameters.Tau)**2)
        
    rot_angle_term = []
    for rot_angle in Template_parameters.rot_angle:
        rot_angle_term.append(((rot_angle - Avg_template_parameters.rot_angle) / Std_template_parameters.rot_angle)**2)
    
    Quality_metric = []
    
    for idx in range(len(flux_term)):
        Quality_metric.append(flux_term[idx] + f_measure_term[idx] + radius_term[idx] + sigma_term[idx] + slash_term[idx] + slash_angle_term[idx] + tau_term[idx] + rot_angle_term[idx])
        
    return Quality_metric

if __name__ == "__main__":
    
    EHT_array = "ngEHT"

    Reconstruction_file_list, Chi2_file_list = Scan_trough_Ehtim_results(Simulation_case = "Ehtim\\Ehtim_Output_Data\\Sim_Paper_2\\run_2\\M87_Wormhole_a_0_redshift_2", EHT_array = EHT_array, Frequency = 230)
    Fit_param_file_list = Scan_trough_VIDA_results(Simulation_case = "VIDA\\VIDA_Output_Data\\Sim_Paper_2\\run_2\\M87_Wormhole_a_0_redshift_2", EHT_array = EHT_array, Frequency = 230)

    Reconstruction_parameters, Template_parameters = Get_reconstruction_parameters(Reconstruction_file_list, Chi2_file_list, Fit_param_file_list)

    Avg_reconstruction_parameters, Avg_template_parameters = Process_parameters(Reconstruction_parameters = Reconstruction_parameters, Template_parameters = Template_parameters, func = average)
    Std_reconstruction_parameters, Std_template_parameters = Process_parameters(Reconstruction_parameters = Reconstruction_parameters, Template_parameters = Template_parameters, func = std)

    Quality_metric = Compute_quality_metric(Reconstruction_parameters, Avg_reconstruction_parameters, Std_reconstruction_parameters, Template_parameters, Avg_template_parameters, Std_template_parameters)

    Sim_names = deepcopy(Reconstruction_parameters.Simulation_name)

    sorted_idx = argsort(Quality_metric) 
    
    Quality_metric = array(Quality_metric)[sorted_idx]
    Sim_names = array(Sim_names)[sorted_idx]
    
    for Q_metric, Sim_name in zip(Quality_metric, Sim_names):
        print(Q_metric, Sim_name)
    

    params = {"ytick.color" : "black",
                "xtick.color" : "black",
                "axes.labelcolor" : "black",
                "axes.edgecolor" : "black",
                "text.usetex" : True,
                "font.family" : "serif",
                "font.serif" : ["Computer Modern Serif"]}
        
    plt.rcParams.update(params)

    """ ================ Template distribution figure ================ """

    Template_fig = plt.figure(figsize = (8, 6))
    Radius_plot = Template_fig.add_subplot(321)
    Radius_plot.hist(Template_parameters.Radius, bins = 'auto', edgecolor = 'black', align = "mid") 
    Radius_plot.set_xlabel(r"$r_0$ [$\mu$arcsec]")

    Sigma_plot = Template_fig.add_subplot(322)
    Sigma_plot.hist(Template_parameters.Sigma, bins = 'auto', edgecolor = 'black', align = "mid") 
    Sigma_plot.set_xlabel(r"$\sigma$ [$\mu$arcsec]")

    Tau_plot = Template_fig.add_subplot(323)
    Tau_plot.hist(Template_parameters.Tau, bins = 'auto', edgecolor = 'black', align = "mid") 
    Tau_plot.set_xlabel(r"$\tau$ [-]")

    Rot_angle_plot = Template_fig.add_subplot(324)
    Rot_angle_plot.hist(Template_parameters.rot_angle, bins = 'auto', edgecolor = 'black', align = "mid") 
    Rot_angle_plot.set_xlabel(r"$\xi_\tau$ [rad]")

    Slash_plot = Template_fig.add_subplot(325)
    Slash_plot.hist(Template_parameters.slash, bins = 'auto', edgecolor = 'black', align = "mid") 
    Slash_plot.set_xlabel(r"$s$ [-]")

    Slash_angle_plot = Template_fig.add_subplot(326)
    Slash_angle_plot.hist(Template_parameters.slash_angle, bins = 'auto', edgecolor = 'black', align = "mid") 
    Slash_angle_plot.set_xlabel(r"$\xi_s$ [rad]")
    
    Template_fig.tight_layout()
    
    """ ================ Centroid distribution figure ================ """
    
    Centroid_fig = plt.figure(figsize = (8, 8))
    Gridspec = GridSpec(3, 3)

    Main_subplot = Centroid_fig.add_subplot(Gridspec[1:3, :2])
    Main_subplot.invert_xaxis()   

    x_histogram = Centroid_fig.add_subplot(Gridspec[0, :2], sharex = Main_subplot)
    y_histogram = Centroid_fig.add_subplot(Gridspec[1:3, 2], sharey = Main_subplot)
    
    Main_subplot.scatter(Template_parameters.x0, Template_parameters.y0, marker = '.')
    Main_subplot.set(xlabel = r"Centroid X [$\mu$arcsec]", ylabel = r"Centroid Y [$\mu$arcsec]")

    x_histogram.hist(Template_parameters.x0, bins = "auto",align = 'mid', edgecolor = 'black')
    x_histogram.set(ylabel = 'count')

    y_histogram.hist(Template_parameters.y0, bins = "auto", orientation = 'horizontal', align = 'mid', edgecolor = 'black')
    y_histogram.set(xlabel = 'count')
    
    Centroid_fig.tight_layout()
    
    """ ================ Reconstruction distribution figure ================ """

    Reconstruction_fig = plt.figure(figsize = (6, 6))
    f_measure_plot = Reconstruction_fig.add_subplot(221)
    f_measure_plot.hist(Reconstruction_parameters.f_measure, bins = 'auto', edgecolor = 'black', align = "mid") 
    f_measure_plot.set_xlabel(r"$f$ [-]")

    Flux_plot = Reconstruction_fig.add_subplot(222)
    Flux_plot.hist(Reconstruction_parameters.flux, bins = 'auto', edgecolor = 'black', align = "mid") 
    Flux_plot.set_xlabel(r"Flux [Jy]")

    Chi_amp_plot = Reconstruction_fig.add_subplot(223)
    Chi_amp_plot.hist(Reconstruction_parameters.chi2_amp, bins = 'auto', edgecolor = 'black', align = "mid") 
    Chi_amp_plot.set_xlabel(r"$\chi^2_{amp}$ [-]")

    Chi_cphase_plot = Reconstruction_fig.add_subplot(224)
    Chi_cphase_plot.hist(Reconstruction_parameters.chi2_phase, bins = 'auto', edgecolor = 'black', align = "mid") 
    Chi_cphase_plot.set_xlabel(r"$\chi^2_{cphase}$ [-]")

    Reconstruction_fig.tight_layout()

    plt.show()
    