from Support_functions.Parsers import VIDA_params_Parser, ehtim_Parser, Units_class
from Support_functions.Image_processing import generate_general_gaussian_template, get_template_pixel_mask, get_brigness_depression_ratio
import numpy as np
import matplotlib.pyplot as plt
import csv

file_name_vida = "\\fit_params"

Radius = []
Sigma = []
Tau = []
rot_angle = []
slash = []
slash_angle = []
x0 = []
y0 = []
div = []
f = []
Simulation_name = []
flux = []
chi2 = []
chi2_fail = []
sim_good = []
template = []

fail = []

for simple_multiplier in [1, 10, 100]: #  10, 100

    for tv_multiplier in [1, 10, 100]: #  10, 100
        
        for l1_multiplier in [0, 1, 10, 100]: # [0, 1, 10, 100]
            
            for fwhm in [50, 60, 40]:
                
                for total_flux in [0.5, 0.6, 0.7]:
                    
                    for FOV_mult in [1, 1.2]:
                
                        Simulation_name.append("Sim_Paper_1_follow_up/run_3/Kerr/Kerr_a_0.5_EHT2022_simple_{}_tv_{}_l1_{}_fwhm_{}_tot_flux_{}_fov_mult_{}".format(simple_multiplier, tv_multiplier, l1_multiplier, fwhm, total_flux, FOV_mult))

for sim_name in Simulation_name:
    
    vida_path = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\Utilities\\VIDA\\VIDA_Output_Data\\"
    ehtim_path = "C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\Utilities\\Ehtim\\Ehtim_Output_Data\\"
    
    try:
            
        with open(ehtim_path + sim_name + "\\Chi2.csv", 'r') as file:
            csvreader = csv.reader(file, delimiter = " ")
            _ = csvreader.__next__()
                
            chi2_amp = float(csvreader.__next__()[3].split("|")[0])
            chi2_phse = float(csvreader.__next__()[3].split("|")[0])
            
        if max(chi2_amp, chi2_phse) < 2:

            VIDA_parser = VIDA_params_Parser(vida_path + sim_name + file_name_vida)
            
            EHTIM_parser_230 = ehtim_Parser(ehtim_path + sim_name + "\\Results_blur")
            data_to_plot_ehtim_230,_ = EHTIM_parser_230.get_plottable_ehtim_data()
            # EHTIM_parser_345 = ehtim_Parser(ehtim_path + sim_name + "\\Results_345_blur")
            # data_to_plot_ehtim_345,_ = EHTIM_parser_345.get_plottable_ehtim_data()
            
            chi2.append(max(chi2_amp, chi2_phse))
            sim_good.append(sim_name)
            
            flux.append(EHTIM_parser_230.get_total_flux())
            data_to_plot_ehtim = data_to_plot_ehtim_230

            Units = Units_class()
            
            axes_limits     = np.array([(limit) for limit in EHTIM_parser_230.WINDOW_LIMITS ]) * Units.MEGA
            Ehtim_image_FOV = np.abs(axes_limits[0] - axes_limits[1])
            
            Radius.append(VIDA_parser.template_params["Gaussian_1"]["d0"] / 2)
            Sigma.append(VIDA_parser.template_params["Gaussian_1"]["sigma"])
            Tau.append(VIDA_parser.template_params["Gaussian_1"]["tau"])
            rot_angle.append(VIDA_parser.template_params["Gaussian_1"]["rot_angle"])
            slash.append(VIDA_parser.template_params["Gaussian_1"]["slash"])
            slash_angle.append(VIDA_parser.template_params["Gaussian_1"]["slash_angle"])
            x0.append(VIDA_parser.template_params["Gaussian_1"]["x0"])
            y0.append(VIDA_parser.template_params["Gaussian_1"]["y0"])
            div.append(VIDA_parser.template_params["Gaussian_1"]["div"])
            template = generate_general_gaussian_template(EHTIM_parser_230.X_PIXEL_COUNT, VIDA_parser.template_params["Gaussian_1"], Ehtim_image_FOV)
            
            ring_mask, dark_spot_mask = get_template_pixel_mask(VIDA_parser.template_params, Ehtim_image_FOV, EHTIM_parser_230.X_PIXEL_COUNT)
            f.append(get_brigness_depression_ratio(ring_mask, dark_spot_mask, data_to_plot_ehtim))
            
        else:
            chi2_fail.append(sim_name)  
        
    except:
        
        fail.append(sim_name)
        
print("Failed to find the following simulations:")

for fail1 in fail:
    print(fail1)
    
print("Simulations with chi2 > 2:")

for fail1 in chi2_fail:
    print(fail1)
    
print("\n")
        
Avg_flux = np.average(flux)
Flux_std = np.std(flux)
print("=========================================================================================")
print("Average Flux        = {}, with std = {}, Min = {}, Max = {}".format(round(Avg_flux,4), round(Flux_std,4), round(min(flux),4), round(max(flux),4)))
print("=========================================================================================")
Avg_radius = np.average(Radius)
Radius_std = np.std(Radius)
print("Average Radius      = {}, with std = {}, Min = {}, Max = {}".format(round(Avg_radius,4), round(Radius_std,4), round(min(Radius),4), round(max(Radius),4)))
print("=========================================================================================")
Avg_sigma = np.average(Sigma)
Sigma_std = np.std(Sigma)
print("Average Sigma       = {}, with std = {}, Min = {}, Max = {}".format(round(Avg_sigma, 4), round(Sigma_std, 4), round(min(Sigma),4), round(max(Sigma),4)))
print("=========================================================================================")
Avg_tau = np.average(Tau)
tau_std = np.std(Tau)
print("Average Tau         = {}, with std = {}, Min = {}, Max = {}".format(round(Avg_tau,4), round(tau_std,4), round(min(Tau),4), round(max(Tau),4)))
print("=========================================================================================")
rot_angle = np.array(rot_angle)
rot_angle[rot_angle < 0] = rot_angle[rot_angle < 0] + np.pi
Avg_rot_angle = np.average(rot_angle)
Rot_angle_std = np.std(rot_angle)
print("Average Rot angle   = {}, with std = {}, Min = {}, Max = {}".format(round(Avg_rot_angle,4), round(Rot_angle_std,4), round(min(rot_angle),4), round(max(rot_angle),4)))
print("=========================================================================================")
Avg_slash = np.average(slash)
Slash_std = np.std(slash)
print("Average Slash       = {}, with std = {}, Min = {}, Max = {}".format(round(Avg_slash,4), round(Slash_std,4), round(min(slash),4), round(max(slash),4)))
print("=========================================================================================")
Avg_slash_angle = np.average(slash_angle)
slash_angle_std = np.std(slash_angle)
print("Average Slash angle = {}, with std = {}, Min = {}, Max = {}".format(round(Avg_slash_angle,4), round(slash_angle_std,4),  round(min(slash_angle),4), round(max(slash_angle),4)))
print("=========================================================================================")
Avg_x0 = np.average(x0)
x0_std = np.std(x0)
print("Average x0          = {}, with std = {}, Min = {}, Max = {}".format(round(Avg_x0,4), round(x0_std,4), round(min(x0),4), round(max(x0),4)))
print("=========================================================================================")
Avg_y0 = np.average(y0)
y0_std = np.std(y0)
print("Average y0          = {}, with std = {}, Min = {}, Max = {}".format(round(Avg_y0,4), round(y0_std,4), round(min(y0),4), round(max(y0),4)))
print("=========================================================================================")
Avg_div = np.average(div)
div_std = np.std(div)
print("Average divergence  = {}, with std = {}, Min = {}, Max = {}".format(round(Avg_div,4), round(div_std,4), round(min(div),4), round(max(div),4)))
print("=========================================================================================")
Avg_f = np.average(f)
f_std = np.std(f)
print("Average f           = {}, with std = {}, Min = {}, Max = {}".format(round(Avg_f,4), round(f_std,4), round(min(f),4), round(max(f),4)))
print("=========================================================================================")
Avg_Chi = np.average(chi2)
Chi_std = np.std(chi2)
print("Average Chi2        = {}, with std = {}, Min = {}, Max = {}".format(round(Avg_Chi,4), round(Chi_std,4), round(min(chi2),4), round(max(chi2),4)))
print("=========================================================================================")
print("Number of reconstructions with Chi2 < 2 = {}".format(len(f)))
print("=========================================================================================")

quality_metric = []
centroid_dist = []

for r, sigma, slash1, slash_angle1, tau1, rot_angle1, flux1, f1, x, y in zip(Radius, Sigma, slash, slash_angle, Tau, rot_angle, flux, f, x0, y0):
    
    quality_metric.append(((flux1 - Avg_flux) / Avg_flux)**2 +
                          ((r - Avg_radius) / Avg_radius)**2 + 
                          ((sigma - Avg_sigma) / Avg_sigma)**2 + 
                          ((slash1 - Avg_slash) / Avg_slash)**2 + 
                          ((slash_angle1 - Avg_slash_angle) / Avg_slash_angle)**2 + 
                          ((tau1 - Avg_tau) / Avg_tau)**2 +
                          ((rot_angle1 - Avg_rot_angle) / Avg_rot_angle)**2 +
                          ((f1 - Avg_f) / Avg_f)**2 )
    
    centroid_dist.append(np.sqrt(x**2 + y**2))
    
       
sort_idx = np.array(centroid_dist).argsort()
sorted_centroid_dist= np.array(centroid_dist)[sort_idx]
sorted_metric = np.array(quality_metric)[sort_idx]
sorted_sims = np.array(sim_good)[sort_idx]
sorted_chi2 = np.array(chi2)[sort_idx]
sorted_f = np.array(f)[sort_idx]

print("Sims that are closest to the average:")
for sim, chi, metric, ff, rr in zip(sorted_sims, sorted_chi2, sorted_metric, sorted_f, sorted_centroid_dist):
    print(sim.split("/")[-1], "Chi2 = {}".format(chi), "Quality metric = {}".format(metric), "f = {}, r = {}".format(ff, rr))

params = {"ytick.color" : "black",
              "xtick.color" : "black",
              "axes.labelcolor" : "black",
              "axes.edgecolor" : "black",
              "text.usetex" : True,
              "font.family" : "serif",
              "font.serif" : ["Computer Modern Serif"]}
    
plt.rcParams.update(params)

Fig_1 = plt.figure(1)
Plot_1 = Fig_1.add_subplot(321)
Plot_1.hist(Radius, bins = 'auto', edgecolor = 'black', align = "mid") 
Plot_1.set_xlabel(r"$r_0$ [$\mu$arcsec]")

Plot_2 = Fig_1.add_subplot(322)
Plot_2.hist(Sigma, bins = 'auto', edgecolor = 'black', align = "mid") 
Plot_2.set_xlabel(r"$\sigma$ [$\mu$arcsec]")

Plot_4 = Fig_1.add_subplot(323)
Plot_4.hist(Tau, bins = 'auto', edgecolor = 'black', align = "mid") 
Plot_4.set_xlabel(r"$\tau$ [-]")

Plot_5 = Fig_1.add_subplot(324)
Plot_5.hist(rot_angle, bins = 'auto', edgecolor = 'black', align = "mid") 
Plot_5.set_xlabel(r"$\xi_\tau$ [rad]")

Plot_4 = Fig_1.add_subplot(325)
Plot_4.hist(slash, bins = 'auto', edgecolor = 'black', align = "mid") 
Plot_4.set_xlabel(r"$s$ [-]")

Plot_5 = Fig_1.add_subplot(326)
Plot_5.hist(f, bins = 'auto', edgecolor = 'black', align = "mid") 
Plot_5.set_xlabel(r"$\xi_s$ [rad]")

plt.show()
import matplotlib.gridspec as gridspec

fig = plt.figure(figsize=(8,8))
gs  = gridspec.GridSpec(3, 3)

ax_main = plt.subplot(gs[1:3, :2])
ax_main.invert_xaxis()   

ax_xDist = plt.subplot(gs[0, :2],sharex=ax_main)
ax_yDist = plt.subplot(gs[1:3, 2],sharey=ax_main)
  
ax_main.scatter(x0, y0, marker='.')
ax_main.set(xlabel = r"Centroid X [$\mu$arcsec]", ylabel = r"Centroid Y [$\mu$arcsec]")

ax_xDist.hist(x0, bins = "auto",align = 'mid', edgecolor = 'black')
# ax_xDist.set(ylabel = 'count')

ax_yDist.hist(y0, bins = "auto", orientation = 'horizontal', align = 'mid', edgecolor = 'black')
# ax_yDist.set(xlabel = 'count')

plt.show()
