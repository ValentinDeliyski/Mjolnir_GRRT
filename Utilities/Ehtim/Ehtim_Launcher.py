
#%%

import os
import sys
import ehtim as eh
from multiprocessing import Pool
from numpy import array, savetxt, maximum

""" Add the parent directory of this file to the search path, 
    so this script can be ran from the "Utilities" folder """
parent_directory = os.path.abspath('...')
sys.path.append(parent_directory)

def run_single_reconstruction(Input_simulation: str, EHT_array: str, Base_reconstruction_params: list, Parameter_modifiers: tuple[float, float, float, float, float, float]) -> None:
                                
    simple_multiplier, tv_multiplier, l1_multiplier, FWHM, Total_flux, FOV_mult = Parameter_modifiers
    
    Simulation_name = ("Sim_Paper_2/run_2/{}/{}".format(Input_simulation, EHT_array) 
                    + "_simple_{}_tv_{}_l1_{}_fwhm_{}_tot_flux_{}_fov_mult_{}".format(simple_multiplier, tv_multiplier, l1_multiplier, FWHM, Total_flux, FOV_mult))
    
    """ This checks weather results for the simulation exist already. If yes, then it skips this reconstruction.
        If not, it then makes a directory for this simulation if there isn't already one. """
    if not os.path.isfile(parent_directory + 'Ehtim/Ehtim_Output_Data/' + Simulation_name + '/Results_blur.txt'):        
        if not os.path.isdir(parent_directory + 'Ehtim/Ehtim_Output_Data/' + Simulation_name + '/'):
            os.makedirs(parent_directory + 'Ehtim/Ehtim_Output_Data/' + Simulation_name + '/')    
    else:
        return    
    
    print("Running " + Simulation_name + "...")

    """ Load the Mjølnir output image and the telescope array. """
    Ray_Tracer_image = eh.image.load_txt(parent_directory + 'Ehtim/Ehtim_Input_Data/{}/Wormhole_data_for_ehtim_230.csv'.format(Input_simulation))
    Telescope_array  = eh.array.load_txt(parent_directory + 'Ehtim/Ehtim_Input_Data/arrays/{}.txt'.format(EHT_array))

    """ Set the simulated observation parameters - these are the example defaults, because I have no idea how to vary them... """
    Integration_time        = 60      # [s]
    Scan_advance_time       = 120    # [s]
    Observation_start_time  = 0      # [hr GMST]
    Observation_stop_time   = 24     # [hr GMST]
    Observation_bandwidth   = 4e9    # [Hz]
    Include_gain_errors     = True   # Enables telescope gain errors
    Include_phase_errors    = True   # Enables telescope phase calibration errors
    Include_SgrA_scattering = True  # Include Sgr A scattering
                            
    Simulated_observation = Ray_Tracer_image.observe(array    = Telescope_array, 
                                                     tint     = Integration_time, 
                                                     tadv     = Scan_advance_time, 
                                                     tstart   = Observation_start_time, 
                                                     tstop    = Observation_stop_time, 
                                                     bw       = Observation_bandwidth,
                                                     sgrscat  = Include_SgrA_scattering, 
                                                     ampcal   = Include_gain_errors, 
                                                     phasecal = Include_phase_errors, 
                                                     ttype    = 'fast')
                            
    """ These are leftovers from the example script - Im not sure how to interpret them... """
    Simulated_beam_parameters = Simulated_observation.fit_beam() # Fitted beam parameters (fwhm_maj, fwhm_min, theta) in [Rad]
    Observation_Resolution    = Simulated_observation.res()      # Nominal array resolution = 1 / longest_baseline
    print("Clean beam parameters: ", Simulated_beam_parameters)
    print("Nominal Resolution: ", Observation_Resolution)

    """ ================================== Generate an initial (prior) image ================================== """
    Linear_pixel_count = 128
    Field_of_view = FOV_mult * Ray_Tracer_image.fovx()
                            
    """ Set the total flux for the initialization image. """
    Init_image_total_flux = Total_flux
                            
    if Init_image_total_flux == None:
        Init_image_total_flux = Ray_Tracer_image.total_flux()
                            
    Initial_image_fwhm = FWHM * eh.RADPERUAS # Gaussian size in microarcssec
                            
    """ Create an empty initialization image, that gets properties added to it in the next function calls. """
    Initial_empty_image = eh.image.make_square(obs  = Simulated_observation, 
                                               npix = Linear_pixel_count, 
                                               fov  = Field_of_view)
                            
    """ Add gaussians ontop of the "flat". Adding the gaussian parameters here so the inputs to the beamaparams variable are more readable. """
    Gauss_prior_x_0   = 0
    Gauss_prior_y_0   = 0
    Gauss_prior_theta = 0
                            
    Initial_flat = Initial_empty_image.add_flat(flux = Init_image_total_flux,
                                                pol  = None)
                            
    Initial_Gaussian = Initial_flat.add_gauss(flux = Init_image_total_flux, 
                                              beamparams = (Initial_image_fwhm, 
                                                            Initial_image_fwhm, 
                                                            Gauss_prior_theta, 
                                                            Gauss_prior_x_0, 
                                                            Gauss_prior_y_0),
                                              pol = None)
    
    """ NOTE: This is not done in the multifrequency example, but done in others and I dont know why... """
    # Averaging_time = 600 # [s]
    # Simulated_observation.add_amp(avg_time = Averaging_time)
    # Simulated_observation.add_cphase(avg_time = Averaging_time)

    #==========================================#
    #               Imager runs                #
    #==========================================#
                            
    Mjolnir_total_flux = Ray_Tracer_image.total_flux() 
        
    for run_idx, (simple_coeff, tv2_coeff, flux_coeff, cm_coff, l1_coeff, amp_coeff, cphase_coeff, max_iter, blur_factor) in enumerate(zip(*Base_reconstruction_params)):

        if run_idx == 0:
                                    
            """ This is the initial imager run. """
            Imager_instance = eh.imager.Imager(obs_in   = Simulated_observation, 
                                               init_im  = Initial_Gaussian, 
                                               prior_im = Initial_Gaussian, 
                                               flux     = Mjolnir_total_flux,
                                               data_term = {'amp':   amp_coeff,
                                                            'cphase': cphase_coeff},
                                               reg_term = {'simple': simple_multiplier * simple_coeff,
                                                           'tv2':    tv_multiplier * tv2_coeff, 
                                                           'flux':   flux_coeff, 
                                                           'cm':     cm_coff, 
                                                           "l1":     l1_multiplier * l1_coeff},
                                               maxit = max_iter, 
                                               ttype ='fast', 
                                               stop  = 1e-10)
                                    
            Imager_instance.make_image(show_updates = False,
                                    pol   = None,
                                    grads = True,
                                    mf    = False)
        else:
            
            Image_output = Imager_instance.out_last()
            
            Imager_instance.init_next  = Image_output.blur_circ(Observation_Resolution * blur_factor)
            Imager_instance.prior_next = Imager_instance.init_next
            
            Imager_instance.dat_term_next = {'amp':   amp_coeff,
                                             'cphase': cphase_coeff}
            
            Imager_instance.reg_term_next = {'simple': simple_multiplier * simple_coeff,
                                             'tv2':    tv_multiplier * tv2_coeff,
                                             'flux':   flux_coeff,
                                             'cm':     cm_coff, 
                                             'l1':     l1_multiplier * l1_coeff}
            
            Imager_instance.maxi_next = max_iter
            Imager_instance.make_image(show_updates = False,
                                       pol   = None,
                                       grads = True,
                                       mf    = False)

    #=========================================#
    #              Final Outputs              #
    #=========================================#

    Final_output = Imager_instance.out_last()
    # Final_output.display(cbar_unit = ['Tb'])

    Final_output_blur = Final_output.blur_gauss(Simulated_beam_parameters, 0.5)
    # Final_output_blur.display(cbar_unit = ['Tb'])

    # Export the visibility data to uvfits/text
    Simulated_observation.save_txt(parent_directory + 'Ehtim/Ehtim_Output_Data/' + Simulation_name + '/obs.txt')    # exports a text file with the visibilities
    Simulated_observation.save_uvfits(parent_directory + 'Ehtim/Ehtim_Output_Data/' + Simulation_name + '/obs.uvp') # exports a UVFITS file modeled on template.UVP

    Final_output.save_txt(parent_directory + 'Ehtim/Ehtim_Output_Data/' + Simulation_name + '/Results.txt')
    Final_output.save_fits(parent_directory + 'Ehtim/Ehtim_Output_Data/' + Simulation_name + '/Results.fits')

    Final_output_blur.save_txt(parent_directory + 'Ehtim/Ehtim_Output_Data/' + Simulation_name + '/Results_blur.txt')
    Final_output_blur.save_fits(parent_directory + 'Ehtim/Ehtim_Output_Data/' + Simulation_name + '/Results_blur.fits')

    obs_chi_amp    = Simulated_observation.chisq(Final_output, ttype = 'fast', dtype = 'camp')
    obs_chi_cphase = Simulated_observation.chisq(Final_output, ttype = 'fast', dtype = 'cphase')

    msg1 = ("chi2 camp = {}".format(round(obs_chi_amp, 5)))
    msg2 = ("chi2 cphase = {}".format(round(obs_chi_cphase, 5)))

    msg_len = maximum(len(msg1), len(msg2))

    msg = "|" + "=" * msg_len + "|" + "\n"

    if (len(msg1) < msg_len):
        msg = msg + "|" + msg1 + " " * (msg_len - len(msg1)) + "|"
    else:
        msg = msg + "|" + msg1 + "|"

    if (len(msg2) < msg_len):
        msg = msg  + "\n" + "|" + msg2 + " " * (msg_len - len(msg2)) + "|"
    else:
        msg = msg  + "\n" + "|" + msg2 + "|"

    msg = msg + "\n" + "|" + "=" * msg_len + "|"

    with open(parent_directory + 'Ehtim/Ehtim_Output_Data/' + Simulation_name + '/Chi2.csv', 'w') as chi2_file:
             savetxt(chi2_file, array([msg]), delimiter = " ", fmt = "%s")
  
def run_multifrequncy_reconstruction(Input_simulation: str, EHT_array: str, Base_reconstruction_params: list, Parameter_modifiers: tuple[float, float, float, float, float, float]) -> None:
    
    simple_multiplier, tv_multiplier, l1_multiplier, FWHM, Total_flux, FOV_mult = Parameter_modifiers
    
    Simulation_name = ("Sim_Paper_2/run_2/{}/{}".format(Input_simulation, EHT_array) 
                      + "_simple_{}_tv_{}_l1_{}_fwhm_{}_tot_flux_{}_fov_mult_{}".format(simple_multiplier, tv_multiplier, l1_multiplier, FWHM, Total_flux, FOV_mult))
    
    """ This checks weather results for the simulation exist already. If yes, then it skips this reconstruction.
        If not, it then makes a directory for this simulation if there isn't already one. """
    if not os.path.isfile(parent_directory + 'Ehtim/Ehtim_Output_Data/' + Simulation_name + '/Results_blur_230.txt'):        
        if not os.path.isdir(parent_directory + 'Ehtim/Ehtim_Output_Data/' + Simulation_name + '/'):
            os.makedirs(parent_directory + 'Ehtim/Ehtim_Output_Data/' + Simulation_name + '/')    
    else:
        return    
    
    print("Running " + Simulation_name + "...")

    """ Load the Mjølnir output image and the telescope array. """
    Ray_Tracer_image_230_GHz = eh.image.load_txt(parent_directory + 'Ehtim/Ehtim_Input_Data/{}/Wormhole_data_for_ehtim_230.csv'.format(Input_simulation))
    Ray_Tracer_image_345_GHz = eh.image.load_txt(parent_directory + 'Ehtim/Ehtim_Input_Data/{}/Wormhole_data_for_ehtim_345.csv'.format(Input_simulation))
    
    """ Set the observation frequency for the two images that Mjolnir spits out. """
    Ray_Tracer_image_230_GHz.rf = 230e9 # [ Hz ]
    Ray_Tracer_image_345_GHz.rf = 345e9 # [ Hz ]
    
    """ Load the ngEHT array parameters for each observing frequency. """
    Telescope_array_230_GHz  = eh.array.load_txt(parent_directory + 'Ehtim/Ehtim_Input_Data/arrays/{}_230.txt'.format(EHT_array))
    Telescope_array_345_GHz  = eh.array.load_txt(parent_directory + 'Ehtim/Ehtim_Input_Data/arrays/{}_345.txt'.format(EHT_array))
    
    """ Set the simulated observation parameters - these are the example defaults, because I have no idea how to vary them... """
    Integration_time        = 120    # [s]
    Scan_advance_time       = 600    # [s]
    Observation_start_time  = 0      # [hr GMST]
    Observation_stop_time   = 24     # [hr GMST]
    Observation_bandwidth   = 2e9    # [Hz]
    Include_gain_errors     = True   # Enables telescope gain errors
    Include_phase_errors    = True   # Enables telescope phase calibration errors
    Include_SgrA_scattering = True  # Include Sgr A scattering
                            
    Simulated_observation_230_GHz = Ray_Tracer_image_230_GHz.observe(array    = Telescope_array_230_GHz, 
                                                                     tint     = Integration_time, 
                                                                     tadv     = Scan_advance_time, 
                                                                     tstart   = Observation_start_time, 
                                                                     tstop    = Observation_stop_time, 
                                                                     bw       = Observation_bandwidth,
                                                                     sgrscat  = Include_SgrA_scattering, 
                                                                     ampcal   = Include_gain_errors, 
                                                                     phasecal = Include_phase_errors, 
                                                                     ttype    = 'fast')
    
    Simulated_observation_345_GHz = Ray_Tracer_image_345_GHz.observe(array    = Telescope_array_345_GHz, 
                                                                     tint     = Integration_time, 
                                                                     tadv     = Scan_advance_time, 
                                                                     tstart   = Observation_start_time, 
                                                                     tstop    = Observation_stop_time, 
                                                                     bw       = Observation_bandwidth,
                                                                     sgrscat  = Include_SgrA_scattering, 
                                                                     ampcal   = Include_gain_errors, 
                                                                     phasecal = Include_phase_errors, 
                                                                     ttype    = 'fast')
    
    Observation_list = [Simulated_observation_230_GHz, Simulated_observation_345_GHz]
                           
    """ These are leftovers from the example script - Im not sure how to interpret them... """
    Simulated_beam_parameters_230_GHz = Simulated_observation_230_GHz.fit_beam() # Fitted beam parameters (fwhm_maj, fwhm_min, theta) in [Rad]
    Observation_Resolution_230_GHz    = Simulated_observation_230_GHz.res()      # Nominal array resolution = 1 / longest_baseline

    Simulated_beam_parameters_345_GHz = Simulated_observation_345_GHz.fit_beam() # Fitted beam parameters (fwhm_maj, fwhm_min, theta) in [Rad]
    Observation_Resolution_345_GHz    = Simulated_observation_345_GHz.res()      # Nominal array resolution = 1 / longest_baseline

    """ ================================== Generate an initial (prior) image ================================== """                        
    Linear_pixel_count = 128
    Field_of_view = FOV_mult * Ray_Tracer_image_230_GHz.fovx()
    
    """ Set the total flux for both initialization images (even though the 345 GHz one ends up not getting used). """
    Init_image_total_flux_230_GHz = Total_flux
    Init_image_total_flux_345_GHz = Total_flux
                            
    if Init_image_total_flux_230_GHz == None:
        Init_image_total_flux_230_GHz = Ray_Tracer_image_230_GHz.total_flux()
        Init_image_total_flux_345_GHz = Ray_Tracer_image_345_GHz.total_flux()
                            
    Initial_image_fwhm = FWHM * eh.RADPERUAS # Gaussian size in microarcssec
    
    """ Create an empty initialization image, that gets properties added to it in the next function calls. """
    Initial_empty_image_230_GHz = eh.image.make_square(obs  = Simulated_observation_230_GHz, 
                                                       npix = Linear_pixel_count, 
                                                       fov  = Field_of_view)
    
    Initial_empty_image_345_GHz = eh.image.make_square(obs  = Simulated_observation_345_GHz, 
                                                       npix = Linear_pixel_count, 
                                                       fov  = Field_of_view)

                            
    """ Add an initial "flat" field to the initialization images, with a specified total flux. """
    Initial_flat_230_GHz = Initial_empty_image_230_GHz.add_flat(flux = Init_image_total_flux_230_GHz,
                                                                pol  = None)
    
    Initial_flat_345_GHz = Initial_empty_image_345_GHz.add_flat(flux = Init_image_total_flux_345_GHz,
                                                                pol  = None)
                                            
    """ Add gaussians ontop of the "flat". Adding the gaussian parameters here so the inputs to the beamaparams variable are more readable. """
    Gauss_prior_x_0   = 0
    Gauss_prior_y_0   = 0
    Gauss_prior_theta = 0
                
    Initial_Gaussian_230_GHz = Initial_flat_230_GHz.add_gauss(flux = Init_image_total_flux_230_GHz, 
                                                              beamparams = (Initial_image_fwhm, 
                                                                            Initial_image_fwhm, 
                                                                            Gauss_prior_theta, 
                                                                            Gauss_prior_x_0, 
                                                                            Gauss_prior_y_0),
                                                              pol = None)

    Initial_Gaussian_345_GHz = Initial_flat_345_GHz.add_gauss(flux = Init_image_total_flux_345_GHz, 
                                                              beamparams = (Initial_image_fwhm, 
                                                                            Initial_image_fwhm, 
                                                                            Gauss_prior_theta, 
                                                                            Gauss_prior_x_0, 
                                                                            Gauss_prior_y_0),
                                                              pol = None)
    
    """ NOTE: This is not done in the multifrequency example, but done in others and I dont know why... """
    # Averaging_time = 600 # [s]
    # Simulated_observation_230_GHz.add_amp(avg_time = Averaging_time)
    # Simulated_observation_230_GHz.add_cphase(avg_time = Averaging_time)
    
    # Simulated_observation_345_GHz.add_amp(avg_time = Averaging_time)
    # Simulated_observation_345_GHz.add_cphase(avg_time = Averaging_time)

    #==========================================#
    #               Imager runs                #
    #==========================================#

    for run_idx, (simple_coeff, tv2_coeff, flux_coeff, cm_coff, l1_coeff, amp_coeff, cphase_coeff, max_iter, blur_factor) in enumerate(zip(*Base_reconstruction_params)):

        if run_idx == 0:
                                    
            """ This is the initial imager run. """
            Imager_instance = eh.imager.Imager(obs_in   = Observation_list, 
                                               init_im  = Initial_Gaussian_230_GHz, 
                                               prior_im = Initial_Gaussian_230_GHz, 
                                               flux     = Init_image_total_flux_230_GHz,
                                               data_term = {'amp':   amp_coeff,
                                                            'cphase': cphase_coeff},
                                               reg_term = {'simple': simple_multiplier * simple_coeff,
                                                           'tv2':    tv_multiplier * tv2_coeff, 
                                                           'flux':   flux_coeff, 
                                                           'cm':     cm_coff, 
                                                           "l1":     l1_multiplier * l1_coeff},
                                               maxit = max_iter, 
                                               ttype = 'fast', 
                                               stop  = 1e-10)
                                    
            Imager_instance.make_image(show_updates = False,
                                       pol   = None,
                                       grads = True,
                                       mf    = True)
        else:
            
            Image_output = Imager_instance.out_last()
            
            Imager_instance.init_next  = Image_output.blur_circ(Observation_Resolution_230_GHz * blur_factor)
            Imager_instance.prior_next = Imager_instance.init_next
            
            Imager_instance.dat_term_next = {'amp':   amp_coeff,
                                             'cphase': cphase_coeff}
            
            Imager_instance.reg_term_next = {'simple': simple_multiplier * simple_coeff,
                                             'tv2':    tv_multiplier * tv2_coeff,
                                             'flux':   flux_coeff,
                                             'cm':     cm_coff, 
                                             'l1':     l1_multiplier * l1_coeff}
            
            Imager_instance.maxi_next = max_iter
            Imager_instance.make_image(show_updates = False,
                                       pol   = None,
                                       grads = True,
                                       mf    = True)

    #=========================================#
    #              Final Outputs              #
    #=========================================#

    Final_output = Imager_instance.out_last()
    
    """ ====================================== Pre-Blur Results ====================================== """
    
    Final_output_230_GHz = Final_output.get_image_mf(Ray_Tracer_image_230_GHz.rf)
    Final_output_345_GHz = Final_output.get_image_mf(Ray_Tracer_image_345_GHz.rf)
    
    Final_spectral_index = Final_output.copy()
    Final_spectral_index.imvec = Final_output.specvec
    
    """ ====================================== Post-Blur Results ===================================== """
    
    Final_output_230_GHz_blur = Final_output_230_GHz.blur_gauss(Simulated_beam_parameters_230_GHz, 0.5)
    Final_output_345_GHz_blur = Final_output_345_GHz.blur_gauss(Simulated_beam_parameters_345_GHz, 0.5)
    
    Final_output_blur = Final_output.blur_gauss(Simulated_beam_parameters_230_GHz, 0.5)
    Final_spectral_index_blur = Final_output_blur.copy()
    Final_spectral_index_blur.imvec = Final_output_blur.specvec

    Simulated_observation_230_GHz.save_uvfits(parent_directory + 'Ehtim/Ehtim_Output_Data/' + Simulation_name + '/obs_230.uvp')
    Simulated_observation_230_GHz.save_txt(parent_directory + 'Ehtim/Ehtim_Output_Data/' + Simulation_name + '/obs_230.txt')
    Simulated_observation_345_GHz.save_uvfits(parent_directory + 'Ehtim/Ehtim_Output_Data/' + Simulation_name + '/obs_345.uvp')
    Simulated_observation_345_GHz.save_txt(parent_directory + 'Ehtim/Ehtim_Output_Data/' + Simulation_name + '/obs_345.txt')

    Final_output_230_GHz.save_txt(parent_directory + 'Ehtim/Ehtim_Output_Data/' + Simulation_name + '/Results_230.txt')
    Final_output_230_GHz.save_fits(parent_directory + 'Ehtim/Ehtim_Output_Data/' + Simulation_name + '/Results_230.fits')
    Final_output_345_GHz.save_txt(parent_directory + 'Ehtim/Ehtim_Output_Data/' + Simulation_name + '/Results_345.txt')
    Final_output_345_GHz.save_fits(parent_directory + 'Ehtim/Ehtim_Output_Data/' + Simulation_name +'/Results_345.fits')

    Final_output_230_GHz_blur.save_txt(parent_directory + 'Ehtim/Ehtim_Output_Data/' + Simulation_name + '/Results_blur_230.txt')
    Final_output_230_GHz_blur.save_fits(parent_directory + 'Ehtim/Ehtim_Output_Data/' + Simulation_name + '/Results_blur_230.fits')
    Final_output_345_GHz_blur.save_txt(parent_directory + 'Ehtim/Ehtim_Output_Data/' + Simulation_name + '/Results_blur_345.txt')
    Final_output_345_GHz_blur.save_fits(parent_directory + 'Ehtim/Ehtim_Output_Data/' + Simulation_name + '/Results_blur_345.fits')

    Final_spectral_index.save_txt(parent_directory + 'Ehtim/Ehtim_Output_Data/' + Simulation_name + '/Results_specIDX_no_blur.txt')
    Final_spectral_index.save_fits(parent_directory + 'Ehtim/Ehtim_Output_Data/' + Simulation_name + '/Results_specIDX_no_blur.fits')

    Final_spectral_index_blur.save_txt(parent_directory + 'Ehtim/Ehtim_Output_Data/' + Simulation_name + '/Results_specIDX_blur.txt')
    Final_spectral_index_blur.save_fits(parent_directory + 'Ehtim/Ehtim_Output_Data/' + Simulation_name + '/Results_specIDX_blur.fits')

    obs345_chi_amp = Simulated_observation_345_GHz.chisq(Final_output_345_GHz, ttype = 'fast', dtype = 'camp')
    obs230_chi_amp = Simulated_observation_230_GHz.chisq(Final_output_230_GHz, ttype = 'fast', dtype = 'camp')

    obs345_chi_cphase = Simulated_observation_345_GHz.chisq(Final_output_345_GHz, ttype = 'fast', dtype = 'cphase')
    obs230_chi_cphase = Simulated_observation_230_GHz.chisq(Final_output_230_GHz, ttype = 'fast', dtype = 'cphase')

    for frequency_idx, (chi2_amp, chi2_cphase) in enumerate(zip([obs230_chi_amp, obs345_chi_amp], [obs230_chi_cphase, obs345_chi_cphase])):

        msg1 = ("chi2 camp = {}".format(round(chi2_amp, 5)))
        msg2 = ("chi2 cphase = {}".format(round(chi2_cphase, 5)))

        msg_len = maximum(len(msg1), len(msg2))

        msg = "|" + "=" * msg_len + "|" + "\n"

        if (len(msg1) < msg_len):
            msg = msg + "|" + msg1 + " " * (msg_len - len(msg1)) + "|"
        else:
            msg = msg + "|" + msg1 + "|"

        if (len(msg2) < msg_len):
            msg = msg  + "\n" + "|" + msg2 + " " * (msg_len - len(msg2)) + "|"
        else:
            msg = msg  + "\n" + "|" + msg2 + "|"

        msg = msg + "\n" + "|" + "=" * msg_len + "|"
        
        match frequency_idx:
    
            case 0:
                filename = "Chi2_230.csv"
            
            case _:
                filename = "Chi2_345.csv"

        with open(parent_directory + 'Ehtim/Ehtim_Output_Data/' + Simulation_name + '/' + filename, 'w') as chi2_file:
                savetxt(chi2_file, array([msg]), delimiter = " ", fmt = "%s")
  
  
    return
  
def Run_reconstruction_sweep(Input_simulatiaon, EHT_array, Base_reconstruction_params) -> None:

    for simple_multiplier in [1, 10, 100]:

        for tv2_multiplier in [1, 10, 100]:
            
            for l1_multiplier in [0, 1, 10, 100]:
                
                for FWHM in [80]:
                    
                    for Total_flux in [None]:
                        
                        for FOV_mult in [1, 1.2]:
                            
                            Parameter_modifiers = (simple_multiplier, tv2_multiplier, l1_multiplier, FWHM, Total_flux, FOV_mult)

                            if EHT_array == "ngEHT":
                                
                                run_multifrequncy_reconstruction(Input_simulation = Input_simulatiaon,
                                                                 EHT_array = EHT_array,
                                                                 Base_reconstruction_params = Base_reconstruction_params,
                                                                 Parameter_modifiers = Parameter_modifiers)
                            else:
                                
                                run_single_reconstruction(Input_simulation = Input_simulatiaon,
                                                          EHT_array = EHT_array,
                                                          Base_reconstruction_params = Base_reconstruction_params,
                                                          Parameter_modifiers = Parameter_modifiers)

if __name__ == "__main__":
    
    Simple_coeffs = [1,   1,  1,   1]
    TV2_coeffs    = [1,   50, 100, 500]
    Flux_coeffs   = [100, 50, 10,  1]
    CM_coeffs     = [100, 50, 10,  1]
    L1_coeffs     = [1,   1,  1,   1]

    Amp_coeffs    = [100, 100, 100, 100]
    Cphase_coeffs = [200, 75, 50, 100]

    Max_iterations = [1000, 3000, 4000, 4000]
    Gauss_blur_factor = [1, 0.75, 0.5, 0.33]
    
    Base_reconstruction_params = (Simple_coeffs, TV2_coeffs, Flux_coeffs, CM_coeffs, L1_coeffs, Amp_coeffs, Cphase_coeffs, Max_iterations, Gauss_blur_factor)
       
    # Input_simulations = [["Sgr_A_Wormhole_a_0.5_redshift_0", "EHT2017"],
    #                      ["Sgr_A_Wormhole_a_0.5_redshift_1", "EHT2017"],
    #                      ["Sgr_A_Wormhole_a_0.5_redshift_2", "EHT2017"],
    #                      ["Sgr_A_Wormhole_a_0.9_redshift_0", "EHT2017"],
    #                      ["Sgr_A_Wormhole_a_0.9_redshift_1", "EHT2017"],
    #                      ["Sgr_A_Wormhole_a_0.9_redshift_2", "EHT2017"],
    #                      ["Sgr_A_Wormhole_a_0_redshift_0", "EHT2017"],
    #                      ["Sgr_A_Wormhole_a_0_redshift_1", "EHT2017"],
    #                      ["Sgr_A_Wormhole_a_0_redshift_2", "EHT2017"]]
                         
    Input_simulations = [["Sgr_A_Wormhole_a_0.5_redshift_0", "EHT2022"],
                         ["Sgr_A_Wormhole_a_0.5_redshift_1", "EHT2022"],
                         ["Sgr_A_Wormhole_a_0.5_redshift_2", "EHT2022"],
                         ["Sgr_A_Wormhole_a_0.9_redshift_0", "EHT2022"],
                         ["Sgr_A_Wormhole_a_0.9_redshift_1", "EHT2022"],
                         ["Sgr_A_Wormhole_a_0.9_redshift_2", "EHT2022"],
                         ["Sgr_A_Wormhole_a_0_redshift_0", "EHT2022"],
                         ["Sgr_A_Wormhole_a_0_redshift_1", "EHT2022"],
                         ["Sgr_A_Wormhole_a_0_redshift_2", "EHT2022"],
                         
                         ["Sgr_A_Wormhole_a_0.5_redshift_0", "EHT2017"],
                         ["Sgr_A_Wormhole_a_0.5_redshift_1", "EHT2017"],
                         ["Sgr_A_Wormhole_a_0.5_redshift_2", "EHT2017"],
                         ["Sgr_A_Wormhole_a_0.9_redshift_0", "EHT2017"],
                         ["Sgr_A_Wormhole_a_0.9_redshift_1", "EHT2017"],
                         ["Sgr_A_Wormhole_a_0.9_redshift_2", "EHT2017"],
                         ["Sgr_A_Wormhole_a_0_redshift_0", "EHT2017"],
                         ["Sgr_A_Wormhole_a_0_redshift_1", "EHT2017"],
                         ["Sgr_A_Wormhole_a_0_redshift_2", "EHT2017"]]

    for Input_simulation in Input_simulations:
        Input_simulation.append(Base_reconstruction_params) 
    
    with Pool(10) as pool:
        pool.starmap(Run_reconstruction_sweep, Input_simulations)
    
    
    # Run_reconstruction_sweep("Sgr_A_Wormhole_a_0.5_redshift_0", "EHT2017", Base_reconstruction_params)

    