#%%
import os
import numpy as np
import ehtim as eh

for array in ["EHT2017", "EHT2022"]:

    Simulation_name = "Wormhole_a_0.9_alpha_1/{}".format(array)

    # Load the image and the array
    im  = eh.image.load_txt('Ehtim_Input_Data/Wormhole_a_0.9_alpha_1/Wormhole_data_for_ehtim_230.csv')
    eht = eh.array.load_txt('Ehtim_Input_Data/arrays/{}.txt'.format(array))

    # Look at the image
    im.display(cbar_unit = ['Tb'])

    # Observe the image
    # tint_sec is the integration time in seconds, and tadv_sec is the advance time between scans
    # tstart_hr is the GMST time of the start of the observation and tstop_hr is the GMST time of the end
    # bw_hz is the  bandwidth in Hz
    # sgrscat=True blurs the visibilities with the Sgr A* scattering kernel for the appropriate image frequency
    # ampcal and phasecal determine if gain variations and phase errors are included
    tint_sec = 5
    tadv_sec = 30
    tstart_hr = 0
    tstop_hr = 24
    bw_hz = 4e9
    obs = im.observe(eht, tint_sec, tadv_sec, tstart_hr, tstop_hr, bw_hz,
                    sgrscat = False, ampcal = True, phasecal = True, ttype = 'fast')

    # Resolution
    beamparams = obs.fit_beam() # fitted beam parameters (fwhm_maj, fwhm_min, theta) in radians
    res = obs.res()             # nominal array resolution, 1/longest baseline
    print("Clean beam parameters: ", beamparams)
    print("Nominal Resolution: ", res)

    # Generate an image prior
    npix = 128
    fov = 1.2 * im.fovx()
    zbl = im.total_flux()        # total flux
    prior_fwhm = 60 * eh.RADPERUAS # Gaussian size in microarcssec
    emptyprior = eh.image.make_square(obs, npix, fov)
    flatprior = emptyprior.add_flat(zbl)
    gaussprior = emptyprior.add_gauss(zbl, (prior_fwhm, prior_fwhm, 0, 0, 0))

    # Average the closure quantities and add them to the obsdata object
    avg_time = 600
    obs.add_bispec(avg_time = avg_time)
    obs.add_amp(avg_time = avg_time)
    obs.add_cphase(avg_time = avg_time)
    obs.add_camp(avg_time = avg_time)
    obs.add_logcamp(avg_time = avg_time)

    #==========================================#
    #               Imager run 1               #
    #==========================================#

    flux = zbl
    imgr  = eh.imager.Imager(obs, gaussprior, gaussprior, flux,
                            data_term={'amp':100, 'cphase':200},
                            reg_term={'simple':1,'tv2':1,'flux':100,'cm':100},
                            maxit = 1000, ttype='fast', stop = 1e-10)
    imgr.make_image(show_updates=False)

    #==========================================#
    #               Imager run 2               #
    #==========================================#

    out = imgr.out_last()
    imgr.init_next = out.blur_circ(res * 3 / 4)
    imgr.prior_next = imgr.init_next
    imgr.dat_term_next = {'amp':100, 'cphase':100 * 0.75}
    imgr.reg_term_next = {'simple':1,'tv2':50,'flux':50,'cm':50}
    imgr.maxi_next  = 3000
    imgr.make_image(show_updates=False)

    #==========================================#
    #               Imager run 3               #
    #==========================================#

    out = imgr.out_last()
    imgr.init_next = out.blur_circ(res / 2)
    imgr.prior_next = imgr.init_next
    imgr.dat_term_next = {'amp': 100, 'cphase': 100 * 0.5}
    imgr.reg_term_next = {'simple': 1,'tv2': 100,'flux': 10,'cm': 10}
    imgr.maxi_next  = 4000
    imgr.make_image(show_updates=False)

    #==========================================#
    #               Imager run 4               #
    #==========================================#

    out = imgr.out_last()
    imgr.init_next = out.blur_circ(res / 3)
    imgr.prior_next = imgr.init_next
    imgr.dat_term_next = {'amp': 100, 'cphase': 100}
    imgr.reg_term_next  = {'simple': 1,'tv2': 500,'flux': 1,'cm': 1}
    imgr.maxi_next  = 4000
    imgr.make_image(show_updates=False)

    #=========================================#
    #              Final Outputs              #
    #=========================================#

    out = imgr.out_last()
    out.display(cbar_unit = ['Tb'])

    out_blur = out.blur_gauss(beamparams, 0.5)
    out_blur.display(cbar_unit = ['Tb'])

    if not os.path.exists('Ehtim_Output_Data/' + Simulation_name + '/'):
        os.makedirs('Ehtim_Output_Data/' + Simulation_name + '/')

    # Export the visibility data to uvfits/text
    obs.save_txt('Ehtim_Output_Data/' + Simulation_name + '/obs.txt')    # exports a text file with the visibilities
    obs.save_uvfits('Ehtim_Output_Data/' + Simulation_name + '/obs.uvp') # exports a UVFITS file modeled on template.UVP

    out.save_txt('Ehtim_Output_Data/' + Simulation_name + '/Results.txt')
    out.save_fits('Ehtim_Output_Data/' + Simulation_name + '/Results.fits')

    out_blur.save_txt('Ehtim_Output_Data/' + Simulation_name + '/Results_blur.txt')
    out_blur.save_fits('Ehtim_Output_Data/' + Simulation_name + '/Results_blur.fits')

    obs_chi_amp    = obs.chisq(out, ttype = 'fast', dtype = 'amp')
    obs_chi_cphase = obs.chisq(out, ttype = 'fast', dtype = 'cphase')

    msg1 =  ("chi2 amp = {}".format(np.round(obs_chi_amp, 5)))
    msg2 = ("chi2 cphase = {}".format(np.round(obs_chi_cphase, 5)))

    msg_len = np.maximum(len(msg1), len(msg2))

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

    with open('Ehtim_Output_Data/' + Simulation_name + '/Chi2.csv', 'w') as chi2_file:
            np.savetxt(chi2_file, np.array([msg]), delimiter = " ", fmt = "%s")
# %%
