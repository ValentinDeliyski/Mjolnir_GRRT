#%%
import os
import numpy as np
import ehtim as eh

Simulation_name = "Wormhole_a_0.9_alpha_1/ngEHT"

# Load the image and the array
Sim_image_230 = eh.image.load_txt('Ehtim_Input_Data/Wormhole_a_0.9_alpha_1/Wormhole_data_for_ehtim_230.csv')
Sim_image_345 = eh.image.load_txt('Ehtim_Input_Data/Wormhole_a_0.9_alpha_1/Wormhole_data_for_ehtim_345.csv')
Tlescope_array_230 = eh.array.load_txt('Ehtim_Input_Data/arrays/ngEHT_230.txt')
Tlescope_array_345 = eh.array.load_txt('Ehtim_Input_Data/arrays/ngEHT_345.txt')

Sim_image_230.rf = 230e9 # [ Hz ]
Sim_image_345.rf = 345e9 # [ Hz ]

# Observe the image at two different frequencies
tint_sec = 120
tadv_sec = 600
tstart_hr = 0
tstop_hr = 24
bw_hz = 2e9
obs230 = Sim_image_230.observe(Tlescope_array_230, tint_sec, tadv_sec, tstart_hr, tstop_hr, bw_hz,
                               sgrscat = False, ampcal = True, phasecal = True, ttype = 'fast')
obs345 = Sim_image_345.observe(Tlescope_array_345, tint_sec, tadv_sec, tstart_hr, tstop_hr, bw_hz,
                               sgrscat = False, ampcal = True, phasecal = True, ttype = 'fast')
obslist = [obs230, obs345]

# Resolution
beamparams230 = obs230.fit_beam() # fitted beam parameters (fwhm_maj, fwhm_min, theta) in radians
res230 = obs230.res() # nominal array resolution, 1/longest baseline
beamparams345 = obs345.fit_beam() # fitted beam parameters (fwhm_maj, fwhm_min, theta) in radians
res345 = obs345.res() # nominal array resolution, 1/longest baseline

# Generate an image prior
npix = 128
fov = 1.2 * Sim_image_230.fovx()
zbl = Sim_image_230.total_flux() # total flux
prior_fwhm = 60 * eh.RADPERUAS   # Gaussian size in microarcssec
emptyprior_230 = eh.image.make_square(obs230, npix, fov)
flatprior_230  = emptyprior_230.add_flat(zbl)
gaussprior_230 = emptyprior_230.add_gauss(zbl, (prior_fwhm, prior_fwhm, 0, 0, 0))

emptyprior_345 = eh.image.make_square(obs345, npix, fov)
flatprior_345  = emptyprior_345.add_flat(zbl)
gaussprior_345 = emptyprior_345.add_gauss(zbl, (prior_fwhm, prior_fwhm, 0, 0, 0))

imgr  = eh.imager.Imager(obslist, gaussprior_230, gaussprior_230, zbl,
                         data_term = {'amp':100, 'cphase':200},
                         reg_term = {'simple':1,'tv2':1,'flux':100,'cm':100},
                         maxit = 1000, ttype = 'fast', stop = 1e-10)
imgr.make_image(mf = True, show_updates = False)
out = imgr.out_last()

out = imgr.out_last()
imgr.init_next = out.blur_circ(res230 * 3 / 4)
imgr.prior_next = imgr.init_next
imgr.dat_term_next = {'amp':100, 'cphase':100 * 0.75}
imgr.reg_term_next ={'simple':1,'tv2':50,'flux':50,'cm':50}
imgr.maxi_next = 3000
imgr.make_image(mf = True, show_updates = False)

#==========================================#
#               Imager run 3               #
#==========================================#

out = imgr.out_last()
imgr.init_next = out.blur_circ(res230 / 2)
imgr.prior_next = imgr.init_next
imgr.dat_term_next = {'amp': 100, 'cphase': 100 * 0.5}
imgr.reg_term_next = {'simple': 1,'tv2': 100,'flux': 10,'cm': 10}
imgr.maxi_next = 4000
imgr.make_image(mf = True, show_updates = False)

#==========================================#
#               Imager run 4               #
#==========================================#

out = imgr.out_last()
imgr.init_next = out.blur_circ(res230 / 3)
imgr.prior_next = imgr.init_next
imgr.dat_term_next = {'amp': 100, 'cphase': 100}
imgr.reg_term_next = {'simple': 1,'tv2': 500,'flux': 1,'cm': 1}
imgr.maxi_next = 4000
imgr.make_image(mf = True, show_updates = False)

#=========================================#
#              Final Outputs              #
#=========================================#

out = imgr.out_last()

#=========== Pre-Blur Results ===========#

out230_mf = out.get_image_mf(230.e9)
out345_mf = out.get_image_mf(345.e9)

print("|=================== Pre-Blur Spectral Index ===================|")

out_specIDX_no_blur = out.copy()
out_specIDX_no_blur.imvec = out.specvec
out_specIDX_no_blur.display(cfun='jet')

#=========== Post-Blur Results ===========#

out230_mf_blur = out230_mf.blur_gauss(beamparams230, 0.5)
out345_mf_blur = out345_mf.blur_gauss(beamparams345, 0.5)

print("|=================== Post-Blur Spectral Index ===================|")

out_blur = out.blur_gauss(beamparams230, 0.5)
out_specIDX_blur = out_blur.copy()
out_specIDX_blur.imvec = out_blur.specvec
out_specIDX_blur.display(cfun = 'jet')

out230_mf_blur.display(cbar_unit = ['Tb'])
out345_mf_blur.display(cbar_unit = ['Tb'])

if not os.path.exists('Ehtim_Output_Data/' + Simulation_name + '/'):
    os.makedirs('Ehtim_Output_Data/' + Simulation_name + '/')

obs230.save_uvfits('Ehtim_Output_Data/' + Simulation_name + '/obs_230.uvp')
obs230.save_txt('Ehtim_Output_Data/' + Simulation_name + '/obs_230.txt')
obs345.save_uvfits('Ehtim_Output_Data/' + Simulation_name + '/obs_345.uvp')
obs345.save_txt('Ehtim_Output_Data/' + Simulation_name + '/obs_345.txt')

out230_mf.save_txt('Ehtim_Output_Data/' + Simulation_name + '/Results_230.txt')
out230_mf.save_fits('Ehtim_Output_Data/' + Simulation_name + '/Results_230.fits')
out345_mf.save_txt('Ehtim_Output_Data/' + Simulation_name + '/Results_345.txt')
out345_mf.save_fits('Ehtim_Output_Data/' + Simulation_name + '/Results_345.fits')

out230_mf_blur.save_txt('Ehtim_Output_Data/' + Simulation_name + '/Results_230_blur.txt')
out230_mf_blur.save_fits('Ehtim_Output_Data/' + Simulation_name + '/Results_230_blur.fits')
out345_mf_blur.save_txt('Ehtim_Output_Data/' + Simulation_name + '/Results_345_blur.txt')
out345_mf_blur.save_fits('Ehtim_Output_Data/' + Simulation_name + '/Results_345_blur.fits')

out_specIDX_no_blur.save_txt('Ehtim_Output_Data/' + Simulation_name + '/Results_specIDX_no_blur.txt')
out_specIDX_no_blur.save_fits('Ehtim_Output_Data/' + Simulation_name + '/Results_specIDX_no_blur.fits')

out_specIDX_blur.save_txt('Ehtim_Output_Data/' + Simulation_name + '/Results_specIDX_blur.txt')
out_specIDX_blur.save_fits('Ehtim_Output_Data/' + Simulation_name + '/Results_specIDX_blur.fits')

obs345_chi_amp = obs345.chisq(out345_mf, ttype = 'fast', dtype = 'amp')
obs230_chi_amp = obs230.chisq(out230_mf, ttype = 'fast', dtype = 'amp')

obs345_chi_cphase = obs345.chisq(out345_mf, ttype = 'fast', dtype = 'cphase')
obs230_chi_cphase = obs230.chisq(out230_mf, ttype = 'fast', dtype = 'cphase')

msg1 = ("chi2 amp,    345 = {} | chi2 amp,    230 = {}".format(np.round(obs345_chi_amp, 5), np.round(obs230_chi_amp, 5)))
msg2 = ("chi2 cphase, 345 = {} | chi2 cphase, 230 = {}".format(np.round(obs345_chi_cphase, 5), np.round(obs230_chi_cphase, 5)))

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
