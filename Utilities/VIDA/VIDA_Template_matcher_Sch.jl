read(run(`powershell cls`), String)

using Plots
using VIDA
using InteractiveUtils
using DelimitedFiles

# Simulation_name = "Wormhole_a_0_alpha_0\\230GHz\\Sim_Results\\Ehtim_2022"
File_name = "Superposition.fits"
for simple_multiplier in [1, 100, 10]

    for tv_multiplier in [100, 10, 1]
        
        for l1_multiplier in [0, 10, 1, 100]
            
            for fwhm in [50, 60, 40]
                
                for total_flux in [0.5, 0.7, 0.6]
                    
                    for FOV_mult in [1, 1.2]

                        if FOV_mult > 1.1

                            Simulation_name = "Sim_Paper_1_follow_up/run_3/Sch/Sch_ngEHT_simple_" * string(simple_multiplier) * "_tv_" * string(tv_multiplier) * "_l1_" * string(l1_multiplier) * "_fwhm_" * string(fwhm) *  "_tot_flux_" * string(total_flux) *  "_fov_mult_" * string(FOV_mult)

                        else

                            Simulation_name = "Sim_Paper_1_follow_up/run_3/Sch/Sch_ngEHT_simple_" * string(simple_multiplier) * "_tv_" * string(tv_multiplier) * "_l1_" * string(l1_multiplier) * "_fwhm_" * string(fwhm) *  "_tot_flux_" * string(total_flux) *  "_fov_mult_" * string(trunc(Int,FOV_mult))

                        end

                        Ehtim_Image = 0

                        try
                            Ehtim_Image = load_fits("C:\\Users\\Valur\\Documents\\Repos\\Mjolnir_GRRT\\Utilities\\Ehtim\\Ehtim_Output_Data\\" * Simulation_name * "\\" * File_name);
                        catch e
                            continue
                        end

                        if isfile("VIDA_Output_Data\\" * Simulation_name * "\\fit_params_superposition.csv")
                            continue
                        end
                        
                        display("Parsing EHTIM results in " * Simulation_name * "...")
    
                        display("Setting up minimization problem...")
                    
                        bh = Bhattacharyya(Ehtim_Image);
                        
                        lower = GeneralGaussianRing(r0 =  0.4,
                                                    σ  =  0.1,
                                                    τ  =  0.01,
                                                    ξτ = -π,
                                                    s  =  0.01,
                                                    ξs = -π,
                                                    x0 = -80.0,
                                                    y0 = -80.0) 
                    
                        upper = GeneralGaussianRing(r0 = 80.0,
                                                    σ  = 30.0,
                                                    τ  = 0.999,
                                                    ξτ = π,
                                                    s  = 0.999,
                                                    ξs = π,
                                                    x0 = 80.0,
                                                    y0 = 80.0) 
                    
                        initial = GeneralGaussianRing(r0 = 20.0,
                                                    σ  = 5.0,
                                                    τ  = 0.2,
                                                    ξτ = 0.78,
                                                    s  = 0.5,
                                                    ξs = 0.78,
                                                    x0 = 0.0,
                                                    y0 = 0.0)
                    
                        prob = ExtractProblem(bh, initial, lower, upper);
                    
                        display("Running minimizer...")
                    
                        optfilt, divmin = extractor(prob, BBO(maxevals = 70000, tracemode = :silent))
                    
                        Final_plot = plot(triptic(Ehtim_Image, optfilt))
                    
                        display("Saving results...")
                    
                        if !ispath("VIDA_Output_Data\\" * Simulation_name * "\\")
                            mkpath("VIDA_Output_Data\\" * Simulation_name * "\\")
                        end
                    
                        savefig(Final_plot, "VIDA_Output_Data\\" * Simulation_name * "\\VIDA_plot_superposition")
                        display(Final_plot)
                    
                        writedlm("VIDA_Output_Data\\" * Simulation_name * "\\fit_params_superposition.csv", 
                                (optfilt.r0, 
                                optfilt.σ,
                                optfilt.τ,
                                optfilt.ξτ,
                                optfilt.s,
                                optfilt.ξs,
                                optfilt.x0,
                                optfilt.y0, 
                                divmin))
                    
                        display("Finished!")
                    
                    end
                end
            end
        end
    end
end

    
