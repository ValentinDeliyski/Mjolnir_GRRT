read(run(`powershell cls`), String)

using Plots
using VIDA
using InteractiveUtils
using DelimitedFiles

function searchdir(path, key) 
    
    return filter(x->occursin(key, x), readdir(path))

end

function run_VIDA_template_fitter(Simulation_run_root_dir, Simulation, Results_dir_common)

    Ehtim_results_root_folder, Reconstructions, _ = first(walkdir(Simulation_run_root_dir * Simulation))

    for Reconstruction in Reconstructions

        Reconstruction_results = searchdir(Ehtim_results_root_folder * "/" * Reconstruction, "Results_blur")

        for Reconstruction_result in Reconstruction_results

            if endswith(Reconstruction_result, ".fits")
                
                Ehtim_results_path = Ehtim_results_root_folder * "/" * Reconstruction * "/" * Reconstruction_result

                Ehtim_Image = 0
                try
                    Ehtim_Image = load_fits(Ehtim_results_root_folder *  "/" * Reconstruction * "/" * Reconstruction_result)
                catch e
                    flush(stdout)
                    print("Could not parse the EHTIM results in " * Ehtim_results_path * "..." * "\n")
                    continue
                end

                Results_directory = Results_dir_common * "/" * Simulation * "/" * Reconstruction
                Result_file_suffix = split(split(Reconstruction_result, ".fits")[1], "Results_blur")[2]

                if isfile(Results_directory *  "\\fit_params" * Result_file_suffix * ".csv")
                    continue
                end

                flush(stdout)
                print("Parsing EHTIM results in " * Ehtim_results_path * "..." * "\n")
                # display("Setting up minimization problem...")
            
                BH_divergence = Bhattacharyya(Ehtim_Image);
                    
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
            
                prob = ExtractProblem(BH_divergence, initial, lower, upper);
            
                # display("Running minimizer...")
                Optimized_template_params, Optimal_BH_divergence = extractor(prob, BBO(maxevals = 70000, tracemode = :silent))
            
                # Final_plot = plot(triptic(Ehtim_Image, Optimized_template_params))
            
                # display("Saving results...")

                if !ispath(Results_directory)
                    mkpath(Results_directory)
                end
  
                # savefig(Final_plot, Results_directory * "\\VIDA_plot" * Result_file_suffix)
            
                writedlm(Results_directory *  "\\fit_params" * Result_file_suffix * ".csv", 
                        (Optimized_template_params.r0, 
                         Optimized_template_params.σ,
                         Optimized_template_params.τ,
                         Optimized_template_params.ξτ,
                         Optimized_template_params.s,
                         Optimized_template_params.ξs,
                         Optimized_template_params.x0,
                         Optimized_template_params.y0, 
                         Optimal_BH_divergence))
            
                # display("Finished!")
                
            end
    
        end

    end
    
end

Results_dir_common =  @__DIR__() * "/VIDA_Output_Data/Sim_paper_2/run_2"

Simulation_run_root_dir, Sim_folders, _ = first(walkdir("C:/Users/Valur/Documents/Repos/Mjolnir_GRRT/Utilities/Ehtim/Ehtim_Output_Data/Sim_Paper_2/run_2/"))

Threads.@threads for Simulation in Sim_folders

    run_VIDA_template_fitter(Simulation_run_root_dir, Simulation, Results_dir_common)

end




