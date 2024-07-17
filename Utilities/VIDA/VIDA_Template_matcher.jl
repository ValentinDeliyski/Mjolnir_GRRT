read(run(`powershell cls`), String)

using Plots
using VIDA
using InteractiveUtils
using DelimitedFiles

Simulation_name = "Wormhole_ngEHT"
File_name = "Results_345_blur.fits"

display("Parsing EHTIM results...")

Ehtim_Image = load_fits(".\\Ehtim\\Ehtim_Output_Data\\" * Simulation_name * "\\" * File_name);

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

if !ispath("VIDA\\VIDA_Output_Data\\" * Simulation_name * "\\")
    mkpath("VIDA\\VIDA_Output_Data\\" * Simulation_name * "\\")
end

savefig(Final_plot, "VIDA\\VIDA_Output_Data\\" * Simulation_name * "\\VIDA_plot")
display(Final_plot)

writedlm("VIDA\\VIDA_Output_Data\\" * Simulation_name * "\\fit_params.csv", 
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