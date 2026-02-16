using Plots
# using PyCall

# need a separate Python helper file for
# centered lower figure as colon operator 
# not implemented fully implemented in PyCall
# gs = pyimport("matplotlib.gridspec")
# pushfirst!(PyVector(pyimport("sys")."path"), joinpath("assets"))
# plothelper = pyimport("plothelper")


function plotmodel(modelvect::Vector{RheoModel}; ymaxG=nothing, ymaxJ=nothing)

    dt = 0.025
    colplot = ["red", "#dc7633", "#229954", "#d4ac0d", "blue"]

    # 1️⃣ Create time-only datasets
    tstrain = timeline(t_start=0.0, t_end=10.0, step=dt)
    tstress = timeline(t_start=0.0, t_end=10.0, step=dt)

    # 2️⃣ Apply default strain/stress functions
    dϵ = strainfunction(tstrain, hstep())
    dσ = stressfunction(tstress, hstep()) - stressfunction(tstress, hstep(offset=5.0))

    # 3️⃣ Initialize Plots layout
    plt = plot(layout=(3,1), size=(700,900))

    # --- Relaxation & Creep ---
    for i in 1:length(modelvect)
        # Relaxation modulus
        dG = modelpredict(dϵ, modelvect[i])
        plot!(plt[1], dG.t, dG.σ, color=colplot[i], label="model $i", xgrid=true, ygrid=true)
        plot!(plt[1], [-5*dt,0,dG.t[1]], [0.0,0.0,dG.σ[1]], color=colplot[i], label="")

        # Creep modulus
        dJ = modelpredict(dσ, modelvect[i,1])
        plot!(plt[2], dJ.t, dJ.ϵ, color=colplot[i], label="model $i", xgrid=true, ygrid=true)
        plot!(plt[2], [-5*dt,0,dJ.t[1]], [0.0,0.0,dJ.ϵ[1]], color=colplot[i], label="")
    end

    # Set y-limits safely
    if !isnothing(ymaxG)
        ylims!(plt[1], -dt, ymaxG)
    end
    if !isnothing(ymaxJ)
        ylims!(plt[2], -dt, ymaxJ)
    end

    xlabel!(plt[1], "Time")
    ylabel!(plt[1], "Stress")
    xlabel!(plt[2], "Time")
    ylabel!(plt[2], "Strain")
    title!(plt[1], "Relaxation response")
    title!(plt[2], "Creep response")

    # --- Frequency response ---
    dω = frequencyspec(ω_start=1e-2, ω_end=1e2, logstep=log10(1e2/1e-2)/30)
    for i in 1:length(modelvect)
        dGp = dynamicmodelpredict(dω, modelvect[i])
        plot!(plt[3], dGp.ω, dGp.Gp, color=colplot[i], label="Gp model $i",
              xscale=:log10, yscale=:log10, xgrid=true, ygrid=true)
        plot!(plt[3], dGp.ω, dGp.Gpp, color=colplot[i], linestyle=:dash, label="Gpp model $i",
              xscale=:log10, yscale=:log10)
    end
    xlabel!(plt[3], "Frequency")
    ylabel!(plt[3], "Storage (—) and Loss (- -) moduli")
    title!(plt[3], "Frequency response")

    return plt
end

