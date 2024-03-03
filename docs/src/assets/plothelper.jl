using PyPlot
# using PyCall

# need a separate Python helper file for
# centered lower figure as colon operator 
# not implemented fully implemented in PyCall
# gs = pyimport("matplotlib.gridspec")
# pushfirst!(PyVector(pyimport("sys")."path"), joinpath("assets"))
# plothelper = pyimport("plothelper")

function plotmodel(modelvect; ymaxG = nothing, ymaxJ = nothing)

    dt = 0.025
    # Create a time only dataset
    dϵ = timeline(t_start = 0.0, t_end = 10.0, step = dt)
    dσ = timeline(t_start = 0.0, t_end = 10.0, step = dt)

    # calculate a strain/stress data by appling a function of time (by defalut a unit step otherwise substitute hstep(amp = 2.))
    dϵ = strainfunction(dϵ, hstep())
    dσ = stressfunction(dσ, hstep()) - stressfunction(dσ, hstep(offset = 5.0))

    colplot = ["red", "#dc7633", "#229954", "#d4ac0d", "blue"]

    fig, axs = subplots(3,1, figsize=(7,20))
    ax1, ax2, ax3 = axs
    # spec.update(wspace = 0.25,
    #             hspace = 0.25,
    #             top = 0.96,
    #             bottom = 0.07,
    #             left = 0.07,
    #             right = 0.98)

    for i = 1:1:length(modelvect)
        #Relaxation modulus
        dG = modelpredict(dϵ, modelvect[i])
        ax1.plot(dG.t, dG.σ, color=colplot[i])
        ax1.plot([-5*dt,0,dG.t[1]],[0.0, 0.0,dG.σ[1]], color=colplot[i])

        # Creep modulus
        dJ = modelpredict(dσ, modelvect[i,1])
        ax2.plot(dJ.t, dJ.ϵ, color=colplot[i])
        ax2.plot([-5*dt,0,dJ.t[1]],[0.0, 0.0,dJ.ϵ[1]], color=colplot[i])
    end

    # Plot settings

    ax1.set_xlabel("Time", fontsize = 10)
    ax1.set_ylabel("Stress", fontsize = 10)
    ax1.tick_params("both", labelsize = 9)
    ax1.set_title("Relaxation response", fontsize = 10)
    isnothing(ymaxJ) ? ax1.set_ylim(bottom = -dt, top = ymaxG) : ax1.set_ylim(bottom = -dt)
    ax1.tick_params("both", labelsize = 9)
    ax1.tick_params(axis="x", which="minor", bottom=true)
    ax1.grid(alpha = 0.2, axis="y", which="major")
    ax1.grid(alpha = 0.2, axis="x", which="major")

    ax2.set_xlabel("Time", fontsize = 10)
    ax2.set_ylabel("Strain", fontsize = 10)
    ax2.tick_params("both", labelsize = 9)
    ax2.set_title("Creep response", fontsize = 10)
    isnothing(ymaxJ) ? ax2.set_ylim(bottom = -dt, top = ymaxJ) : ax2.set_ylim(bottom = -dt)
    ax2.tick_params("both", labelsize = 9)
    ax2.tick_params(axis="x", which="minor", bottom=true)
    ax2.grid(alpha = 0.2, axis="y", which="major")
    ax2.grid(alpha = 0.2, axis="x", which="major")

    dω = frequencyspec(ω_start = 1e-2, ω_end = 1e2, logstep = log10(1e2/1e-2)/30)

    for i = 1:1:length(modelvect)
        # Storage and Loss moduli
        d_Gp = dynamicmodelpredict(dω, modelvect[i])
        ax3.loglog(d_Gp.ω, d_Gp.Gp, color=colplot[i], "-")
        ax3.loglog(d_Gp.ω, d_Gp.Gpp, color=colplot[i], "--")
    end

    ax3.set_xlabel("Frequency", fontsize = 10)
    ax3.set_ylabel("Storage (—) and Loss (- -) moduli", fontsize = 10)
    ax3.tick_params("both", labelsize = 9)
    ax3.grid(alpha = 0.2, axis="y", which="major")
    ax3.grid(alpha = 0.2, axis="x", which="major")
    ax3.set_title("Frequency response", fontsize = 10)

    return fig

end
