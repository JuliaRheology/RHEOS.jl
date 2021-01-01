using PyPlot

function plotmodel(modelvect; ymaxG = nothing, ymaxJ = nothing)

    dt = 0.01
    # Create a time only dataset
    dϵ = timeline(t_start = 0.0, t_end = 10.0, step = dt)
    dσ = timeline(t_start = 0.0, t_end = 10.0, step = dt)

    # calculate a strain/stress data by appling a function of time (by defalut a unit step otherwise substitute hstep(amp = 2.))
    dϵ = strainfunction(dϵ, hstep())
    dσ = stressfunction(dσ, hstep()) - stressfunction(dσ, hstep(offset = 5.0))

    colplot = ["red", "#dc7633", "#229954", "#d4ac0d", "blue"]

    fig, ax = subplots(1,2, figsize=(15,7))

    for i = 1:1:length(modelvect)
        #Relaxation modulus
        dG = modelpredict(dϵ, modelvect[i])
        ax[1].plot(dG.t, dG.σ, color=colplot[i])
        ax[1].plot([-5*dt,0,dG.t[1]],[0.0, 0.0,dG.σ[1]], color=colplot[i])

        # Creep modulus
        dJ = modelpredict(dσ, modelvect[i,1])
        ax[2].plot(dJ.t, dJ.ϵ, color=colplot[i])
        ax[2].plot([-5*dt,0,dJ.t[1]],[0.0, 0.0,dJ.ϵ[1]], color=colplot[i])
    end

    # Plot settings

    ax[1].set_xlabel("Time", fontsize = 14);
    ax[1].set_ylabel("Stress", fontsize = 14);
    ax[1].tick_params("both", labelsize=12);
    ax[1].set_title("Relaxation response", fontsize = 14);
    isnothing(ymaxJ) ? ax[1].set_ylim(bottom = -dt, top = ymaxG) : ax[1].set_ylim(bottom = -dt)
    ax[1].tick_params("both", labelsize=12)
    ax[1].tick_params(axis="x", which="minor", bottom=true)
    ax[1].grid(b=true, alpha = 0.2, axis="y", which="major");
    ax[1].grid(b=true, alpha = 0.2, axis="x", which="major");

    ax[2].set_xlabel("Time", fontsize = 14);
    ax[2].set_ylabel("Strain", fontsize = 14);
    ax[2].tick_params("both", labelsize=12);
    ax[2].set_title("Creep response", fontsize = 14);
    isnothing(ymaxJ) ? ax[2].set_ylim(bottom = -dt, top = ymaxJ) : ax[2].set_ylim(bottom = -dt)
    ax[2].tick_params("both", labelsize=12)
    ax[2].tick_params(axis="x", which="minor", bottom=true)
    ax[2].grid(b=true, alpha = 0.2, axis="y", which="major");
    ax[2].grid(b=true, alpha = 0.2, axis="x", which="major");

    dω = frequencyspec()
    fig, ax = subplots(1,1, figsize=(5,5))
    for i = 1:1:length(modelvect)
        # Storage and Loss moduli
        d_Gp = dynamicmodelpredict(dω, modelvect[i])
        ax.loglog(d_Gp.ω, d_Gp.Gp, color=colplot[i], "-")
        ax.loglog(d_Gp.ω, d_Gp.Gpp, color=colplot[i], "--")
    end
    ax.set_xlabel("Frequency", fontsize = 10);
    ax.set_ylabel("Storage (—) and Loss (- -) moduli", fontsize = 10);
    ax.tick_params("both", labelsize=8);
    ax.grid(b=true, alpha = 0.2, axis="y", which="major");
    ax.grid(b=true, alpha = 0.2, axis="x", which="major");
    ax.set_title("Frequency response", fontsize = 10);

end
