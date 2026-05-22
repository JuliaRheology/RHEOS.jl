
PowerLawPlateau =  RheoModelClass(
        # Model name
        name="burgerliquid",
        # Model parameters,
        p = (:Gᵩ, :G₀, :τ, :α),
        # Relaxation modulus
        G = quote
                Gᵩ + (G₀ - Gᵩ)./(1 + t/τ).^(α)
            end,
                #TODO: Placeholder eq
        equation = (ϵ =((1.0,1.0),), σ =((1.0,1.0),)),
        # Network
        info= "                "
        )
