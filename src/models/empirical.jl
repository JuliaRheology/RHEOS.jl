
PowerLawPlateau =  RheoModelClass(
        # Model name
        name="burgerliquid",
        # Model parameters,
        p = [:Gᵩ, :G₀, :τ, :α],
        # Relaxation modulus
        G = quote
                Gᵩ + (G₀ - Gᵩ)./(1 + t/τ).^(α)
            end,
        # Network
        info= "                "
        )
