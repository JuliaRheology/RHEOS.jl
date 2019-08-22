#!/usr/bin/env julia

Fract_Maxwell = RheoModelClass(
          # Model name
          name="fractmaxwell",
          # Model parameters,
          p = [:cₐ, :a, :cᵦ, :β],
          # Relaxation modulus
          G = quote
                 cᵦ*t^(-β)*mittleff(a - β, 1 - β, -cᵦ*t^(a - β)/cₐ)
              end,
          # Creep modulus
          J = quote
                t^(a)/(cₐ*gamma(1 + a)) + t^(β)/(cᵦ*gamma(1 + β))
              end,
          # Storage modulus
          Gp = quote
                 denominator = (cₐ*ω^a)^2 + (cᵦ*ω^β)^2 + 2*(cₐ*ω^a)*(cᵦ*ω^β)*cos((a-β)*π/2)
                 numerator = ((cᵦ*ω^β)^2)*(cₐ*ω^a)*cos(a*π/2) + ((cₐ*ω^a)^2)*(cᵦ*ω^β)*cos(β*π/2)
                 numerator/denominator
               end,
         # Loss modulus
          Gpp = quote
                  denominator = (cₐ*ω^a)^2 + (cᵦ*ω^β)^2 + 2*(cₐ*ω^a)*(cᵦ*ω^β)*cos((a-β)*π/2)
                  numerator = ((cᵦ*ω^β)^2)*(cₐ*ω^a)*sin(a*π/2) + ((cₐ*ω^a)^2)*(cᵦ*ω^β)*sin(β*π/2)
                  numerator/denominator
                end,
          # Constraints
          constraint = quote
                   all([ (a<1) & (a>0)
                         (β<1) & (β>0)
                         -a+β < 0     ])
                  end,
          # Network
          info= "
             ___╱╲__________╱╲____
                ╲╱          ╲╱
                  cₐ,a         cᵦ, β
                 "
          )

FractS_Maxwell = RheoModelClass(
        # Model name
        name="fractmaxwell_spring",
        # Model parameters,
        p = [:cₐ, :a, :k],
        # Relaxation modulus
        G = quote
                # special case when α==0.0, 2 springs in series
                a==0.0 ? 1.0/(1.0/cₐ + 1.0/k) : k*mittleff(a, -k*t^a/cₐ) # normal case otherwise
            end,
        # Creep modulus
        J = quote
              t^(a)/(cₐ*gamma(1 + a)) + 1/k
            end,
        # Storage modulus
        Gp = quote
                denominator = (cₐ*ω^a)^2 + k^2 + 2*(cₐ*ω^a)*k*cos(a*π/2)
                numerator = k^2*(cₐ*ω^a)*cos(a*π/2) + (cₐ*ω^a)^2*k
                numerator/denominator
             end,
       # Loss modulus
        Gpp = quote
                denominator = (cₐ*ω^a)^2 + k^2 + 2*(cₐ*ω^a)*k*cos(a*π/2)
                numerator = k^2*(cₐ*ω^a)*sin(a*π/2)
                numerator/denominator
              end,
        # Constraints
        constraint = quote
                 (a<1) & (a>0)
                end,
        # Network
        info= "
           ___╱╲_________╱╲  ╱╲  ╱╲  ________
              ╲╱           ╲╱  ╲╱  ╲╱
                cₐ,a               k
               "
        )


FractD_Maxwell = RheoModelClass(
          # Model name
          name="fractmaxwell_dashpot",
          # Model parameters,
          p = [:η, :cᵦ, :β],
          # Relaxation modulus
          G = quote
                cᵦ*t^(-β)*mittleff(1 - β, 1 - β, -cᵦ*t^(1 - β)/η)
              end,
          # Creep modulus
          J = quote
                t/η + t^β/(cᵦ*gamma(1 + β))
              end,
          # Storage modulus
          Gp = quote
                  denominator = (η*ω)^2 + (cᵦ*ω^β)^2 + 2*(η*ω)*(cᵦ*ω^β)*cos((1-β)*π/2)
                  numerator = ((η*ω)^2)*(cᵦ*ω^β)*cos(β*π/2)
                  numerator/denominator
               end,
         # Loss modulus
          Gpp = quote
                  denominator = (η*ω)^2 + (cᵦ*ω^β)^2 + 2*(η*ω)*(cᵦ*ω^β)*cos((1-β)*π/2)
                  numerator = ((cᵦ*ω^β)^2)*(η*ω) + ((η*ω)^2)*(cᵦ*ω^β)*sin(β*π/2)
                  numerator/denominator
                end,
          # Constraints
          constraint = quote
                   (β<1) & (β>0)
                  end,
          # Network
          info= "
                  ___
              _____| |_________╱╲____
                  _|_|         ╲╱
                    η            cᵦ, β
                 "
          )


Maxwell = RheoModelClass(
        # Model name
        name="maxwell",
        # Model parameters,
        p = [:η, :k],
        # Relaxation modulus
        G = quote
              k*exp(-k*t/η)
            end,
        # Creep modulus
        J = quote
              t/η + 1/k
            end,
        # Storage modulus
        Gp = quote
                denominator = 1 + η^2*ω^2/k^2
                numerator = η^2*ω^2/k
                numerator/denominator
             end,
       # Loss modulus
        Gpp = quote
                denominator = 1 + η^2*ω^2/k^2
                numerator = η*ω
                numerator/denominator
              end,
        # Network
        info= "
                ___
            _____| |________╱╲  ╱╲  ╱╲  ___
                _|_|          ╲╱  ╲╱  ╲╱
                  η                  k
               "
        )
