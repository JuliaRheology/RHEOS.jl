#!/usr/bin/env julia

Fract_Zener = RheoModelClass(
          # Model name
          name="frac_zener",
          # Model parameters,
          p = [ :cₐ, :a, :cᵦ, :β, :cᵧ, :γ],
          # Relaxation modulus
          G = quote
                 cᵦ*t^(-β)*mittleff(a - β, 1 - β, -cᵦ*t^(a - β)/cₐ) + cᵧ*t^(-γ)/gamma(1 - γ)
              end,
          # Creep modulus
          J = quote
                  Ĵ(s) = (1/s)*(cₐ*s^a + cᵦ*s^β)/(cₐ*s^a*cᵦ*s^β + cᵧ*s^γ*(cₐ*s^a + cᵦ*s^β))
                  InverseLaplace.talbot(s -> Ĵ(s), t)
              end,
          # Storage modulus
          Gp = quote
                   denominator = (cₐ*ω^a)^2 + (cᵦ*ω^β)^2 + 2*(cₐ*ω^a)*(cᵦ*ω^β)*cos((a - β)*π/2)
                   numerator = ((cᵦ*ω^β)^2)*(cₐ*ω^a)*cos(a*π/2) + ((cₐ*ω^a)^2)*(cᵦ*ω^β)*cos(β*π/2)
                   numerator/denominator + cᵧ*ω^γ*cos(γ*π/2)
               end,
         # Loss modulus
          Gpp = quote
                   denominator = (cₐ*ω^a)^2 + (cᵦ*ω^β)^2 + 2*(cₐ*ω^a)*(cᵦ*ω^β)*cos((a - β)*π/2)
                   numerator = ((cᵦ*ω^β)^2)*(cₐ*ω^a)*sin(a*π/2) + ((cₐ*ω^a)^2)*(cᵦ*ω^β)*sin(β*π/2)
                   numerator/denominator + cᵧ*ω^γ*sin(γ*π/2)
                end,
          # Constraints
          constraint = quote
                   all([   (a<1) & (a>0)
                           (β<1) & (β>0)
                            -a+β < 0
                           (γ<1) & (γ>0)] )
                  end,
          # Network
          info= "

                  ______╱╲__________╱╲______
                 |      ╲╱          ╲╱      |
          _______|      cₐ,a         cᵦ, β  |_______
                 |                          |
                 |____________╱╲____________|
                              ╲╱
                              cᵧ, γ
                     "
          )


FractSLS_Zener = RheoModelClass(
        # Model name
        name="fracsls_Zener",
        # Model parameters,
        p = [ :cₐ, :a, :kᵦ, :kᵧ],
        # Relaxation modulus
        G = quote
                kᵦ*mittleff(a, -(kᵦ/cₐ)*t^a) + kᵧ
            end,
        # Creep modulus
        J = quote
                Ĵ(s) = (1/s)*(cₐ*s^a + kᵦ)/(cₐ*s^a*kᵦ + kᵧ*(cₐ*s^a + kᵦ))
                 InverseLaplace.talbot(s -> Ĵ(s), t)
            end,
        # Storage modulus
        Gp = quote
                denominator = (cₐ*ω^a)^2 + kᵦ^2 + 2*(cₐ*ω^a)*kᵦ*cos(a*π/2)
                numerator = kᵦ^2*(cₐ*ω^a)*cos(a*π/2) + ((cₐ*ω^a)^2)*kᵦ
                numerator/denominator + kᵧ
             end,
       # Loss modulus
        Gpp = quote
                denominator = (cₐ*ω^a)^2 + kᵦ^2 + 2*(cₐ*ω^a)*kᵦ*cos(a*π/2)
                numerator = kᵦ^2*(cₐ*ω^a)*sin(a*π/2)
                numerator/denominator
              end,
        # Constraints
        constraint = quote
                 (a<1) & (a>0)
                end,
        # Network
        info= "

                _____╱╲_______╱╲  ╱╲  ╱╲  _____
               |     ╲╱         ╲╱  ╲╱  ╲╱     |
        _______|       cₐ,a              kᵦ    |_______
               |                               |
               |__________╱╲  ╱╲  ╱╲  _________|
                            ╲╱  ╲╱  ╲╱
                                 kᵧ
                   "
        )


SLS_Zener = RheoModelClass(
          # Model name
          name="SLS_Zener",
          # Model parameters,
          p = [ :η, :kᵦ, :kᵧ],
          # Relaxation modulus
          G = quote
                 kᵧ + kᵦ*exp(-t*kᵦ/η)
              end,
          # Creep modulus
          J = quote
                  c₀ = 1/kᵧ
                  c₁ = kᵦ/(kᵧ*(kᵧ + kᵦ))
                  τᵣ = η*(kᵧ + kᵦ)/(kᵧ*kᵦ)
                  c₀ - c₁*exp(-t/τᵣ)
              end,
          # Storage modulus
          Gp = quote
                  τ = η/kᵦ
                  denominator = 1 + τ^2*ω^2
                  numerator = ω^2*τ^2*kᵦ
                  numerator/denominator + kᵧ
               end,
         # Loss modulus
          Gpp = quote
                  τ = η/kᵦ
                  denominator = 1 + τ^2*ω^2
                  numerator = ω*τ*kᵦ
                  numerator/denominator
                end,
          # Network
          info= "
                      ___
                  _____| |________╱╲  ╱╲  ╱╲  ___
                 |    _|_|          ╲╱  ╲╱  ╲╱   |
          _______|      η                  kᵦ    |_______
                 |                               |
                 |__________╱╲  ╱╲  ╱╲  _________|
                              ╲╱  ╲╱  ╲╱
                                   kᵧ
                     "
          )

FractJeffreys_Zener = RheoModelClass(
                    # Model name
                    name="fjeff_Zener",
                    # Model parameters,
                    p = [ :ηₐ, :cᵦ, :β, :ηᵧ],
                    # Relaxation modulus
                    G = quote
                            diracterm = t!=0.0 ? 0.0 : Inf
                            cᵦ*t^(-β)*mittleff(1-β, 1-β, -cᵦ*t^(1-β)/ηₐ) + ηᵧ*diracterm
                        end,
                    # Creep modulus
                    J = quote
                            Ĵ(s) = (1/s)*(ηₐ*s + cᵦ*s^β)/(ηₐ*s*cᵦ*s^β + ηᵧ*s*(ηₐ*s + cᵦ*s^β))
                            InverseLaplace.talbot(s -> Ĵ(s), t)
                        end,
                    # Storage modulus
                    Gp = quote
                            denominator = (ηₐ*ω)^2 + (cᵦ*ω^β)^2 + 2*(ηₐ*ω)*(cᵦ*ω^β)*cos((1 - β)*π/2)
                            numerator = ((ηₐ*ω)^2)*(cᵦ*ω^β)*cos(β*π/2)
                            # cosine floating point error work-around
                            if β==1.0
                                0.0
                            else
                                numerator/denominator
                            end
                         end,
                   # Loss modulus
                    Gpp = quote
                            denominator = (ηₐ*ω)^2 + (cᵦ*ω^β)^2 + 2*(ηₐ*ω)*(cᵦ*ω^β)*cos((1 - β)*π/2)
                            numerator = ((cᵦ*ω^β)^2)*(ηₐ*ω) + ((ηₐ*ω)^2)*(cᵦ*ω^β)*sin(β*π/2)
                            numerator/denominator + ηᵧ*ω
                          end,
                    # Constraints
                    constraint = quote
                             (β<1) & (β>0)
                            end,
                    # Network
                    info= "

                                    ___
                            _________| |_________╱╲________
                           |        _|_|         ╲╱        |
                    _______|          ηₐ            cᵦ, β  |_______
                           |              ___              |
                           |_______________| |_____________|
                                          _|_|
                                             ηᵧ
                               "
                    )


Jeffreys_Zener = RheoModelClass(
                # Model name
                name="jeffreys_Zener",
                # Model parameters,
                p = [:ηₐ, :k, :ηᵧ],
                # Relaxation modulus
                G = quote
                        diracterm = t!=0.0 ? 0.0 : Inf
                        k*exp(-k*t/ηₐ) + ηᵧ*diracterm
                    end,
                # Creep modulus
                J = quote
                        Ĵ(s) = (1/s)*(ηₐ*s + k)/(ηₐ*s*k + ηᵧ*s*(ηₐ*s + k))
                        InverseLaplace.talbot(s -> Ĵ(s), t)
                    end,
                # Storage modulus
                Gp = quote
                        denominator = (ηₐ*ω)^2 + k^2
                        numerator = ((ηₐ*ω)^2)*k
                        numerator/denominator
                     end,
               # Loss modulus
                Gpp = quote
                        denominator = (ηₐ*ω)^2 + k^2
                        numerator = (k^2)*(ηₐ*ω)
                        numerator/denominator + ηᵧ*ω
                      end,
                # Network
                info= "

                              ___
                        _______| |_______╱╲  ╱╲  ╱╲  ___
                       |      _|_|         ╲╱  ╲╱  ╲╱  |
                _______|          ηₐ            k      |_______
                       |              ___              |
                       |_______________| |_____________|
                                      _|_|
                                         ηᵧ
                           "
                )

FractSolid = RheoModelClass(
        # Model name
        name="fractsolid",
        # Model parameters,
        p = [:η, :cᵦ, :β, :k],
        # Relaxation modulus
        G = quote
              k + cᵦ*t^(-β)*mittleff(1 - β, 1 - β, -cᵦ*(t^(1 - β))/η)
            end,
        # Creep modulus
        J = quote
              a = η/cᵦ
              b = k*η/cᵦ
              Ĵ(s) = (1/s^2)*(1+a*s^(1-β))/(η+k/s+b*s^(-β))
              InverseLaplace.talbot(s -> Ĵ(s), t)
            end,
        # Storage modulus
        Gp = quote
               denominator = (η*ω)^2 + (cᵦ*ω^β)^2
               numerator = ((η*ω)^2)*(cᵦ*ω^β)*cos(β*π/2)
               numerator/denominator + k
             end,
       # Loss modulus
        Gpp = quote
                denominator = (η*ω)^2 + (cᵦ*ω^β)^2
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
                  _____| |__________╱╲__________
                 |    _|_|          ╲╱          |
             ___ |      η              cᵦ, β    |___
                 |                              |
                 |__________╱╲  ╱╲  ╱╲  ________|
                              ╲╱  ╲╱  ╲╱
                                k
               "
        )
