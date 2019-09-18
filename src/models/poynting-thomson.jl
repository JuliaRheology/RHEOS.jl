#!/usr/bin/env julia

Fract_PT = RheoModelClass(
          # Model name
          name="fPT",
          # Model parameters,
          p = [:cₐ, :a, :cᵦ, :β, :cᵧ, :γ],
          # Relaxation modulus
          G = quote
                  G(s) = (1/s)*(cᵧ*s^γ *(cₐ * s^a + cᵦ * s^β))/(cᵧ*s^γ* + cₐ*s^a + cᵦ*s^β)
                  InverseLaplace.talbot(s -> G(s), t)
              end,
          # Creep modulus
          J = quote
                  t^a/cₐ* mittleff(a - β, 1 + a, -cᵦ*t^(a-β)/cₐ) + t^γ / (cᵧ * gamma(1 + γ))
              end,
          # Storage modulus
          Gp = quote
                  denominator = (cₐ * ω^a)^2 + (cᵦ * ω^β)^2 + (cᵧ * ω^γ)^2 + 2* (cₐ * ω^a) * (cᵦ * ω^β) * cos((a-β)*π/2) + 2* (cₐ * ω^a) * (cᵧ * ω^γ) * cos((a-γ)*π/2) + 2* (cᵦ * ω^β) * (cᵧ * ω^γ)* cos((β-γ)*π/2)
                  numerator = (cᵧ * ω^γ) * cos(γ*π/2) * ((cₐ * ω^a)^2 + (cᵦ * ω^β)^2) + (cᵧ * ω^γ)^2 * ((cₐ * ω^a)*cos(a*π/2)+(cᵦ * ω^β)*cos(β*π/2)) + (cₐ * ω^a) * (cᵦ * ω^β) * (cᵧ * ω^γ) * (cos((a-β-γ)*π/2)+cos((β-a-γ)*π/2))
                  numerator/denominator
               end,
         # Loss modulus
          Gpp = quote
                  denominator = (cₐ * ω^a)^2 + (cᵦ * ω^β)^2 + (cᵧ * ω^γ)^2 + 2* (cₐ * ω^a) * (cᵦ * ω^β) * cos((a-β)*π/2) + 2* (cₐ * ω^a) * (cᵧ * ω^γ) * cos((a-γ)*π/2) + 2* (cᵦ * ω^β) * (cᵧ * ω^γ)* cos((β-γ)*π/2)
                  numerator = (cᵧ * ω^γ) * sin(γ*π/2) * ((cₐ * ω^a)^2 + (cᵦ * ω^β)^2) + (cᵧ * ω^γ)^2 * ((cₐ * ω^a)*sin(a*π/2)+(cᵦ * ω^β)*sin(β*π/2)) + (cₐ * ω^a) * (cᵦ * ω^β) * (cᵧ * ω^γ) * (sin((a-β-γ)*π/2)+sin((β-a-γ)*π/2))
                  numerator/denominator
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
                     _________╱╲_________
                    |         ╲╱ cₐ, a   |
                ____|                    |______╱╲____
                    |                    |      ╲╱
                    |_________╱╲_________|        cᵧ, γ
                              ╲╱
                                 cᵦ, β
                                 "
          )





FractSLS_PT = RheoModelClass(
          # Model name
          name="fSLS_PT",
          # Model parameters,
          p = [:cₐ, :a, :kᵦ, :kᵧ],
          # Relaxation modulus
          G = quote
                Zkᵦ = (kᵧ)^2/(kᵦ + kᵧ)
                Zcₐ = cₐ * (kᵧ)^2 / (kᵦ + kᵧ)^2
                Zkᵧ = kᵦ * kᵧ / (kᵦ + kᵧ)
                Zkᵦ*mittleff(a, -(Zkᵦ/Zcₐ)*t^a) + Zkᵧ
              end,
          # Creep modulus
          J = quote
                  1/kᵧ + t^(a)/cₐ * mittleff(a, 1 + a, -kᵦ*t^(a)/cₐ)
              end,
          # Storage modulus
          Gp = quote
                  Zkᵦ = (kᵧ)^2/(kᵦ + kᵧ)
                  Zcₐ = cₐ * (kᵧ)^2 / (kᵦ + kᵧ)^2
                  Zkᵧ = kᵦ * kᵧ / (kᵦ + kᵧ)
                  denominator = (Zcₐ*ω^a)^2 + Zkᵦ^2 + 2*(Zcₐ*ω^a)*Zkᵦ*cos(a*π/2)
                  numerator = Zkᵦ^2*(Zcₐ*ω^a)*cos(a*π/2) + ((Zcₐ*ω^a)^2)*Zkᵦ
                  numerator/denominator + Zkᵧ
               end,
         # Loss modulus
          Gpp = quote
                  Zkᵦ = (kᵧ)^2/(kᵦ + kᵧ)
                  Zcₐ = cₐ * (kᵧ)^2 / (kᵦ + kᵧ)^2
                  Zkᵧ = kᵦ * kᵧ / (kᵦ + kᵧ)
                  denominator = (Zcₐ*ω^a)^2 + Zkᵦ^2 + 2*(Zcₐ*ω^a)*Zkᵦ*cos(a*π/2)
                  numerator = Zkᵦ^2*(Zcₐ*ω^a)*sin(a*π/2)
                  numerator/denominator
                end,
        # Constraints
          constraint = quote
                 (a<1) & (a>0)
                end,
          # Network
          info= "
                     _________╱╲_________
                    |         ╲╱ cₐ, a   |
                ____|                    |______╱╲  ╱╲  ╱╲  ____
                    |                    |        ╲╱  ╲╱  ╲╱
                    |____╱╲  ╱╲  ╱╲  ____|                   kᵧ
                           ╲╱  ╲╱  ╲╱
                                     kᵦ
                                 "
          )

SLS_PT = RheoModelClass(
          # Model name
          name="SLS_PT",
          # Model parameters,
          p = [:η, :kᵦ, :kᵧ],
          # Relaxation modulus
          G = quote
                Zkᵦ = (kᵧ)^2/(kᵦ + kᵧ)
                Zη = η * (kᵧ)^2 / (kᵦ + kᵧ)^2
                Zkᵧ = kᵦ * kᵧ / (kᵦ + kᵧ)
                Zkᵧ + Zkᵦ*exp(-t*Zkᵦ/Zη)
              end,
          # Creep modulus
          J = quote
                  1/kᵧ + (1 - exp(-kᵦ*t/η))/kᵦ
              end,
          # Storage modulus
          Gp = quote
                  Zkᵦ = (kᵧ)^2/(kᵦ + kᵧ)
                  Zη = η * (kᵧ)^2 / (kᵦ + kᵧ)^2
                  Zkᵧ = kᵦ * kᵧ / (kᵦ + kᵧ)
                  τ = Zη/Zkᵦ
                  denominator = 1 + τ^2*ω^2
                  numerator = ω^2*τ^2*Zkᵦ
                  numerator/denominator + Zkᵧ
               end,
         # Loss modulus
          Gpp = quote
                  Zkᵦ = (kᵧ)^2/(kᵦ + kᵧ)
                  Zη = η * (kᵧ)^2 / (kᵦ + kᵧ)^2
                  Zkᵧ = kᵦ * kᵧ / (kᵦ + kᵧ)
                  τ = Zη/Zkᵦ
                  denominator = 1 + τ^2*ω^2
                  numerator = ω*τ*Zkᵦ
                  numerator/denominator
                end,
          # Network
          info= "
                             ___
                     _________| |________
                    |        _|_| η      |
                ____|                    |______╱╲  ╱╲  ╱╲  ____
                    |                    |        ╲╱  ╲╱  ╲╱
                    |____╱╲  ╱╲  ╱╲  ____|                   kᵧ
                           ╲╱  ╲╱  ╲╱
                                     kᵦ
                                 "
          )


FractJeffreys_PT = RheoModelClass(
              # Model name
              name="fjeff_PT",
              # Model parameters,
              p = [ :ηₐ, :cᵦ, :β, :ηᵧ],
              # Relaxation modulus
              G = quote
                      Zηₐ = (ηᵧ)^2/(ηₐ + ηᵧ)
                      Zcᵦ = cᵦ * (ηᵧ)^2 / (ηₐ + ηᵧ)^2
                      Zηᵧ = ηₐ * ηᵧ / (ηₐ + ηᵧ)
                      diracterm = t!=0.0 ? 0.0 : Inf
                      Zcᵦ*t^(-β)*mittleff(1-β, 1-β, -Zcᵦ*t^(1-β)/Zηₐ) + Zηᵧ*diracterm
                  end,
              # Creep modulus
              J = quote
                     t/ηₐ * mittleff(1-β, 2, -cᵦ*t^(1-β)/ηₐ) + t/ηᵧ
                  end,
              # Storage modulus
              Gp = quote
                      Zηₐ = (ηᵧ)^2/(ηₐ + ηᵧ)
                      Zcᵦ = cᵦ * (ηᵧ)^2 / (ηₐ + ηᵧ)^2
                      Zηᵧ = ηₐ * ηᵧ / (ηₐ + ηᵧ)
                      denominator = (Zηₐ*ω)^2 + (Zcᵦ*ω^β)^2 + 2*(Zηₐ*ω)*(Zcᵦ*ω^β)*cos((1 - β)*π/2)
                      numerator = ((Zηₐ*ω)^2)*(Zcᵦ*ω^β)*cos(β*π/2)
                      # cosine floating point error work-around
                      if β==1.0
                          0.0
                      else
                          numerator/denominator
                      end
                   end,
             # Loss modulus
              Gpp = quote
                      Zηₐ = (ηᵧ)^2/(ηₐ + ηᵧ)
                      Zcᵦ = cᵦ * (ηᵧ)^2 / (ηₐ + ηᵧ)^2
                      Zηᵧ = ηₐ * ηᵧ / (ηₐ + ηᵧ)
                      denominator = (Zηₐ*ω)^2 + (Zcᵦ*ω^β)^2 + 2*(Zηₐ*ω)*(Zcᵦ*ω^β)*cos((1 - β)*π/2)
                      numerator = ((Zcᵦ*ω^β)^2)*(Zηₐ*ω) + ((Zηₐ*ω)^2)*(Zcᵦ*ω^β)*sin(β*π/2)
                      numerator/denominator + Zηᵧ*ω
                    end,
              # Constraints
              constraint = quote
                       (β<1) & (β>0)
                      end,
              # Network
              info= "
                                      ___
                              _________| |________
                             |        _|_| ηₐ     |        ___
                         ____|                    |_________| |_____
                             |                    |        _|_| ηᵧ
                             |_________╱╲_________|
                                       ╲╱
                                          cᵦ, β
                                          "
              )


Jeffreys_PT = RheoModelClass(
          # Model name
          name="jeffreys_PT",
          # Model parameters,
          p = [:ηₐ, :k, :ηᵧ],
          # Relaxation modulus
          G = quote
              Zηₐ = (ηᵧ)^2/(ηₐ + ηᵧ)
              Zk = k * (ηᵧ)^2 / (ηₐ + ηᵧ)^2
              Zηᵧ = ηₐ * ηᵧ / (ηₐ + ηᵧ)
              diracterm = t!=0.0 ? 0.0 : Inf
              Zk*exp(-Zk*t/Zηₐ) + Zηᵧ*diracterm
          end,
          # Creep modulus
          J = quote
                (- exp(-k*t/ηₐ) + 1 )/k + t/ηᵧ
              end,
          Gp = quote
                  Zηₐ = (ηᵧ)^2/(ηₐ + ηᵧ)
                  Zk = k * (ηᵧ)^2 / (ηₐ + ηᵧ)^2
                  Zηᵧ = ηₐ * ηᵧ / (ηₐ + ηᵧ)
                  denominator = (Zηₐ*ω)^2 + Zk^2
                  numerator = ((Zηₐ*ω)^2)*Zk
                  numerator/denominator
               end,
         # Loss modulus
          Gpp = quote
                  Zηₐ = (ηᵧ)^2/(ηₐ + ηᵧ)
                  Zk = k * (ηᵧ)^2 / (ηₐ + ηᵧ)^2
                  Zηᵧ = ηₐ * ηᵧ / (ηₐ + ηᵧ)
                  denominator = (Zηₐ*ω)^2 + Zk^2
                  numerator = (Zk^2)*(Zηₐ*ω)
                  numerator/denominator + Zηᵧ*ω
                end,
          # Network
          info= "
                             ___
                     _________| |________
                    |        _|_| ηₐ     |        ___
                ____|                    |_________| |_____
                    |                    |        _|_| ηᵧ
                    |____╱╲  ╱╲  ╱╲  ____|
                           ╲╱  ╲╱  ╲╱
                                     k
                                 "
          )
