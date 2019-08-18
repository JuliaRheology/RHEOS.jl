#!/usr/bin/env julia
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
