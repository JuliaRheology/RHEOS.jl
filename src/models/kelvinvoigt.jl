#!/usr/bin/env julia

Fract_KelvinVoigt =  RheoModelClass(
        # Model name
        name="fractKV",
        # Model parameters,
        p = [:cₐ, :a, :cᵦ, :β],
        # Relaxation modulus
        G = quote
              cₐ*t^(-a)/gamma(1 - a) + cᵦ*t^(-β)/gamma(1 - β)
            end,
        # Creep modulus
        J = quote
              (t^(a)/cₐ)*mittleff(a - β, 1 + a, -cᵦ*t^(a - β)/cₐ)
            end,
        # Storage modulus
        Gp = quote
                # cosine floating point error work-around
                if a!=1.0 && β!=1.0
                    cₐ*ω^a*cos(a*π/2) + cᵦ*ω^β*cos(β*π/2)
                elseif a==1.0 && β!=1.0
                    cᵦ*ω^β*cos(β*π/2)
                elseif a==1.0 && β==1.0
                    0.0
                end
             end,
        # Loss modulus
        Gpp = quote
                cₐ*ω^a*sin(a*π/2) + cᵦ*ω^β*sin(β*π/2)
              end,
        # Constraints
        constraint = quote
                 all([  (a<1) & (a>0)
                         (β<1) & (β>0)
                          -a+β < 0])
                end,
        # Network
        info= "
                ________ ╱╲ ________
               |         ╲╱  cₐ, a  |
           ____|                    |____
               |                    |
               |________ ╱╲ ________|
                         ╲╱  cᵦ, β
                "
        )


FractS_KelvinVoigt =  RheoModelClass(
        # Model name
        name="fractSpringKV",
        # Model parameters,
        p = [:cₐ, :a, :k],
        # Relaxation modulus
        G = quote
              cₐ*t^(-a)/gamma(1 - a) + k
            end,
        # Creep modulus
        J = quote
              (t^(a)/cₐ)*mittleff(a, 1 + a, -k*t^a/cₐ)
            end,
        # Storage modulus
        Gp = quote
                if a!=1.0
                    cₐ*ω^a*cos(a*π/2) + k
                else
                    k
                end
             end,
        # Loss modulus
        Gpp = quote
                cₐ*ω^a*sin(a*π/2)
              end,
        # Constraints
        constraint = quote
                 (a<1) & (a>0)
                end,
        # Network
        info= "
                ________ ╱╲ ________
               |         ╲╱  cₐ, a  |
           ____|                    |____
               |                    |
               |____╱╲  ╱╲  ╱╲  ____|
                      ╲╱  ╲╱  ╲╱
                                k
                "
        )

FractD_KelvinVoigt =  RheoModelClass(
        # Model name
        name="fractDashpotKV",
        # Model parameters,
        p = [:η, :cᵦ, :β],
        # Relaxation modulus
        G = quote
              t!=0.0 ? cᵦ*t^(-β)/gamma(1 - β) : Inf
            end,
        # Creep modulus
        J = quote
              (t/η)*mittleff(1 - β, 1 + 1, -cᵦ*t^(1.0 - β)/η)
            end,
        # Storage modulus
        Gp = quote
                # cosine floating point error work-around
                if β==1.0
                    0.0
                else
                    cᵦ*ω^β*cos(β*π/2)
                end
             end,
        # Loss modulus
        Gpp = quote
                η*ω + cᵦ*ω^β*sin(β*π/2)
              end,
        # Constraints
        constraint = quote
                 (β<1) & (β>0)
                end,
        # Network
        info= "
                        ___
                _________| |________
               |        _|_| η      |
           ____|                    |____
               |                    |
               |________ ╱╲ ________|
                         ╲╱  cᵦ, β
                "
                )

KelvinVoigt =  RheoModelClass(
        # Model name
        name="KV",
        # Model parameters,
        p = [:η, :k],
        # Relaxation modulus
        G = quote
              t!=0.0 ? k : Inf
            end,
        # Creep modulus
        J = quote
              (1 - exp(-k*t/η))/k
            end,
        # Storage modulus
        Gp = quote
                k
             end,
        # Loss modulus
        Gpp = quote
                η*ω
              end,
        # Network
        info= "
                        ___
                _________| |________
               |        _|_| η      |
           ____|                    |____
               |                    |
               |____╱╲  ╱╲  ╱╲  ____|
                      ╲╱  ╲╱  ╲╱
                                k
                "
                )
