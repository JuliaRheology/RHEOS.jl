#!/usr/bin/env julia

Springpot =  RheoModelClass(
        # Model name
        name="springpot",
        # Model parameters,
        p = [:cᵦ, :β],
        # Relaxation modulus
        G = quote
              cᵦ*t^(-β)/gamma(1 - β)
            end,
        # Creep modulus
        J = quote
              (t^β)/(cᵦ*gamma(1 + β))
            end,
        # Storage modulus
        Gp = quote
                cᵦ*(ω^β)*cos(π*β/2)
             end,
        # Loss modulus
        Gpp = quote
                cᵦ*(ω^β)*sin(π*β/2)
              end,
        # Constraints
        constraint = quote
                (0<β<1)
                end,
        # Network
        info= "
                ____ ╱╲ ____
                     ╲╱  cᵦ, β
                "
        )

Spring =  RheoModelClass(
        # Model name
        name="spring",
        # Model parameters,
        p = [:k],
        # Relaxation modulus
        G = quote
              k
            end,
        # Creep modulus
        J = quote
              1/k
            end,
        # Storage modulus
        Gp = quote
                k
             end,
        # Loss modulus
        Gpp = quote
                0.0
              end,
        # Network description
        description = (type="basic",),
        # Network
        info= "
                ___╱╲  ╱╲  ╱╲  ________
                     ╲╱  ╲╱  ╲╱  k
                "
        )

Dashpot =  RheoModelClass(
        # Model name
        name="dashpot",
        # Model parameters,
        p = [:η],
        # Relaxation modulus
        G = quote
                η
              end,
        # Creep modulus
        J = quote
              t/η
            end,
        # Storage modulus
        Gp = quote
                0.0
             end,
        # Loss modulus
        Gpp = quote
                η*ω
              end,
        # Network description
        description = (type="basic",),
        # Network
        info= "
                 ___
             _____| |_____
                 _|_|
                     η
                ",
        use_G_integral = true
        
        )


