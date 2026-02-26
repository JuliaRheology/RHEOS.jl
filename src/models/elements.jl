#!/usr/bin/env julia

Springpot =  RheoModelClass(
        # Model name
        name="springpot",
        # Model parameters,
        p = (:cᵦ, :β),
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

        #Differential equation
        equation = (ϵ =((:cᵦ,:β),), σ =((1.0,0.0),)),
        
        # Constraints
        constraint = [quote
                        β-1
                end,
                quote
                        -β
                end
                ],

        # constraint = quote
        #         (0<β<1)
        #         end,

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
        p = (:k,),
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
        
              #TODO: Placeholder eq
        equation = (ϵ =((1.0,1.0),), σ =((1.0,1.0),)),

        
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
        p = (:η,),
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
        equation = (ϵ =((:η,1.0),), σ =((1.0,0.0),)),

        # Network
        info= "
                 ___
             _____| |_____
                 _|_|
                     η
                "
        )


