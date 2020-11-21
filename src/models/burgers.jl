
BurgersLiquid =  RheoModelClass(
        # Model name
        name="burgerliquid",
        # Model parameters,
        p = [:η₁ , :k₁, :η₂, :k₂],
        # Relaxation modulus
        G = quote
              p1 = η₁/k₁ + η₁/k₂ + η₂/k₂
              p2 = η₁ * η₂ / (k₁ * k₂)
              q1 = η₁
              q2 = η₁ * η₂ / k₂
              A = sqrt(p1^2 - 4*p2)
              r1 = (p1-A)/(2*p2)
              r2 = (p1+A)/(2*p2)
              ((q1-q2*r1)*exp(-r1*t) - (q1-q2*r2)* exp(-r2*t))/A
            end,
        # Creep modulus
        J = quote
              1/k₁ + t/ η₁ + 1/k₂ * (1- exp(-k₂*t/η₂))
            end,
        # Storage modulus
        Gp = quote
                p1 = η₁/k₁ + η₁/k₂ + η₂/k₂
                p2 = η₁ * η₂ / (k₁ * k₂)
                q1 = η₁
                q2 = η₁ * η₂ / k₂
                numerator = p1 * q1 * ω^2 - q2 * ω^2 * (1-p2*ω^2)
                denominator = p1^2 * ω^2 + (1-p2*ω^2)^2
                numerator/denominator
             end,
        # Loss modulus
        Gpp = quote
                p1 = η₁/k₁ + η₁/k₂ + η₂/k₂
                p2 = η₁ * η₂ / (k₁ * k₂)
                q1 = η₁
                q2 = η₁ * η₂ / k₂
                numerator = p1 * q2 * ω^3 + q1 * ω * (1-p2*ω^2)
                denominator = p1^2 * ω^2 + (1-p2*ω^2)^2
                numerator/denominator
              end,
        # Network description
        description = (type = "series", components=((:Dashpot, (:η₁,)), (:Spring, (:k₁,)), (:KelvinVoigt, (:η₂, :k₂)))),
        # Network
        info= "

                                                         ___
                                                 _________| |________
                 ___                            |        _|_| η₂     |
             _____| |________╱╲  ╱╲  ╱╲  _______|                    |____
                 _|_|          ╲╱  ╲╱  ╲╱       |                    |
                    η₁                  k₁      |____╱╲  ╱╲  ╱╲  ____|
                                                       ╲╱  ╲╱  ╲╱
                                                                 k₂

                "
        )
