
SLS2 = RheoModelClass(
        # Model name
        name="SLS2",
        # Model parameters,
        p = [:G₀, :G₁, :η₁, :G₂, :η₂],
        # Relaxation modulus
        G = quote
              G₀ + G₁*exp(-t*G₁/η₁) + G₂*exp(-t*G₂/η₂)
            end,
        # Creep modulus
        J = quote
              tau1 = G₁/η₁
              tau2 = G₂/η₂
              Jbar(s) = (1/s^2)*(G₀/s+G₁*1/(s+tau1)+G₂*1/(s+tau2))
              invLaplace(s -> Jbar(s), t)
            end,
        # Network
        info= "
                                       ___
                 ___╱╲  ╱╲  ╱╲  ________| |_____
                 |    ╲╱  ╲╱  ╲╱       _|_|     |
                 |      G₂              η₂      |
                 |                     ___      |
              ___|__╱╲  ╱╲  ╱╲  ________| |_____|___
                 |    ╲╱  ╲╱  ╲╱       _|_|     |
                 |      G₁              η₁      |
                 |                              |
                 |__________╱╲  ╱╲  ╱╲  ________|
                              ╲╱  ╲╱  ╲╱
                                G₀
               "
        )


#
#  Generate array of parameter symbols and expression of the relaxation modulus
#  for the generalised maxwell model
#

function G_GSLS(n)

    params=[:k0]
    G_terms=[:(k0)]

    for i=1:n
        k=Symbol("k_"*string(i))
        η=Symbol("η_"*string(i))
        e=:($k * exp.( (- $k/$η) * t))

        params=vcat(params, [k,η])
        G_terms=vcat(G_terms,e)
    end
    G_Expr=Expr(:call,:+,G_terms...)

    return (params, G_Expr)
end


# (p,e)=G_GSLS(3)


# cd("/home/alexandre/.julia/dev/RHEOS/src/"); include("RHEOS.jl"); using Main.RHEOS; m=RheoModel(SLS2,(G₀=1., G₁=2., η₁=3., G₂=4., η₂=5.))
