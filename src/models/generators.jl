
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


(p,e)=G_GSLS(3)


# Cool replacement function inspired from
# https://stackoverflow.com/questions/29778698/julia-lang-metaprogramming-turn-expression-into-function-with-expression-depend
# This probably could go in its own module as this is very useful.
#

expr_replace!(ex, s, v) = ex
expr_replace!(ex::Symbol, s, v) = s == ex ? v : ex

function expr_replace!(ex::Expr, s, v)
    for i=1:length(ex.args)
        ex.args[i] = expr_replace!(ex.args[i], s, v)
    end
    return ex
end

expr_replace(ex, s, v) = expr_replace!(copy(ex), s, v)

function expr_replace(ex, nt)
    e = copy(ex)
    k=keys(nt)
    v=values(nt)
    for i in 1:length(nt)
        expr_replace!(e, k[i], v[i])
    end
    return e
end

expr_replace(e,(k_0=0.5,k_1=1.123,k_2=3.14))


# cd("/home/alexandre/.julia/dev/RHEOS/src/"); include("RHEOS.jl"); using Main.RHEOS; m=RheoModel(SLS2,(G₀=1., G₁=2., η₁=3., G₂=4., η₂=5.))
