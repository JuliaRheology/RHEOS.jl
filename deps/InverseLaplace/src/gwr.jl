# Valkó, P.P. and Abate, J.
# Comparison of Sequence Accelerators for the Gaver Method of Numerical Laplace Transform Inversion
# Computers and Mathematics with Application,  Vol. 48 (Iss.3-40) 2004 pp. 629-636
# Gaver Wynn rho method

# type GWR <: AbstractILt
#     func::Function
#     Nterms::Int
# end

# GWR(func::Function) = GWR(func, 16)

"""
    gwr(func::Function, t::AbstractFloat, M::Integer=16)

Evaluate the inverse Laplace transform of `func` at the point `t`. Use `M` terms in the algorithm.
For `typeof(t)` is `Float64`, the default for `M` is `16`. For `BigFloat` the default is `64`.

If `BigFloat` precision is larger than default, try increasing `M`.

# Example

```jldoctest
julia> InverseLaplace.gwr( s -> 1/s^3,  3.0)
4.499985907607361
```

!!! note

    This function uses the Gaver-Wynn rho method.
    It evaluates `func` only for real arguments.

"""
function gwr(func, t, M)
    Dt = typeof(t)
    bM = convert(Dt,M)
    tau = log(convert(Dt,2))/t
    broken = false
    Fi = Array{Dt}(2 * M)
    for i in 1: 2 * M
        Fi[i] = func(i * tau)
    end
    M1 = M
    G0 = zeros(Dt,M1+1)
    for n in 1:M
        sm = zero(Dt)
        bn = convert(Dt,n)
        for i in 0:n
            bi = convert(Dt,i)
            sm += binomial(big(n),big(i)) * (-1)^i * Fi[n+i]
        end
        G0[n] = tau * factorial(2*bn)/(factorial(bn)*factorial(bn-1)) * sm
    end
    Gm = zeros(Dt,M1+1)
    Gp = zeros(Dt,M1+1)
    best = G0[M1]
    for k in 0:M1-2
        for n in (M1-2-k):-1:0
            expr = G0[n+2] - G0[n+1]
            if expr == 0
                broken = true
                break
            end
            expr = Gm[n+2] + (k+1)/expr
            Gp[n+1] = expr
            if isodd(k) && n == M1 - 2 - k
                best = expr
            end
        end
        if broken break end
        for n in 0:(M1-k)
            Gm[n+1] = G0[n+1]
            G0[n+1] = Gp[n+1]
        end
    end
    best
end

gwr(func, t::Float64) = gwr(func,t,16)
gwr(func, t::BigFloat) = gwr(func,t,64)
gwr(func, t::Integer, args...) = gwr(func,BigFloat(t), args...)
gwr(func, t::Rational, args...) = gwr(func,float(t), args...)
