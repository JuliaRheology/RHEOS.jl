using Compat
import Optim

@compat abstract type AbstractWeeks <: AbstractILt end

function Base.show{T<:AbstractWeeks}(io::IO, w::T)
    print(io, string(typeof(w)), "(Nterms=", w.Nterms, ",sigma=", w.sigma, ",b=", w.b,')')
end

#### Weeks

type Weeks <: AbstractWeeks
    func::Function
    Nterms::Int
    sigma::Float64
    b::Float64
    coefficients::Array{Float64,1}
end

function _get_coefficients(func, Nterms, sigma, b)
    a0 = real(_wcoeff(func,Nterms,sigma,b))
    a = a0[Nterms+1:2*Nterms]
end

_get_coefficients(w::Weeks) = _get_coefficients(w.func, w.Nterms, w.sigma, w.b)

_set_coefficients(w::Weeks) = (w.coefficients = _get_coefficients(w))

"""
   w::Weeks = Weeks(func::Function, Nterms::Integer=64, sigma=1.0, b=1.0)

return `w`, which estimates the inverse Laplace transform of `func` with
the Weeks algorithm. `w(t)` evaluates the transform at `t`. The accuracy depends on the choice
of `sigma` and `b`, with the optimal choices depending on `t`.

The call to `Weeks` that creates `w` is expensive relative to evaluation via `w(t)`.

# Example

Compute the inverse transform of the transform of `cos` at argument `pi/2`.
```
julia> ft = Weeks(s -> s/(s^2+1), 80);

julia> ft(pi/2)
0.0
```
"""
Weeks(func::Function, Nterms::Integer=64, sigma=1.0, b=1.0) =  Weeks(func,Nterms,sigma,b,_get_coefficients(func,Nterms,sigma,b))

function eval_ilt(w::Weeks, t)
    L = _laguerre(w.coefficients,2*w.b*t)
    L * exp((w.sigma-w.b)*t)
end

function eval_ilt(w::Weeks, t::AbstractArray)
    L = _laguerre(w.coefficients,2*w.b*t)
    [ L1 * exp((w.sigma-w.b)*t1) for (t1,L1) in zip(t,L)]
end

function optimize(w::Weeks, t)
    (w.sigma, w.b) = _optimize_sigma_and_b(w.func, t, w.Nterms, 0.0, 30, 30)
    w.coefficients = _get_coefficients(w)
    w
end

"""
    optimize(w::AbstractWeeks, t, Nterms)

optimize the parameters of the inverse Laplace transform `w` at the
argument `t`. If `Nterms` is ommitted, the current value of `w.Nterms`
is retained.

The accuracy of the Weeks algorithm depends strongly on `t`. For some ranges
of `t`, the accuracy is relatively insensitive to the parameters. For other
values of `t`, even using optimized parameters results in estimates of
the inverse transform that are extremely inaccurate.

`optimize` is expensive in CPU time and allocation, it performs nested single-parameter
optimization over two parameterss.
"""
function optimize(w::AbstractWeeks, t, N)
    w.Nterms = N
    optimize(w,t)
end

"""
    opteval{T<:AbstractWeeks}(w::T, t, Nterms)

estimate an inverse Laplace transform at argument `t` using `w` after
optimizing the parameters for `t`. If `Nterms` is omitted, then the
current value of `w.Nterms` is used.

Use `Weeks` or `WeeksErr` to create `w`.
"""
function opteval(w::AbstractWeeks, t, N)
    optimize(w,t, N)
    w(t)
end

opteval(w::AbstractWeeks, t) = opteval(w,t,w.Nterms)

"""
    setparameters{T<:AbstractWeeks}(w::T, sigma, b, Nterms)

Set the parameters for the inverse Laplace transform object `w` and recompute
the internal data. Subsequent calls `w(t)` will use these parameters. If `Nterms`
or both `Nterms` and `b` are omitted, then their current values are retained.
"""
function setparameters(w::AbstractWeeks, sigma, b, Nterms)
    (w.sigma, w.b, w.Nterms) = (sigma,b,Nterms)
    _set_coefficients(w)
    w
end

setparameters(w,sigma,b) = setparameters(w,sigma,b,w.Nterms)
setparameters(w,sigma) = setparameters(w,sigma,w.b)

function setNterms(w::AbstractWeeks, Nterms::Integer)
    w.Nterms = Nterms
    _set_coefficients(w)
    w
end

#### WeeksErr

type WeeksErr <: AbstractWeeks
    func::Function
    Nterms::Int
    sigma::Float64
    b::Float64
    coefficients::Array{Float64,1}
    sa1::Float64
    sa2::Float64
end

function _get_coefficients_and_params(func, Nterms, sigma, b)
    M = 2 * Nterms
    a0 = real(_wcoeff(func,M,sigma,b))
    a1 = a0[2*Nterms+1:3*Nterms]
@compat   sa1 = sum(abs.(a1))
@compat   sa2 = sum(abs.(@view a0[3*Nterms+1:4*Nterms]))
    (a1,sa1,sa2)
end

_get_coefficients_and_params(w::WeeksErr) = _get_coefficients_and_params(w.func, w.Nterms, w.sigma, w.b)

_get_coefficients(w::WeeksErr) = _get_coefficients_and_params(w.func, w.Nterms, w.sigma, w.b)

_set_coefficients(w::WeeksErr) = (w.coefficients, w.sa1, w.sa2) = _get_coefficients_and_params(w)

function optimize(w::WeeksErr, t)
    (w.sigma, w.b) = _optimize_sigma_and_b(w.func, t, w.Nterms, 0.0, 30, 30)
    _set_coefficients(w)
#    (w.coefficients, w.sa1, w.sa2) = _get_coefficients_and_params(w)
    w
end

"""
   w::WeeksErr = WeeksErr(func::Function, Nterms::Integer=64, sigma=1.0, b=1.0)

return `w`, which estimates the inverse Laplace transform of `func` via the Weeks algorithm.
`w(t)` returns a tuple containing the inverse transform at `t` and an error estimate. The accuracy of the
inversion depends on the choice of `sigma` and `b`.

# Example

Compute the inverse transform of the transform of `cos`, and an error estimate, at argument `pi/2` using `80` terms.
```
julia> ft = WeeksErr(s -> s/(s^2+1), 80);

julia> ft(pi/2)
(0.0,3.0872097665938698e-15)
```
This estimate is more accurate than `cos(pi/2)`.
```
julia> ft(pi/2)[1] - cos(pi/2)
-6.123233995736766e-17

julia> ft(pi/2)[1] - 0.0         # exact value
0.0

julia> ft(pi/2)[1] - cospi(1/2)  # cospi is more accurate
0.0
```
"""
function WeeksErr(func::Function, Nterms::Integer=64, sigma=1.0, b=1.0)
    params = _get_coefficients_and_params(func, Nterms, sigma, b)
    WeeksErr(func,Nterms,sigma,b,params...)
end

function eval_ilt(w::WeeksErr, t)
    L = _laguerre(w.coefficients,2*w.b*t)
    f = L * exp((w.sigma-w.b)*t)
    est = exp(w.sigma*t)*(w.sa2+eps()*w.sa1)
    (f,est)
end

function eval_ilt(w::WeeksErr, t::AbstractArray)
    L = _laguerre(w.coefficients,2*w.b*t)
@compat    f = L .* exp.((w.sigma-w.b)*t)
@compat    est = exp.(w.sigma*t)*(w.sa2+eps()*w.sa1)
    (f,est)
end

for w in (:Weeks, :WeeksErr)
    @eval (w::$(w))(t) = eval_ilt(w,t)
end

##### internal functions

function _wcoeff(F,N,sig,b)
    n = -N:1:N-1
    h = pi/N
    th = h*(n+1//2)
@compat    y = b*cot.(th/2)
    imunit = Complex(zero(eltype(y)),one(eltype(y)))
    s = sig +  imunit * y
    FF0 = map(F, s)
    FF = [ FF1 * (b + imunit * y1) for (FF1,y1) in zip(FF0,y)]
    a = fftshift(fft(fftshift(FF)))/(2*N)
@compat   exp.(Complex(zero(h),-one(h)) * n*h/2) .* a
end

function _laguerre(a::AbstractVector,x::AbstractArray)
    N = length(a) - 1
    unp1 = zeros(x)
    un = a[N+1]*ones(x)
    local unm1
    for n in N:-1:1
#        unm1 =  (1//n)*(2*n-1-x) .* un - n/(n+1)*unp1 + a[n]
        unm1 =  [(1//n)*(2*n-1-x0) * un0 - n/(n+1)*unp10 + a[n] for (x0,un0,unp10) in zip(x,un,unp1)]
        unp1 = un
        un = unm1
    end
    unm1
end

function _laguerre(a::AbstractVector,x)
    N = length(a) - 1
    unp1 = zero(x)
    un = a[N+1]*one(x)
    local unm1
    for n in N:-1:1
        unm1 = (1//n)*(2*n-1-x) * un - n/(n+1)*unp1 + a[n]
        unp1 = un
        un = unm1
    end
    unm1
end

function _optimize_sigma_and_b(F, t, N, sig0, sigmax, bmax)
    sigma_opt = Optim.minimizer(Optim.optimize( sig -> _werr2e(sig, F,t,N, sig0, sigmax, bmax), sig0, sigmax))
    b_opt = Optim.minimizer(Optim.optimize( b -> _werr2t(b, F, N, sigma_opt), 0, bmax))
    (sigma_opt, b_opt)
end

function _werr2e(sig,F,t,N,sig0,sigmax,bmax)
    b = Optim.minimizer(Optim.optimize( (b) -> _werr2t(b, F,N, sig) , 0.0, bmax))
    M = 2*N
    a = _wcoeff(F,M,sig,b)
    a1 = @view a[2*N+1:3*N]
@compat   sa1 = sum(abs.(a1))
    a2 = @view a[3*N+1:4*N]
@compat   sa2 = sum(abs.(a2))
    sig*t + log(sa2+eps()*sa1)
end

function _werr2t(b, F, N, sig)
    M = 2*N
    a = _wcoeff(F,M,sig,b)
@compat  sa2 = sum(abs.( @view a[3*N+1:4*N]))
    log(sa2)
end
