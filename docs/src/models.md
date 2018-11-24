# Models

## Models included in RHEOS
Severa models (i.e. their relaxation modulus, creep modulus and complex modulus) already implemented in RHEOS. Their constructors are listed below. Note that the `params` argument is always optional, as represented by their containment in square barckets. If left blank then the model's default parameters are used. Additional models are forthcoming.

### Elements
```@docs
Spring
DashPot
SpringPot
```

### Maxwell Type
```@docs
FractionalMaxwell
FractionalMaxwellSpring
FractionalMaxwellDashpot
Maxwell
```

### Kelvin-Voigt Type
```@docs
FractionalKelvinVoigt
FractionalKVspring
FractionalKVdashpot
KelvinVoigt
```

### Zener Type
```@docs
FractionalZener
FractionalSLS
FractionalSpecial
SLS
```

## Creating your own model
If you know some (or all) of the moduli for a model that you would like use but hasn't already been implemented in RHEOS, this section will explain how to quickly import these moduli into a [`RheologyModel`](@ref) object for use with other parts of RHEOS. For the sake of example, we will use the relaxation modulus, storage modulus and loss modulus of the Standard Linear Solid model as defined in RHEOS.
```
function G_sls(t::Vector{T}, params::Vector{T}) where T<:Real
    η, kᵦ, kᵧ = params

    G = kᵧ .+ kᵦ*exp.(-t*kᵦ/η)
end

function Gp_sls(ω::Vector{T}, params::Vector{T}) where T<:Real
    η, kᵦ, kᵧ = params

    τ = η/kᵦ

    denominator = 1 .+ (τ^2)*(ω.^2)
    numerator = (ω.^2)*(τ^2)*kᵦ

    Gp = numerator./denominator .+ kᵧ
end

function Gpp_sls(ω::{T}, params::Vector{T}) where T<:Real
    η, kᵦ, kᵧ = params

    τ = η/kᵦ

    denominator = 1 .+ (τ^2)*(ω.^2)
    numerator = ω*τ*kᵦ

    Gpp = numerator./denominator
end
```
Now we have the our moduli defined as Julia functions we can store them, along with some (optional) default parameters, in a [`RheologyModel`](@ref) struct in the following way.
```
our_model = RheologyModel(G = G_sls, Gp = Gp_sls, Gpp = Gpp_sls, params = [1.0, 0.5, 1.0])
```
Now we can fit this model to data and use it to make predictions. Any moduli not included in this final step will default to a `null_modulus` which always returns the array `[-1.0]`.

