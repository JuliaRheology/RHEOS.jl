# InverseLaplace

*Numerical inverse Laplace transform.*

The source repository is [https://github.com/jlapeyre/InverseLaplace.jl](https://github.com/jlapeyre/InverseLaplace.jl).

## Contents

```@contents
```

## Index

```@index
```

## Inverse Laplace transform methods

Constructing these Julia types, corresponding to different numerical methods,
returns a callable object that evaluates the inverse transform at specified points.

```@docs
ILT
Talbot
GWR
Weeks
WeeksErr
```

## Setting parameters

The inverse Laplace tranform routines should not be treated as black boxes. They are prone to instability and can give inaccurate or
wrong results. There are some parameters you can set to try to minimize these problems.

```@docs
setNterms
optimize
opteval
setparameters
```

## Analzying accuracy

```@docs
ILtPair
abserr
TransformPair
iltpair_power
```

## Lower-level interface

Some of the lower-level routines can be called directly, without constructing types defined in `InverseLaplace`.

```@docs
ilt
talbot
gwr
InverseLaplace.talbotarr
```


## References

J.A.C Weideman, *Algorithms for Parameter Selection in the Weeks Method for Inverting the Laplace Transform,
SIAM Journal on Scientific Computing*, Vol. 21, pp. 111-128 **(1999)**


Abate, J. and Valkó, P.P., *Multi-precision Laplace transform inversion
International Journal for Numerical Methods in Engineering*, Vol. 60 (Iss. 5-7) pp 979–993 **(2004)**

Valkó, P.P. and Abate, J.,
*Comparison of Sequence Accelerators for the Gaver Method of Numerical Laplace Transform Inversion*,
Computers and Mathematics with Application,  Vol. 48 (Iss.3-40) pp. 629-636 **(2004)**
