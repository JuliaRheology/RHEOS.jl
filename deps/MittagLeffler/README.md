# MittagLeffler

*Mittag-Leffler function*

[![Build Status](https://travis-ci.org/jlapeyre/MittagLeffler.jl.svg?branch=master)](https://travis-ci.org/jlapeyre/MittagLeffler.jl)

[![Coverage Status](https://coveralls.io/repos/jlapeyre/MittagLeffler.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/jlapeyre/MittagLeffler.jl?branch=master)

[![codecov.io](http://codecov.io/github/jlapeyre/MittagLeffler.jl/coverage.svg?branch=master)](http://codecov.io/github/jlapeyre/MittagLeffler.jl?branch=master)


```julia
mittlefferr(α,β,z,ρ)   # evaluate Mittag-Leffler function with tolerance ρ
mittlefferr(α,z,ρ)     # mittlefferr(α,1,z,ρ)

mittleff(α,β,z)   # evaluate Mittag-Leffler function with tolerance eps()
mittleff(α,z)     # mittleff(α,1,z)
```

Arguments must satisfy `α > 0`, `β` real, `z` real or complex, `ρ>0`.

For `α<1` and/or `abs(z)<1`, accurate, series-only method are used. The series-only methods work
with BigFloat precision for corresponding input types. Some other parameter ranges also use series
or asymptotic methods.

For some arguments, integrals are evaluated with `quadgk`, with no control on errors. Some results
are accurate, others are not.

### Bugs

`mittleff` fails for some arguments. In particular, some of those that evaluate integrals.

### Reference

Rudolfo Gorenflo, Joulia Loutchko and Yuri Loutchko, *Computation of the Mittag-Leffler function and its derivative*,  Fract. Calc. Appl. Anal, **(2002)**
