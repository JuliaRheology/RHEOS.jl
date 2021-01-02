﻿<a name="logo"/>
<div align="center">
<img src="docs/Logo.png" height="130"></img>
</a>
</div>

|**Test Status**|**Documentation**|**Test Coverage**|**Paper**|
|:-------------:|:---------------:|:---------------:|:-------:|
| [![Build Status][travis-img]][travis-url] | [![Stable][docs-sta-img]][docs-sta-url] [![Dev][docs-dev-img]][docs-dev-url] | [![codecov][codecov-img]][codecov-url] [![coveralls][coveralls-img]][coveralls-url] | [![status][joss-img]][joss-url] |

# RHEOS - RHEology, Open-Source
RHEOS, an abbreviation of Rheology Open Source, is a software package written in the Julia programming language that provides tools for analyzing rheological data. Features include:

- Stress/Strain/Time data can be easily be fitted to a viscoelastic model

- G'/G''/Frequency data can easily be fitted to a viscoelastic model

- Many standard and fractional viscoelastic models have already been implemented within RHEOS new ones can easily be added by users

- A fitted model can be used to predict the behaviour of the material under other loading conditions, enabling the fit/predict paradigm of model selection

- Artificial loading conditions can be generated within RHEOS to better understand a model's response

## Installation
1. Install the latest version of Julia
2. From Julia interactive command-line REPL, enter pkg mode by pressing ```]```
3. (Optional) Enable desired Project.toml environment
4. Run the command ```add RHEOS```

## Documentation
If you installed RHEOS using the instructions above then you will have the _latest stable release_ of RHEOS; to access the documentation for this version [click here][docs-sta-url] or the blue `docs/stable` badge at the top of this README page. To access the _latest documentation built directly from the master branch_ [click here][docs-dev-url] or the blue `docs/dev` badge at the top of this README page. 

## Citation
If you use RHEOS in your work, please consider citing the following papers as appropriate:

+ J. L. Kaplan, A. Bonfanti, A. J. Kabla (2019). _RHEOS.jl -- A Julia Package for Rheology Data Analysis_. Journal of Open Source Software, 4(41), 1700, [https://doi.org/10.21105/joss.01700](https://doi.org/10.21105/joss.01700)

+ A. Bonfanti, J. L. Kaplan, G. Charras, A. J. Kabla (2020). _Fractional viscoelastic models for power-law materials_. Soft Matter, 16, 6002-6020, [https://doi.org/10.1039/D0SM00354A](https://doi.org/10.1039/D0SM00354A)

## Embedded Dependencies
#### [MittagLeffler.jl](https://github.com/jlapeyre/MittagLeffler.jl)

## Contributing to RHEOS
If you believe you have found any bugs or invalid behaviour in RHEOS, please feel free to file an issue on this repository. You can also raise an issue if you feel that any part of the documentation needs clarification, or for any feature requests. Even better than just raising an issue, you could both raise an issue and issue a pull request which fixes that issue. Note that meta-documentation on running tests and building documentation locally is available at the [JuliaRheology/RheoHelpDocs](https://github.com/JuliaRheology/RheoHelpDocs) repository. Please be aware that RHEOS is released with a [Contributor Code of Conduct](CONDUCT.md) and by participating in this project you agree to abide by its terms.


## References
+ W. N. Findley, J. S. Lai, K. Onaran (1989). *Creep and Relaxation of Nonlinear Viscoelastic Materials (with an Introduction to Linear Viscoelasticity)*, Dover Publications, New York. 

+ S. G. Johnson. *The NLopt nonlinear-optimization package*, https://github.com/stevengj/nlopt

+ J. Bezanson, A. Edelman, S. Karpinski, V. B. Shah (2017). *Julia: A Fresh Approach to Numerical Computing*, SIAM Review, doi: 10.1137/141000671.

[travis-img]: https://travis-ci.org/JuliaRheology/RHEOS.jl.svg?branch=master
[travis-url]: https://travis-ci.org/JuliaRheology/RHEOS.jl

[docs-sta-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-sta-url]: https://JuliaRheology.github.io/RHEOS.jl/stable

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://JuliaRheology.github.io/RHEOS.jl/dev

[codecov-img]: https://codecov.io/gh/JuliaRheology/RHEOS.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/JuliaRheology/RHEOS.jl

[coveralls-img]: https://coveralls.io/repos/github/JuliaRheology/RHEOS.jl/badge.svg?branch=master
[coveralls-url]: https://coveralls.io/github/JuliaRheology/RHEOS.jl?branch=master

[license-img]: https://img.shields.io/badge/License-MIT-ff69b2.svg?style=flat

[joss-img]: https://joss.theoj.org/papers/553250d815e1990db1b89c742854c71a/status.svg
[joss-url]: https://joss.theoj.org/papers/553250d815e1990db1b89c742854c71a
