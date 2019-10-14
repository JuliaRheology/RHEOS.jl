<a name="logo"/>
<div align="center">
<img src="docs/Logo.png" height="130"></img>
</a>
</div>

Linux: [![Build Status](https://travis-ci.org/JuliaRheology/RHEOS.jl.svg?branch=master)](https://travis-ci.org/JuliaRheology/RHEOS.jl) &nbsp;
Windows: [![Build Status](https://ci.appveyor.com/api/projects/status/github/JuliaRheology/RHEOS.jl?branch=master&svg=true)](https://ci.appveyor.com/project/JuliaRheology/RHEOS-jl) &nbsp;
[![Latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://JuliaRheology.github.io/RHEOS.jl/latest) &nbsp;
[![License](https://img.shields.io/badge/License-MIT-ff69b2.svg?style=flat)](LICENSE.md) &nbsp;
[![status](https://joss.theoj.org/papers/553250d815e1990db1b89c742854c71a/status.svg)](https://joss.theoj.org/papers/553250d815e1990db1b89c742854c71a)

# RHEOS - RHEology, Open-Source
RHEOS, an abbreviation of Rheology Open Source, is a software package written in the Julia programming language that provides tools for analyzing rheological data. Features include:

- Stress/Strain/Time data can be easily be fitted to a viscoelastic model

- G'/G''/Frequency data can easily be fitted to a viscoelastic model

- Many standard and fractional viscoelastic models have already been implemented within RHEOS new ones can easily be added by users

- A fitted model can be used to predict the behaviour of the material under other loading conditions, enabling the fit/predict paradigm of model selection

- Artificial loading conditions can be generated within RHEOS to better understand a model's response

## Installation
1. Install Julia, version 1.2
2. From Julia REPL, enter pkg mode by pressing ```]```
3. (Optional) Enable desired Project.toml environment
4. Run the command ```add "https://github.com/JuliaRheology/RHEOS.jl"```

## Citation
If you use RHEOS in your work, please consider citing the following paper:
**TBA**

## Embedded Dependencies
#### [MittagLeffler.jl](https://github.com/jlapeyre/MittagLeffler.jl)

## Contributing to RHEOS
If you believe you have found any bugs or invalid behaviour in RHEOS, please feel free to file an issue on this repository. You can also raise an issue if you feel that any part of the documentation needs clarification, or for any feature requests. Even better than just raising an issue, you could both raise an issue and issue a pull request which fixes that issue. Note that meta-documentation on running tests and building documentation locally is available at the [JuliaRheology/RheoHelpDocs](https://github.com/JuliaRheology/RheoHelpDocs) repository. Please be aware that RHEOS is released with a [Contributor Code of Conduct](CONDUCT.md) and by participating in this project you agree to abide by its terms.

## References

+ W. N. Findley, J. S. Lai, K. Onaran — *Creep and Relaxation of Nonlinear Viscoelastic Materials (with an Introduction to Linear Viscoelasticity)*, Dover Publications, New York. (1989)

+ S. G. Johnson — *The NLopt nonlinear-optimization package*, https://github.com/stevengj/nlopt

+ J. Bezanson, A. Edelman, S. Karpinski, V. B. Shah — *Julia: A Fresh Approach to Numerical Computing*, SIAM Review, doi: 10.1137/141000671. (2017)
