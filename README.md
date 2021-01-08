<a name="logo"/>
<div align="center">
<img src="https://raw.githubusercontent.com/JuliaRheology/RHEOS.jl/master/docs/src/assets/logo.png" height="130"></img>
</a>
</div>
<!-- delim -->

# RHEOS - RHEology, Open-Source

|**Tests**|**Documentation**|**Coverage**|**Paper**|
|:-------------:|:---------------:|:---------------:|:-------:|
| [![Build Status](https://github.com/JuliaRheology/RHEOS.jl/workflows/stable/badge.svg)](https://github.com/JuliaRheology/RHEOS.jl/actions?query=workflow%3Astable) [![Build Status](https://github.com/JuliaRheology/RHEOS.jl/workflows/nightly/badge.svg)](https://github.com/JuliaRheology/RHEOS.jl/actions?query=workflow%3Anightly) | [![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaRheology.github.io/RHEOS.jl/stable) [![Development Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaRheology.github.io/RHEOS.jl/dev) | [![Code Coverage](https://codecov.io/gh/JuliaRheology/RHEOS.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaRheology/RHEOS.jl) | [![status](https://joss.theoj.org/papers/553250d815e1990db1b89c742854c71a/status.svg)](https://joss.theoj.org/papers/553250d815e1990db1b89c742854c71a) |

RHEOS, an abbreviation of Rheology Open Source, is a software package written in the Julia programming language that provides tools for analyzing rheological data. Features include:

- Stress/Strain/Time data can be easily be fitted to a viscoelastic model

- G'/G''/Frequency data can easily be fitted to a viscoelastic model

- Many standard and fractional viscoelastic models have already been implemented within RHEOS new ones can easily be added by users

- A fitted model can be used to predict the behaviour of the material under other loading conditions, enabling the fit/predict paradigm of model selection

- Artificial loading conditions can be generated within RHEOS to better understand a model's response

## Statement of Need

### Arbitrary stress-strain curves and broad relaxation spectra require advanced software

A large majority of scientists and engineers who undertake rheological experiments fit their data with one or several viscoelastic models in order to classify materials, quantify their behaviour and predict their response to external perturbations.

### Learning about rheology is facilitated by the ability to explore a large database of models

Obtaining intuition for fractional viscoelastic theory can be difficult and learning material is sparse. Tools are needed to support researchers with their exploration of standard and advanced models and how they behave in response to idealised loading conditions, in particular when analytical expressions are difficult to obtain.


### Extracting parameters and comparing models and systems require standardised tools

Because understanding of materials is often dependent on summarising their behaviour with a model, one must be able to test and compare a broad range of models to inform model selection and reliably identify material parameters. There are currently very limited options available in the public domain, and most research groups have to invest significant effort developing custom software. An open-source standardised library of models and fitting algorithms would support the rheology research community and make analysis more systematic, transparent and reproducible.

## Features

RHEOS addresses the issues outlined in the Statement of Need in several ways.

- As well as being able to fit and predict assuming step loading of stress or strain, RHEOS can handle arbitrary loading for non-singular and singular models, and for constant or variable sample rates.

- RHEOS includes an extensive library of both traditional and fractional viscoelastic models. Although this library will satisfy most users, it is also straightforward to add additional models to RHEOS should they need to.

- For intuition-building and model exploration, RHEOS includes signal generation features so that common loading patterns (e.g. step, ramp, stairs) can be applied to unfamiliar models.

- As a convenience to the user, RHEOS also includes easy-to-use CSV importing and exporting functions, as well as a number of preprocessing functions for resampling and smoothing.

All of the above features are linked together in a seamless interface intended to be very approachable for less experienced programmers. The different paradigms of creep, relaxation and oscillatory testing are all accounted for, and models fitted against one type of data can be used to predict against a different type of data. (For instance, fitting against relaxation data and predicting the frequency response spectrum.)

## Documentation
The documentation (accessible via the blue docs button) aims to provide tutorials in different elements of RHEOS. The API section is a comprehensive list of RHEOS types and functions, and brief descriptions of their use. For corrections or further questions, please note the 'Contributing to RHEOS' section below and create an issue on the [GitHub repository](https://github.com/JuliaRheology/RHEOS.jl). **Note that whenever you restart your Julia session you will have to reload RHEOS by typing `using RHEOS`, to avoid repetition this line is not included in every piece of example code.**  Any plotting library can be used but the tutorials within this documentation uses the PyPlot Julia package.

## Installation
1. Install the latest version of Julia
2. From Julia interactive command-line REPL, enter pkg mode by pressing `]`
3. (Optional) Enable desired Project.toml environment
4. Run the command `add RHEOS`

## Documentation
If you installed RHEOS using the instructions above then you will have the _latest stable release_ of RHEOS; to access the documentation for this version [click here][docs-sta-url] or the blue `docs/stable` badge at the top of this README page. To access the _latest documentation built directly from the master branch_ [click here][docs-dev-url] or the blue `docs/dev` badge at the top of this README page. 

## Citation
If you use RHEOS in your work, please consider citing the following papers as appropriate:

+ J. L. Kaplan, A. Bonfanti, A. J. Kabla (2019). _RHEOS.jl -- A Julia Package for Rheology Data Analysis_. Journal of Open Source Software, 4(41), 1700, [https://doi.org/10.21105/joss.01700](https://doi.org/10.21105/joss.01700)

+ A. Bonfanti, J. L. Kaplan, G. Charras, A. J. Kabla (2020). _Fractional viscoelastic models for power-law materials_. Soft Matter, 16, 6002-6020, [https://doi.org/10.1039/D0SM00354A](https://doi.org/10.1039/D0SM00354A)

## Embedded Dependencies
#### [MittagLeffler.jl](https://github.com/jlapeyre/MittagLeffler.jl)

## Contributing to RHEOS
If you believe you have found any bugs or invalid behaviour in RHEOS, please feel free to file an issue on this repository. You can also raise an issue if you feel that any part of the documentation needs clarification, or for any feature requests. Even better than just raising an issue, you could both raise an issue and issue a pull request which fixes that issue.

We generally work using a greatly simplified version of 'git flow'. Small changes and fixes are directly committed or PR'ed in to the `master` branch, as longs as tests pass. Larger developments are generally made in feature branches and then merged to `master` via a pull request, once tests are passing. Stable releases are tagged and sent to the Julia central package repository via JuliaRegistrator and TagBot.

Note that meta-documentation on running tests and building documentation locally is available at the [JuliaRheology/RheoHelpDocs](https://github.com/JuliaRheology/RheoHelpDocs) repository. Please be aware that RHEOS is released with a 'Contributor Code of Conduct' (CONDUCT.md) and by participating in this project you agree to abide by its terms.

## References
+ W. N. Findley, J. S. Lai, K. Onaran (1989). *Creep and Relaxation of Nonlinear Viscoelastic Materials (with an Introduction to Linear Viscoelasticity)*, Dover Publications, New York. 

+ S. G. Johnson. *The NLopt nonlinear-optimization package*, https://github.com/stevengj/nlopt

+ J. Bezanson, A. Edelman, S. Karpinski, V. B. Shah (2017). *Julia: A Fresh Approach to Numerical Computing*, SIAM Review, doi: 10.1137/141000671.