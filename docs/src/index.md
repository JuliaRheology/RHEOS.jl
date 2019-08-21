# RHEOS - RHEology Open Source

## Summary

Rheology is the science of deformation and flow, with a focus on materials that do not exhibit simple linear elastic or viscous Newtonian behaviours. Rheology plays an important role in the empirical characterisation of soft viscoelastic materials commonly found in the food and cosmetics industry, as well as in biology and bioengineering. A broad range of theoretical tools exist to extract material parameters and interpret them thanks to data analysis and/or physical modelling.

RHEOS (RHEology, Open-Source) is a software package designed to make the analysis of rheological data simpler, faster and more reproducible. RHEOS has a particular emphasis on linear rheological models containing fractional derivatives which have demonstrable utility for the modelling of biological materials but have hitherto remained in relative obscurity -- possibly due to their mathematical and computational complexity. RHEOS is written in Julia, which greatly assists achievement of our aims as it provides excellent computational efficiency and approachable syntax. RHEOS is fully documented and has extensive testing coverage.

It should be noted that RHEOS is not an optimisation package. It builds on another optimisation package, NLopt, by adding a large number of abstractions and functionality specific to the exploration of viscoelastic data.

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
The sections in this documentation each aim to provide tutorials in different elements of RHEOS. The [API](@ref) section is a comprehensive list of RHEOS types and functions, and brief descriptions of their use. For corrections or further questions, please create an issue on the [Github repository](https://github.com/JuliaRheology/RHEOS.jl). **Note that whenever you restart your Julia session you will have to reload RHEOS by typing `using RHEOS`, to avoid repetition this line is not included in every piece of example code.**  Any plotting library can be used but the tutorials within this documentation uses the PyPlot Julia package.

**Note that every section in this documentation is a Jupyter notebook (located in RHEOS/examples/) that can be locally modified to get familiar and to explore the functionalities offered by RHEOS**

## Installation
1. Install Julia, version 1.1.1
2. From Julia REPL, enter pkg mode by pressing ```]```
3. (Optional) Enable desired Project.toml environment
4. Run the command ```add "https://github.com/JuliaRheology/RHEOS.jl"```

## Included Dependencies
#### [MittagLeffler.jl](https://github.com/jlapeyre/MittagLeffler.jl)

## Citation
If you use RHEOS in your work, please consider citing the following paper
*TBA*

## References
+ W. N. Findley, J. S. Lai, K. Onaran — *Creep and Relaxation of Nonlinear Viscoelastic Materials (with an Introduction to Linear Viscoelasticity)*, Dover Publications, New York. (1976)

+ S. G. Johnson — *The NLopt nonlinear-optimization package*, http://ab-initio.mit.edu/nlopt

+ J. Bezanson, A. Edelman, S. Karpinski, V. B. Shah — *Julia: A Fresh Approach to Numerical Computing*, SIAM Review, doi: 10.1137/141000671. (2017)
