# RHEOS.jl 

## Overview

RHEOS, an abbreviation of Rheology Open Source, is a software package written in the Julia programming language that provides tools for analyzing rheological data. Features include:

- Stress/Strain/Time data can be easily be fitted to a viscoelastic model

- G'/G''/Frequency data can easily be fitted to a viscoelastic model

- Many standard and fractional viscoelastic models have already been implemented within RHEOS new ones can easily be added by users

- A fitted model can be used to predict the behaviour of the material under other loading conditions, enabling the fit/predict paradigm of model selection

- Artificial loading conditions can be generated within RHEOS to better understand a model's response

## Documentation

The sections in this documentation each aim to provide tutorials in different elements of RHEOS. The [API](@ref) section is a comprehensive list of RHEOS types and functions brief descriptions of their use. For corrections or further questions, please create an issue on the [Github repository](https://github.com/JuliaRheology/RHEOS.jl).

## Installation

1. Install Julia, version 1.0.2
2. From Julia REPL, enter pkg mode by pressing ```]```
3. (Optional) Enable desired Project.toml environment
4. Run the command ```add "https://github.com/JuliaRheology/RHEOS.jl"```

## Included Dependencies
#### [FastConv.jl](https://github.com/aamini/FastConv.jl)

#### [MittagLeffler.jl](https://github.com/jlapeyre/MittagLeffler.jl)

## Citation
If you use RHEOS in your work, please consider citing the following paper
*TBA*

## References

+ W. N. Findley, J. S. Lai, K. Onaran — *Creep and Relaxation of Nonlinear Viscoelastic Materials (with an Introduction to Linear Viscoelasticity)*, Dover Publications, New York. (1976)

+ S. G. Johnson — *The NLopt nonlinear-optimization package*, http://ab-initio.mit.edu/nlopt

+ J. Bezanson, A. Edelman, S. Karpinski, V. B. Shah — *Julia: A Fresh Approach to Numerical Computing*, SIAM Review, doi: 10.1137/141000671. (2017)
