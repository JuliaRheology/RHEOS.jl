<a name="logo"/>
<div align="center">
<img src="docs/Logo.png" height="130"></img>
</a>
</div>

Linux: [![Build Status](https://travis-ci.org/JuliaRheology/RHEOS.jl.svg?branch=master)](https://travis-ci.org/JuliaRheology/RHEOS.jl) &nbsp;
Windows: [![Build Status](https://ci.appveyor.com/api/projects/status/github/JuliaRheology/RHEOS.jl?branch=master&svg=true)](https://ci.appveyor.com/project/JuliaRheology/RHEOS-jl) &nbsp;
[![Latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://JuliaRheology.github.io/RHEOS.jl/latest) &nbsp;
[![License](https://img.shields.io/badge/License-MIT-ff69b2.svg?style=flat)](LICENSE.md)

# RHEOS - RHEology, Open-Source
RHEOS, an abbreviation of Rheology Open Source, is a software package written in the Julia programming language that provides tools for analyzing rheological data. Features include:

- Stress/Strain/Time data can be easily be fitted to a viscoelastic model

- G'/G''/Frequency data can easily be fitted to a viscoelastic model

- Many standard and fractional viscoelastic models have already been implemented within RHEOS new ones can easily be added by users

- A fitted model can be used to predict the behaviour of the material under other loading conditions, enabling the fit/predict paradigm of model selection

- Artificial loading conditions can be generated within RHEOS to better understand a model's response

## Code Example
The code below demonstrates how simple it is to load in data, fit a model and then predict its response using RHEOS.
```
# load in data from CSV using RHEOS convenience function
data = importdata(["stress","strain", "time"], "DataRelaxation.csv")

# fit a Standard Linear Solid model to the data via its relaxation modulus, G
fitted_SLS_model = modelfit(data, SLS(), :G)

# based on the original strain data and the SLS model fitted in the line above
# predict the stress using the SLS relaxation modulus
SLS_predicted = modelpredict(data, fitted_SLS_model, :G)
```

## Documentation

RHEOS is fully documented with comprehensive tutorials that assume minimal knowledge about Julia programming. The documentation can be accessed by [clicking here](https://juliarheology.github.io/RHEOS.jl/latest/) or the blue 'docs' badge near the top of this page (just below the RHEOS logo).

## Installation
1. Install Julia, version 1.0.2
2. From Julia REPL, enter pkg mode by pressing ```]```
3. (Optional) Enable desired Project.toml environment
4. Run the command ```add "https://github.com/JuliaRheology/RHEOS.jl"```

## To do
- [ ] Increase Test Coverage
- [ ] Implement derived singularity approximation
- [ ] Remove "sampling" field from RheologyData, add tests where required
- [ ] Add default constraints for fractional models for convenience in broad parameter searches
- [ ] Add FFT fitting to handle singularities and sidestep Mittag-Leffler bottleneck
- [ ] Test integration-by-parts singularity handling
- [ ] Add Sync Interpolation for going from variable to constant sample rate

## Included Dependencies
#### [FastConv.jl](https://github.com/aamini/FastConv.jl)

#### [MittagLeffler.jl](https://github.com/jlapeyre/MittagLeffler.jl)

## Citation
If you use RHEOS in your work, please consider citing the following paper:
**TBA**

## References

+ W. N. Findley, J. S. Lai, K. Onaran — *Creep and Relaxation of Nonlinear Viscoelastic Materials (with an Introduction to Linear Viscoelasticity)*, Dover Publications, New York. (1976)

+ S. G. Johnson — *The NLopt nonlinear-optimization package*, http://ab-initio.mit.edu/nlopt

+ J. Bezanson, A. Edelman, S. Karpinski, V. B. Shah — *Julia: A Fresh Approach to Numerical Computing*, SIAM Review, doi: 10.1137/141000671. (2017)
