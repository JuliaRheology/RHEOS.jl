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
*A suite of tools for analysing rheology data.*

- Stress/Strain/Time data can be easily be fitted to a viscoelastic model

- G'/G''/Frequency data can easily be fitted to a viscoelastic model

- Many standard and fractional viscoelastic models have already been implemented within RHEOS new ones can easily be added by users

- A fitted model can be used to predict the behaviour of the material under other loading conditions, enabling the fit/predict paradigm of model selection

- Artificial loading conditions can be generated within
RHEOS to better understand a model's response

## Installation

1. Install Julia, version 1.0.1
2. From Julia REPL, type ```Pkg.clone("https://github.com/JuliaRheology/RHEOS.jl.git")```
3. Run ```julia TEMP_INSTALL.jl```, a script located in your RHEOS directory

## To do
- [ ] Add FFT fitting to handle singularities and sidestep Mittag-Leffler bottleneck
- [ ] Documentation
- [ ] Add tests
- [ ] Add Sync Interpolation for going from variable to constant sample rate
- [ ] Set up OSX CI on Travis

## References & Included Dependencies
#### [FastConv.jl](https://github.com/aamini/FastConv.jl)
+ A. Alexander, B. Horn and A. Edelman - *Accelerated Convolutions for Efficient Multi-Scale Time to Contact Computation in Julia*, arXiv preprint arXiv:1612.08825 **(2016)**

#### [MittagLeffler.jl](https://github.com/jlapeyre/MittagLeffler.jl)
+ R. Gorenflo, J. Loutchko and Y. Loutchko - *Computation of the Mittag-Leffler function and its derivative*,  Fract. Calc. Appl. Anal. **(2002)**
