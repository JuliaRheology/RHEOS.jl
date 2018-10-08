# RHEOS.jl 

## Overview

RHEOS, an abbreviation of Rheology Open Source, is a software package written in the Julia programming language that provides tools for analyzing rheological data. Features include:

- G'/G''/Frequency data can easily be fitted to a viscoelastic model

- Many standard and fractional viscoelastic models have already been implemented within RHEOS new ones can easily be added by users

- A fitted model can be used to predict the behaviour of the material under other loading conditions, enabling the fit/predict paradigm of model selection

- Artificial loading conditions can be generated within
RHEOS to better understand a model's response

## Installation

1. Install Julia, version 1.0.1
2. From Julia REPL, enter pkg mode by pressing ```]```
3. (Optional) Enable desired Project.toml environment
4. Run the command ```add "https://github.com/JuliaRheology/RHEOS.jl"```

## Included Dependencies
#### [FastConv.jl](https://github.com/aamini/FastConv.jl)

#### [MittagLeffler.jl](https://github.com/jlapeyre/MittagLeffler.jl)

## References
#### If you use RHEOS in your work, please consider citing the following papers

+ S. G. Johnson -  *The NLopt nonlinear-optimization package*, http://ab-initio.mit.edu/nlopt

+ R. Gorenflo, J. Loutchko and Y. Loutchko - *Computation of the Mittag-Leffler function and its derivative*,  Fract. Calc. Appl. Anal. **(2002)**

+ A. Alexander, B. Horn and A. Edelman - *Accelerated Convolutions for Efficient Multi-Scale Time to Contact Computation in Julia*, arXiv preprint arXiv:1612.08825 **(2016)**