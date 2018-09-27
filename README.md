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
A suite of tools for analysing rheology data. Stress/Strain/Time data can be easily be fitted to a viscoelastic model. (Commonly occuring
standard and fractional viscoelastic models have alread been implemented within RHEOS — new ones can be added by users.) A fitted model can
be used to predict the behaviour of the material under other loading conditions. Artificial loading conditions can be generated within
RHEOS to better understand models' response. 

## Installation

- Install Julia, version 0.6.3
- From Julia REPL, type ```Pkg.clone("https://github.com/JuliaRheology/RHEOS.jl.git")```
- Run ```julia TEMP_INSTALL.jl```, a script located in your RHEOS directory
- (optional) Set environment variable "JULIA_NUM_THREADS" equal to number of processor cores available

## To do
- [ ] Add frequency domain fitting optimization
- [ ] Documentation
- [ ] Add tests
- [ ] Add Sync Interpolation for going from variable to constant sample rate
- [ ] Set up OSX CI on Travis

## References & Included Dependencies
#### [FastConv.jl](https://github.com/aamini/FastConv.jl)
+ A. Alexander, B. Horn and A. Edelman - *Accelerated Convolutions for Efficient Multi-Scale Time to Contact Computation in Julia*, arXiv preprint arXiv:1612.08825 **(2016)**

#### [MittagLeffler.jl](https://github.com/jlapeyre/MittagLeffler.jl)
+ R. Gorenflo, J. Loutchko and Y. Loutchko - *Computation of the Mittag-Leffler function and its derivative*,  Fract. Calc. Appl. Anal. **(2002)**

#### [InverseLaplace.jl](https://github.com/jlapeyre/InverseLaplace.jl)
+ J.A.C Weideman, Algorithms for Parameter Selection in the Weeks Method for Inverting the Laplace Transform, SIAM Journal on Scientific Computing, Vol. 21, pp. 111-128 **(1999)**

+ Abate, J. and Valkó, P.P., Multi-precision Laplace transform inversion International Journal for Numerical Methods in Engineering, Vol. 60 (Iss. 5-7) pp 979–993 **(2004)**

+ Valkó, P.P. and Abate, J., Comparison of Sequence Accelerators for the Gaver Method of Numerical Laplace Transform Inversion, Computers and Mathematics with Application, Vol. 48 (Iss.3-40) pp. 629-636 **(2004)**

#### Fractional Zener, Fractional Maxwell and Fractional Kelvin-Voigt Models
F. Mainardi, G. Spada - [*Creep, relaxation and viscosity properties for basic fractional models in rheology*](https://doi.org/10.1140/epjst/e2011-01387-1), Eur. Phys. J. Spec. Top. **(2011)**
