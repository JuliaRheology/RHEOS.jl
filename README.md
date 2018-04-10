# RHEOS - RHEology, Open-Source
#### A suite of tools for analysing rheology data

*Documentation currently under development*

## Installation

- Install julia, version > 0.6.
- Clone repository into desired location
- Run TEMP_INSTALL.jl  
- Set environment variable "JULIA_NUM_THREADS" equal to number of processor cores available

## References & Included Dependencies
#### [FastConv.jl](https://github.com/aamini/FastConv.jl)
A. Alexander, B. Horn and A. Edelman - *Accelerated Convolutions for Efficient Multi-Scale Time to Contact Computation in Julia*, arXiv preprint arXiv:1612.08825 **(2016)**

#### [MittagLeffler.jl](https://github.com/jlapeyre/MittagLeffler.jl)
R. Gorenflo, J. Loutchko and Y. Loutchko - *Computation of the Mittag-Leffler function and its derivative*,  Fract. Calc. Appl. Anal. **(2002)**

#### Fractional Zener, Fractionl Maxwell and Fractional Kelvin-Voigt Models
F. Mainardi, G. Spada - *Creep, relaxation and viscosity properties for basic fractional models in rheology*, Eur. Phys. J. Spec. Top. **(2011)**

#### Viscoelastic Hertz Model of Spherical Indententation 
M. L. Oyen - [*Spherical Indentation Creep Following Ramp Loading*](https://doi.org/10.1557/JMR.2005.0259), Journal of Materials Research **(2005)**  
J. M. Mattice, A. G. Lau, M. L. Oyen and R. W. Kent - [*Spherical indentation load-relaxation of soft biological tissues*](https://doi.org/10.1557/jmr.2006.0243), Journal of Materials Research **(2006)**
