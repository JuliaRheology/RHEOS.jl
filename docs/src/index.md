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

