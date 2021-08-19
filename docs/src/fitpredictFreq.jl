# # Fitting and Predicting - Frequency Data
#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/fitpredictFreq.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/fitpredictFreq.ipynb)

# RHEOS can also fit models to dynamic mechanical analysis data from oscillatory tests.

using RHEOS

# ## Fitting

# #### Step 1: Loading Experimental Data

# RHEOS has a convenient function to import data from CSV files (more information in [File I/O](@ref) section).

data = importcsv("assets/data_freq.csv", ω_col = 1, Gp_col = 2, Gpp_col = 3)
rheofreqdatatype(data)

# The frequency, storage modulus and loss modulus are stored into a [`RheoFreqData`](@ref) struct.

# #### Step 2: Model Fit

# Let's fit a Fractional Kelvin-Voigt model. The first argument is our data, the second argument tells RHEOS which model to fit ([`RheoModelClass`](@ref)). The function will return a [`RheoModel`](@ref) object, i.e. a [`RheoModelClass`](@ref) with fixed values for the parameters. The function [`dynamicmodelfit`](@ref) possesses all the additional arguments (lower and upper bounds, initial parameters, tolerance) as for the fitting of time data. Note that initial parameters, upper bounds and lower bounds can be partially specified and the rest will be filled in automatically.

## Lower bounds
Lo = (cₐ = 0.0, a = 0.01, cᵦ = 0.0, β = 0.01)
## Upper bounds
Hi = (cₐ = 1e2, a = 0.99, cᵦ = 1e2, β = 0.99)
## Initial parameters
P0 = (cₐ = 1.0, a = 0.36, cᵦ = 1.0, β = 0.03)

FractKV_model = dynamicmodelfit(data, Fract_KelvinVoigt, weights = "log", lo = Lo, hi = Hi, p0 = P0)

# For the fitting process RHEOS relies on the optimistion package [NLopt.jl](https://nlopt.readthedocs.io/en/latest/). RHEOS makes use of a local derivative free algorithm, specifically the Tom Rowan's "Subplex" algorithm originally introduced in his [PhD thesis at The University of Texas at Austin, 1990](https://dl.acm.org/doi/book/10.5555/100816).

# The storage and loss moduli can sometimes occupy different orders of magnitude. This can cause problems during fitting as the optimisation routine will weight errors at the higher orders of magnitude more strongly than those at the lower orders of magnitude. This is a more general problem often faced during multi-objective optimisation problems. RHEOS offers four rescaling options (weight):
# - The `local` option rescales the cost at each point by the point itself.
# - The `log` approach simply re-scales all the storage and loss moduli and their predicted values logarithmically before finding the error between them.
# - A third option in RHEOS divides the cost at each point by the mean value of the storage and loss modulus respectively depending on which cost is being calculated. This can work well but performance is hindered if the storage or loss moduli individually vary over many orders of magnitude (`mean`).
# - The fourth option offered by RHEOS is simply manual weightings provided by the user for each modulus (`manual`).
# In the first two cases, the benefit arises from the fact that the optimisation weighting is rebalanced in favour of smaller values. The logarithmic rescaling method seems to work particularly well but runs into problems if the storage or loss moduli are exactly 0 at any frequency, this is due to the negative singularity of the logarithmic function for the 0 argument. 0 arguments can also cause problems using the `local` option due to division by 0.

# ## Predicting

# RHEOS allows the user to simulate the frequency response of a model (with defined parameters, [`RheoModel`](@ref) struct) to an imposed loading.

# Given an incomplete data set (frequency only) and model with values substituted into parameters ([`RheoModel`](@ref)), return a new "complete" dataset based on the model with the simulated missing variable (storage and loss).

# #### Assessing the Quality of Fit

# The ability of predicting model's response is first exploited to assess the quality of the fits above. The incomplete [`RheoFreqData`](@ref) variable is defined by extracting the frequency from the original data or by defining a new frequency vector.

data_ext = onlyfreq(data)
## Alternatively, a RheoFreqData with only frequency data can be generated as
## data_ext = frequencyspec(ω_start = 1.0e-2, ω_end = 1.0e2, logstep = 0.1)

rheofreqdatatype(data_ext)

# For the prediction, RHEOS' function [`dynamicmodelpredict`](@ref) requires the incomplete data set and a model with fixed parameters ([`RheoModel`](@ref)). For the assessment of the fitting quality the [`RheoModel`](@ref) is the output of the fitting function.

fractKV_predict = dynamicmodelpredict(data_ext, FractKV_model)

## Now we can plot data and model together for comparison
using PyPlot
fig, ax = subplots(1, 1, figsize = (7, 5))
ax.loglog(data.ω, data.Gp, "o", markersize = 5, color = "blue")
ax.loglog(data.ω, data.Gpp, "o", markersize = 5, color = "red")
ax.loglog(fractKV_predict.ω, fractKV_predict.Gp, "--", color = "blue")
ax.loglog(fractKV_predict.ω, fractKV_predict.Gpp, "--", color = "red")
ax.set_xlabel("Frequency")
ax.set_ylabel("Storage and Loss moduli")
#!nb fig #hide

# #### Simulate Different Model Behaviours

# The ability of predicting model's response can be exploited to simulate the behaviour of the model to other external loading conditions. For example, we will explore the response of the fitted model to a creep experiment.

# First we need to define an incomplete [`RheoTimeData`](@ref) struct, which can be achieved via the data generation functions provided in RHEOS (see section [Generating Data](@ref)).

## Define timeline
dσ = timeline(t_end = 10)
## and a step in stress
dσ = stressfunction(dσ, hstep())

## we can now predict the creep response of the Maxwell model 
FractKV_creepPredict = modelsteppredict(dσ, FractKV_model)
## Visualisation of the simulated response
fig, ax = subplots(1, 1, figsize = (7, 5))
ax.plot(FractKV_creepPredict.t, FractKV_creepPredict.ϵ)
ax.set_xlabel("Time")
ax.set_ylabel("Strain")
#!nb fig #hide

# **Reference frequency data**: Deng, Linhong, et al. "Fast and slow dynamics of the cytoskeleton." Nature materials 5.8 (2006): 636.
