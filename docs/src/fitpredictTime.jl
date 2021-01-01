# # Fitting and Predicting - Time Data
#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/fitpredictTime.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/fitpredictTime.ipynb)

# A large majority of scientists and engineers who undertake rheological experiments fit their data with one or several viscoelastic models in order to classify materials, quantify their behaviour and predict their response to external perturbations.

# Standard linear viscoelastic models take the form of an ordinary differential equation between stress σ and strain ϵ. Fitting models and predicting their response in the time domain then requires computing viscoelastic hereditary integrals such as:

# ```math
# \sigma(t) = \int_{0}^t G(t - \tau) \frac{d \epsilon(\tau)}{d \tau} d \tau
# ```

# where G is the relaxation response (stress response to a step in strain) of the material, which can be analytically or numerically defined for each rheological model and its set of rheological parameters. A similar relation exists to calculate the strain from the stress history:

# ```math
# \epsilon(t) = \int_{0}^t J(t - \tau) \frac{d \sigma(\tau)}{d \tau} d \tau
# ```

# where J is the creep modulus of the material. Fitting and predicting behaviour then becomes non-trivial and standardised tools to process are needed.

using RHEOS

# ## Fitting

# #### Step 1: Loading in Experimental Data

# RHEOS has a convenient function to import data from CSV files (more information in [File I/O](@ref) section).

data = importcsv("assets/data_time.csv", t_col = 1, ϵ_col = 2, σ_col = 3)
RheoTimeDataType(data)

# The time, strain and stress data are stored into a [`RheoTimeData`](@ref) struct.

# #### Step 2: Fitting Models

# Let's fit a Maxwell model via its relaxation modulus, G, as our data is from a stress relaxation test, thus the strain is imposed. The first argument is our data, the second argument tells RHEOS which model to fit ([`RheoModelClass`](@ref)) and the final argument tells RHEOS whether to fit the model using a relaxation modulus (*strain* imposed) or creep modulus (*stress* imposed). The function will return a [`RheoModel`](@ref) object, i.e. a [`RheoModelClass`](@ref) with fixed values for the parameters.

maxwell_model = modelfit(data, Maxwell, strain_imposed)

# For the fitting process RHEOS relies on the optimistion package [NLopt.jl](https://nlopt.readthedocs.io/en/latest/). RHEOS makes use of a local derivative free algorithm, specifically the Tom Rowan's "Subplex" algorithm originally introduced in his [PhD thesis at The University of Texas at Austin, 1990](https://dl.acm.org/doi/book/10.5555/100816).

# Next, we will fit a fractional Maxwell (spring) model (the only difference from the above model is that the dash-pot is replaced by a spring-pot). This time we will also add upper and lower bounds on the model parameters. This is highly recommended for fractional models in particular as values less than 0 or greater than 1 for the spring-pot parameter are unphysical and can cause errors in the Mittag-Leffler function used (more information about the fractional models in review TBA). To facilitate the convergence of the fitting algorithm, it is also suggested to define an initial set of parameters (the function uses 0.5 for all parameters not defined). For initial parameters, lower bounds and upper bounds, it is possible to only specify some of the parameters and the rest will be filled in automatically.

## Lower bounds
Lo = (cₐ = 0.0, a = 0.01, k = 0.0)
## Upper bounds
Hi = (cₐ = Inf, a = 0.99, k = Inf)
## Initial parameters
P0 = (cₐ = 0.1, a = 0.2, k = 2.0)

Fmaxwell_model = modelfit(data, FractS_Maxwell, strain_imposed, lo = Lo, hi = Hi, p0 = P0)

# Further optional parameters can be defined in the fitting function, such as the relative tolerance (default parameter of 1e-4) or the finite difference formula to use for the derivative (`"BD"` = Backward Difference (default), `"CD"` = Central Difference). If the user is interested in the model parameters on each optimisation iteration, they can be printed on the terminal by enabling `verbose = true`.

#md # !!! warning "Important"
#md #     Often an ideal step loading is assumed, which simplifies the equations in the first section to σ(t) = G(t)ϵ₀ and ϵ(t) = J(t)σ₀. This enables the modulus to be used directly for the fitting, which vastly reduces the computational burden. For this reason, RHEOS provides a dedicated function [`modelstepfit`](@ref) that requires the same parameters as [`modelfit`](@ref). If this assumption is appropriate for the data then fitting can be sped up greatly by use of this function.

#nb # Important Note: Often an ideal step loading is assumed, which simplifies the equations in the first section to σ(t) = G(t)ϵ₀ and ϵ(t) = J(t)σ₀. This enables the modulus to be used directly for the fitting, which vastly reduces the computational burden. For this reason, RHEOS provides a dedicated function `modelstepfit` that requires the same parameters as `modelfit`. If this assumption is appropriate for the data then fitting can be sped up greatly by use of this function.

# ## Predicting

# RHEOS allows the user to simulate the response of a model (with defined parameters, [`RheoModel`](@ref) struct) to an imposed loading.

# Given an incomplete data set (time data with only either stress or strain missing) and model with values substituted into parameters ([`RheoModel`](@ref)), return a new "complete" dataset based on the model with the simulated missing variable.

# #### Assess Quality of the Fit

# The ability of predicting model's response is first exploited to assess the quality of the fits above. The incomplete [`RheoTimeData`](@ref) variable is defined by extracting the time and the imposed variable (for the current example, the strain).

data_ext = extract(data, strain_only)
## alternatively: data_ext = extract(data, 1);

RheoTimeDataType(data_ext)

# For the prediction, RHEOS' function [`modelpredict`](@ref) requires the incomplete data set and a model with fixed parameters ([`RheoModel`](@ref)). For the assessment of the fitting quality the [`RheoModel`](@ref) required is the output of the fitting function.

maxwell_predict = modelpredict(data_ext, maxwell_model)

## Now we can plot data and model together for comparison
using PyPlot

fig, ax = subplots(1, 1, figsize = (7, 5))
ax.plot(data.t, data.σ, "o", markersize = 5)
ax.plot(maxwell_predict.t, maxwell_predict.σ, color = "red")
ax.set_xlabel("Time")
ax.set_ylabel("Stress")
fig

# #### Simulate different model behaviours

# The ability of predicting model's response can be exploited to simulate the behaviour of the model to other external loading conditions. For example, we will explore the response of the fitted model to a creep experiment.

# First we need to define an incomplete [`RheoTimeData`](@ref) struct, which can be achieved via the data generation functions provided in RHEOS (see section [Generating Data](@ref)).

## Define timeline
dσ = timeline(t_end = 10)
## and a step in stress
dσ = stressfunction(dσ, hstep())

## we can now predict the creep response of the Maxwell model 
maxwell_creepPredict = modelsteppredict(dσ, maxwell_model)

## Init plot figure
fig, ax = subplots(1, 1, figsize = (7, 5))

## Visualisation of the imposed loading
ax.plot(maxwell_creepPredict.t, maxwell_creepPredict.σ, "--")
## Visualisation of the simulated response
ax.plot(maxwell_creepPredict.t, maxwell_creepPredict.ϵ)
ax.set_xlabel("Time")
ax.set_ylabel("Strain")
fig

#md # !!! warning "Important"
#md #     As for the fitting procedure, if an ideal step loading is assumed, RHEOS provides a dedicated function [`modelsteppredict`](@ref) that requires the same parameters as [`modelpredict`](@ref) which significantly reduces the computation time.

#nb # Important Note: As for the fitting procedure, if an ideal step loading is assumed, RHEOS provides a dedicated function modelsteppredict that requires the same parameters as modelpredict which significantly reduces the computating time.

# Similarly, more complex behaviour can be simulated (such as a stairs loading).

## Define timeline
dσ = timeline(t_end = 10)
## and a step in stress
dσ = stressfunction(dσ, stairs())

## we can now predict the stairs response of the Maxwell model 
maxwell_creepPredict = modelpredict(dσ, maxwell_model)
## Plotting
fig, ax = subplots(1, 1, figsize = (7, 5))
ax.plot(maxwell_creepPredict.t, maxwell_creepPredict.σ, "--")
ax.plot(maxwell_creepPredict.t, maxwell_creepPredict.ϵ)
ax.set_xlabel("Time")
ax.set_ylabel("Strain")
fig
