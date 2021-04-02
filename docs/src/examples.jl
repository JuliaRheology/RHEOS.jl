# # Additional Examples
#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/examples.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/examples.ipynb)

using RHEOS

# [PyPlot](https://github.com/JuliaPy/PyPlot.jl) needs to be installed to run these examples and display plots of the data.

using PyPlot

# ## Example 1

# - Loading experimental data from a .csv file
# - Fitting a model
# - Compare original data with fitted model prediction

## Make sure the examples folder is the current directory
## check by typing "pwd()"

## Import data
data = importcsv("assets/example1_data.csv", t_col = 1, ϵ_col = 2, σ_col = 3)

## Plot data
fig, ax = subplots(1, 1, figsize = (7, 5))
ax.plot(data.t, data.σ, ".", color = "green")
ax.plot(data.t, data.ϵ, "-", color = "blue")
ax.set_ylabel("Strain (blue), Stress (green)")
ax.set_xlabel("Time")
#!nb fig #hide
#-

## We now fit a Maxwell model
maxwell_model = modelfit(data, Maxwell, strain_imposed)

# Note that the fitting function requires guidance regarding the type of testing used. It helps optimise the fitting process.

# The data in this example is the stress response to a strain ramp followed by plateau. It therefore corresponds to a strain imposed process.

# We now want to calculate the stress values predicted by the model given the experimental strain data. Let's create a new data set with the strain profile.

maxwell_predict = extract(data, strain_only)
## and calculate the stress based on the model
maxwell_predict = modelpredict(maxwell_predict, maxwell_model)
## Now we can plot data and model together for comparison

## Plot data
fig, ax = subplots(1, 1, figsize = (7, 5))
ax.plot(data.t, data.σ, ".", color = "green")
ax.plot(maxwell_predict.t, maxwell_predict.σ, color = "red")
ax.set_xlabel("Time")
ax.set_ylabel("Stress")
#!nb fig #hide

# ## Example 2

# This script is a slight modification of Example 1 to present to the user the possibility of creating new [`RheoModelClass`](@ref) from an existing one with some of the parameters frozen to specific values. As an example, we fix the spring constant of the model above (k) to 2 and we let RHEOS fit the viscosity η.

Maxwell_springFix = freeze_params(Maxwell, (k = 2, ))
#-

maxwellD_model = modelfit(data, Maxwell_springFix, strain_imposed)
maxwellD_predict = extract(data, strain_only)
## and calculate the stress based on the model
maxwellD_predict = modelpredict(maxwellD_predict, maxwellD_model)
## Now we can plot data and model together for comparison

## Plot data
fig, ax = subplots(1, 1, figsize = (7, 5))
ax.plot(data.t, data.σ, ".", color = "green")
ax.plot(maxwellD_predict.t, maxwellD_predict.σ, color = "red")
ax.set_xlabel("Time")
ax.set_ylabel("Stress")
#!nb fig #hide

# ## Example 3

# This script shows how to use RHEOS to explore the behaviour of various models This involves:
# - Creating a strain function
# - Defining models based on parameter values
# PyPlot needs to be installed to run these examples and display plots of the data.

## Creates a time only dataset
dϵ = timeline(t_end = 10)
## calculates strain data by applying a function of time
dϵ = strainfunction(dϵ, t -> sin(t))

## Plot strain data
fig, ax = subplots(1, 1, figsize = (7, 5))
ax.plot(dϵ.t, dϵ.ϵ, "--b")

## we can now simulate various models based on this strain only dataset
## Let's study the role of the dashpot strength in the MAxwell model
for η in [0.1, 0.3, 1, 3, 10]
    maxwell_model = RheoModel(Maxwell, k = 2., η = η)
    d_maxwell = modelpredict(dϵ, maxwell_model)
    ax.plot(d_maxwell.t, d_maxwell.σ)
end
ax.set_xlabel("Time")
ax.set_ylabel("Stress")
ax.grid("on")
#!nb fig #hide
