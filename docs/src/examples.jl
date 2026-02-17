# # Additional Examples
#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/examples.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/examples.ipynb)

using Plots
using RHEOS


# ## Example 1

# - Loading experimental data from a .csv file
# - Fitting a model
# - Compare original data with fitted model prediction

## Make sure the examples folder is the current directory
## check by typing "pwd()"

## Import data
data = importcsv("assets/example1_data.csv", t_col = 1, ϵ_col = 2, σ_col = 3)

## Plot data
p = plot(data.t, data.σ,
    seriestype = :scatter,
    color = :green,
    label = "Stress",
    xlabel = "Time",
    ylabel = "Strain (blue), Stress (green)",
    legend = :topright,
    markersize = 4, 
    framestyle = :box
)
plot!(p, data.t, data.ϵ, color = :blue, label = "Strain", linewidth = 3)
#!nb p #hide
#-

## We now fit a Maxwell model
maxwell_model = modelfit(data, Maxwell, strain_imposed)

# Note that the fitting function requires guidance regarding the type of testing used. It helps optimise the fitting process.

# The data in this example is the stress response to a strain ramp followed by plateau. It therefore corresponds to a strain imposed process.

# We now want to calculate the stress values predicted by the model given the experimental strain data. Let's create a new data set with the strain profile.

maxwell_predict = onlystrain(data)
## and calculate the stress based on the model
maxwell_predict = modelpredict(maxwell_predict, maxwell_model)
## Now we can plot data and model together for comparison

## Plot data
p = plot(data.t, data.σ,
    seriestype = :scatter, 
    color = :green,
    label = "Data",
    xlabel = "Time",
    ylabel = "Stress",
    legend = :topright,
    markersize = 4, 
    framestyle = :box
)

plot!(p, maxwell_predict.t, maxwell_predict.σ, color = :red, label = "Maxwell fit", linewidth = 3)
#!nb p #hide

# ## Example 2

# This script is a slight modification of Example 1 to present to the user the possibility of creating new [`RheoModelClass`](@ref) from an existing one with some of the parameters frozen to specific values. As an example, we fix the spring constant of the model above (k) to 2 and we let RHEOS fit the viscosity η.

Maxwell_springFix = freeze_params(Maxwell, k = 2)
#-

maxwellD_model = modelfit(data, Maxwell_springFix, strain_imposed)
maxwellD_predict = onlystrain(data)
## and calculate the stress based on the model
maxwellD_predict = modelpredict(maxwellD_predict, maxwellD_model)
## Now we can plot data and model together for comparison

## Plot data
p = plot(data.t, data.σ,
    seriestype = :scatter,      
    color = :green,
    label = "Data",
    xlabel = "Time",
    ylabel = "Stress",
    legend = :topright,
    markersize = 4, 
    framestyle = :box
)

plot!(p, maxwellD_predict.t, maxwellD_predict.σ, color = :red, label = "MaxwellD fit", linewidth = 3)
#!nb p #hide

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
p = plot(dϵ.t, dϵ.ϵ,
    linestyle = :dash,
    color = :blue,
    label = "Input strain",
    xlabel = "Time",
    ylabel = "Stress",
    linewidth = 3,
    legend = :bottomright, 
    framestyle = :box
)

## we can now simulate various models based on this strain only dataset
## Let's study the role of the dashpot strength in the MAxwell model
for η in [0.1, 0.3, 1, 3, 10]
    maxwell_model = RheoModel(Maxwell, k = 2.0, η = η)
    d_maxwell = modelpredict(dϵ, maxwell_model)
    plot!(p, d_maxwell.t, d_maxwell.σ, label = "η = $η", linewidth = 2)
end
plot!(p, grid = true)
#!nb p #hide

# ## Example 4

# This example demonstrates generating a timeline and stress data, fitting multiple models to the data,
# calling the `extractfitdata` function, listing the errors, and determining which model fits the best.


# Generate Timeline
datat = timeline(t_start = 0, t_end = 20.0, step = 0.02)  # Create a timeline from 0 to 20 seconds with a step size of 0.02 seconds

# Generate Stress Data (Ramp & hold)
dramp_stress = stressfunction(datat, ramp(offset = 4.0, gradient = 0.8))  # Generate a ramp stress function with offset 4.0 and gradient 0.8
dhold_stress = dramp_stress - stressfunction(datat, ramp(offset = 5.0, gradient = 0.8))  # Generate a hold stress function by subtracting a shifted ramp

# Define the rheological model and predict
model = RheoModel(SLS_Zener, (η = 1, kᵦ = 1, kᵧ = 1))
data = modelpredict(dhold_stress, model)

# Fit three models to the data
SLS_Zener_model = modelfit(data, SLS_Zener, stress_imposed)
Maxwell_model = modelfit(data, Maxwell, stress_imposed)
BurgersLiquid_model = modelfit(data, BurgersLiquid, stress_imposed)

# Call the extractfitdata function to extract fitting data
extracted_data = extractfitdata(data)

# Determine which model fits best by comparing errors
function find_best_model(extracted_data)
    best_model = ""
    min_error = Inf

    for (model_name, fitdata) in extracted_data
        err = fitdata[1].info.error
        println("Model: $model_name, Total Error: $err")

        if err < min_error
            min_error = err
            best_model = model_name
        end
    end

    println("Best fitting model: $best_model with total error: $min_error")
    return best_model, min_error
end

best_model, min_error = find_best_model(extracted_data)

# Create strain-only data for model predictions
stress_only_data = onlystress(data)
# Get model predictions for plotting
SLS_Zener_predict = modelpredict(stress_only_data, SLS_Zener_model)
Maxwell_predict = modelpredict(stress_only_data, Maxwell_model)
BurgersLiquid_predict = modelpredict(stress_only_data, BurgersLiquid_model)
# Plot data and fitted models
p = scatter(data.t, data.ϵ,
    color = :green,
    label = "Original Data",
    xlabel = "Time",
    ylabel = "Strain",
    legend = :bottomright, 
    markersize = 5, 
    markerstrokewidth = 1,
    grid = true, 
    framestyle = :box
)
## Overlay fitted models
plot!(p, SLS_Zener_predict.t, SLS_Zener_predict.ϵ, color = :red, linestyle = :solid, linewidth = 3, label = "SLS_Zener Model")
plot!(p, Maxwell_predict.t, Maxwell_predict.ϵ, color = :blue, linestyle = :dash, linewidth = 3, label = "Maxwell Model")
plot!(p, BurgersLiquid_predict.t, BurgersLiquid_predict.ϵ, color = :purple, linestyle = :dot, linewidth = 3, label = "BurgersLiquid Model")

#!nb p #hide
