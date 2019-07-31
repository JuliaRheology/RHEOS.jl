#=
Example 1

- Loading experimental data from csv file
- Fitting model
- Comparing original data with fit
=#
using RHEOS

# PyPlot needs to be installed to run these examples
# and display plots of the data.

using PyPlot

# Make sure the examples folder is the current directory
# check by typing "pwd()"

# Import data
data = importcsv("../examples/example1_data.csv", t_col=1, ϵ_col=2, σ_col=3)

# Plot data
plot(data.t,data.ϵ,"-")
plot(data.t,data.σ,".")

# We now fit a Standard Linear Solid model
maxwell_model = modelfit(data, Maxwell, strain_imposed)

#=
Note that the fitting function requires guidance regarding the type of testing used.
It helps optimise the fitting process.

The data in this example is the stress response to a strain ramp followed by plateau.
It therefore corresponds to a strain imposed process.
=#

# We now want to calculate the stress values predicted by the model given the experimental strain data.
# Lets create a new data set with the strain profile

maxwell_predict = extract(data, strain_only)
# and calculate the stress based on the model
maxwell_predict = modelpredict(maxwell_predict, maxwell_model)
# Now we can plot data and model together for comparison
plot(maxwell_predict.t,maxwell_predict.σ)
