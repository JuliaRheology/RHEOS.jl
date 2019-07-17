#=
Example 2

This script shows how to use RHEOS to explore the behaviour of various modeks
This involves:

- Creating a strain function
- Defining models based on parameter values
=#

# PyPlot needs to be installed to run these examples
# and display plots of the data.

using PyPlot



# Creates a time only dataset
dϵ=timeline()
# calculates strain data by applying a function of time
dϵ=strainfunction(dϵ,t->sin(t))

# Plot strain data
plot(dϵ.t,dϵ.ϵ,"--b")


# we can now simulate various models based on this strain only dataset
# Let's study the role of the dashpot strength in the MAxwell model
for η in [0.1, 0.3, 1, 3, 10]
    maxwell_model = RheoModel(Maxwell, k = 2., η = η)
    d_maxwell = modelpredict(dϵ, maxwell_model)
    plot(d_maxwell.t,d_maxwell.σ)
end
