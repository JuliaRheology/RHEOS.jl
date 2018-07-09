include("../src/RHEOS.jl")
using RHEOS
using PyPlot

step = stepgen(50.0, 25.0; stepsize = 0.05, t_trans = 0.051)
repeatedstep = repeatdata(step, 5; t_trans = 0.051)

data = RheologyData(repeatedstep)

fractzener_model = RheologyModel(G_fractzener)

# numerical solution
fractzener_predicted = modelpredict(data, fractzener_model)

# analytical solution
tups = [25.0, 75.0, 125.0, 175.0, 225.0]
tdowns = [50.0, 100.0, 150.0, 200.0, 250.0]

analytic_data = zeros(length(data.t))
for i in 1:5
    params = [0.5, 1.0, 1.0, 0.5]

    up_tbegin_el = closestindex(data.t, tups[i])
    up = zeros(length(data.t))
    up[(up_tbegin_el+1):end] = G_fractzener(data.t[(up_tbegin_el+1):end] - tups[i], params)

    down_tbegin_el = closestindex(data.t, tdowns[i])
    down = zeros(length(data.t))
    down[(down_tbegin_el+1):end] = G_fractzener(data.t[(down_tbegin_el+1):end] - tdowns[i], params)

    analytic_data+=up
    analytic_data-=down

end

fig, ax = subplots()
ax[:plot](data.t, data.σ, "-", label="original")
ax[:plot](fractzener_predicted.t, fractzener_predicted.σ, "-", label="fract zener")
ax[:plot](data.t, analytic_data, "--", label="analytic", alpha=0.6)
ax[:legend](loc="best")
show()