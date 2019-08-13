using RHEOS; using PyPlot

data = importcsv("../examples/example1_data.csv", t_col=1, ϵ_col=2, σ_col=3)

maxwell_model = modelfit(data, Maxwell, strain_imposed)
maxwell_predict = modelpredict(extract(data, strain_only), maxwell_model)

fig, ax = subplots()

ax.plot(data.t, data.ϵ, "--", label="Strain (Data)")
ax.plot(data.t, data.σ, ".", label="Stress (Data)")
ax.plot(maxwell_predict.t, maxwell_predict.σ, label="Stress (Predicted)")

ax.set_xlabel("Time (s)")
ax.set_ylabel("Stress, Strain")
ax.legend()
