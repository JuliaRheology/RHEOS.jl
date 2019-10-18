# # using RHEOS; using PyPlot; using PyCall; @pyimport seaborn as sns;

# # initialise seaborn package with correct parameters
# sns.set(rc=Dict("figure.figsize" => (10.5, 9.5)))
# sns.set_style("ticks")

# # initialise data
# ω_step = 0.01
# ω₀ = 0.0
# ω₊ = 10000.0
# chirp = collect(ω₀:ω_step:ω₊)

# tstep = 0.001
# t₀ = 0.0
# tᵢ = 10.0
# t = collect(t₀:tstep:tᵢ)

# models
# fkv_model = FractionalKelvinVoigt([1.0, 1.0, 1.0, 0.0])
# fkvs_model = FractionalKVspring([1.0, 1.0, 1.0])
# fkvd_model = FractionalKVdashpot([1.0, 1.0, 0.0])
# kelvinvoigt = KelvinVoigt([1.0, 1.0])

# G' for all models at α = 1.0, β = 0.0
# fig, ax = subplots()

# fkv_result = fkv_model.Gp(chirp, fkv_model.parameters)
# fkvs_result = fkvs_model.Gp(chirp, fkvs_model.parameters)
# fkvd_result = fkvd_model.Gp(chirp, fkvd_model.parameters)
# kelvinvoigt_result = kelvinvoigt.Gp(chirp, kelvinvoigt.parameters)

# ax[:loglog](chirp, kelvinvoigt_result, "-", color="black", linewidth=3, zorder=2)
# ax[:loglog](chirp, fkv_result, "--", color="red", linewidth=3, zorder=3)
# ax[:loglog](chirp, fkvs_result, "-.", color="green", linewidth=3, zorder=3)
# ax[:loglog](chirp, kelvinvoigt_result, "o", alpha=0.3, color="blue", linewidth=3, zorder=1)

# ax[:tick_params]("both", labelsize=24)

# ax[:grid]("on", alpha=0.4)
# show()

# G''
# fig, ax = subplots()

# fkv_result = fkv_model.Gpp(chirp, fkv_model.parameters)
# fkvs_result = fkvs_model.Gpp(chirp, fkvs_model.parameters)
# fkvd_result = fkvd_model.Gpp(chirp, fkvd_model.parameters)
# kelvinvoigt_result = kelvinvoigt.Gpp(chirp, kelvinvoigt.parameters)

# ax[:semilogx](chirp, kelvinvoigt_result, "-", color="black", linewidth=3, zorder=2)
# ax[:semilogx](chirp, fkv_result, "--", color="red", linewidth=3, zorder=3)
# ax[:semilogx](chirp, fkvs_result, "-.", color="green", linewidth=3, zorder=4)
# ax[:semilogx](chirp, kelvinvoigt_result, "o", alpha=0.3, color="blue", linewidth=3, zorder=1)

# ax[:tick_params]("both", labelsize=24)
# ax[:grid]("on", alpha=0.4)
# show()

# G
# fig, ax = subplots()

# fkv_result = fkv_model.G(t, fkv_model.parameters)
# fkvs_result = fkvs_model.G(t, fkvs_model.parameters)
# fkvd_result = fkvd_model.G(t, fkvd_model.parameters)
# kelvinvoigt_result = kelvinvoigt.G(t, kelvinvoigt.parameters)

# ax[:plot](t, kelvinvoigt_result, "-", color="black", linewidth=3, zorder=2)
# ax[:plot](t, fkv_result, "--", color="red", linewidth=3, zorder=3)
# ax[:plot](t, fkvs_result, "-.", color="green", linewidth=3, zorder=3)
# ax[:plot](t, kelvinvoigt_result, "o", alpha=0.1, color="blue", linewidth=3, zorder=1)

# ax[:tick_params]("both", labelsize=24)
# ax[:grid]("on", alpha=0.4)
# show()

# fkv_model = FractionalKelvinVoigt([1.0, 0.99999, 1.0, 0.00001])
# fkvs_model = FractionalKVspring([1.0, 0.99999, 1.0])
# fkvd_model = FractionalKVdashpot([1.0, 1.0, 0.00001])
# kelvinvoigt = KelvinVoigt([1.0, 1.0])

# # J
# fig, ax = subplots()

# fkv_result = fkv_model.J(t, fkv_model.parameters)
# fkvs_result = fkvs_model.J(t, fkvs_model.parameters)
# fkvd_result = fkvd_model.J(t, fkvd_model.parameters)
# kelvinvoigt_result = kelvinvoigt.J(t, kelvinvoigt.parameters)

# ax[:plot](t, kelvinvoigt_result, "-", color="black", linewidth=3, zorder=2)
# ax[:plot](t, fkv_result, "--", color="red", linewidth=3, zorder=3)
# ax[:plot](t, fkvs_result, "-.", color="green", linewidth=3, zorder=3)
# ax[:plot](t, kelvinvoigt_result, "o", alpha=0.3, color="blue", linewidth=3, zorder=1)

# ax[:tick_params]("both", labelsize=24)
# ax[:grid]("on", alpha=0.4)
# show()

# using Revise; using RHEOS; using PyPlot; using PyCall; @pyimport seaborn as sns;

# # initialise seaborn package with correct parameters
# sns.set(rc=Dict("figure.figsize" => (10.5, 9.5)))
# sns.set_style("ticks")

# # initialise data
# ω_step = 0.01
# ω₀ = 0.0
# ω₊ = 10000.0
# chirp = collect(ω₀:ω_step:ω₊)

# tstep = 0.001
# t₀ = 0.0
# tᵢ = 10.0
# t = collect(t₀:tstep:tᵢ)

# # models for G', G''
# fraczen = FractionalZener([2.0, 1.0, 1.0, 0.0, 0.5, 0.0])
# fracSLS = FractionalSLS([2.0, 1.0, 1.0, 0.5])
# fracSpec = FractionalSpecial([0.5, 1.0, 0.0, 2.0])
# stanSLS = SLS([0.5, 1.0, 2.0])

# # G' for all models
# fig, ax = subplots()

# fraczen_result = fraczen.Gp(chirp, fraczen.parameters)
# fracSLS_result = fracSLS.Gp(chirp, fracSLS.parameters)
# fracSpec_result = fracSpec.Gp(chirp, fracSpec.parameters)
# stanSLS_result = stanSLS.Gp(chirp, stanSLS.parameters)

# ax[:loglog](chirp, fraczen_result, "-", color="black", linewidth=3, zorder=2)
# ax[:loglog](chirp, fracSLS_result, "--", color="red", linewidth=3, zorder=3)
# ax[:loglog](chirp, fracSpec_result, "-.", color="green", linewidth=3, zorder=3)
# ax[:loglog](chirp, stanSLS_result, "o", alpha=0.3, color="blue", linewidth=3, zorder=1)

# ax[:tick_params]("both", labelsize=24)

# ax[:grid]("on", alpha=0.4)
# show()

# # G''
# fig, ax = subplots()

# fraczen_result = fraczen.Gpp(chirp, fraczen.parameters)
# fracSLS_result = fracSLS.Gpp(chirp, fracSLS.parameters)
# fracSpec_result = fracSpec.Gpp(chirp, fracSpec.parameters)
# stanSLS_result = stanSLS.Gpp(chirp, stanSLS.parameters)

# ax[:loglog](chirp, fraczen_result, "-", color="black", linewidth=3, zorder=2)
# ax[:loglog](chirp, fracSLS_result, "--", color="red", linewidth=3, zorder=3)
# ax[:loglog](chirp, fracSpec_result, "-.", color="green", linewidth=3, zorder=3)
# ax[:loglog](chirp, stanSLS_result, "o", alpha=0.3, color="blue", linewidth=3, zorder=1)

# ax[:tick_params]("both", labelsize=24)
# ax[:grid]("on", alpha=0.4)
# show()

# fraczen = FractionalZener([1.0, 0.999, 1.0, 0.001, 1.0, 0.001])
# fracSLS = FractionalSLS([1.0, 0.999, 1.0, 1.0])
# fracSpec = FractionalSpecial([1.0, 1.0, 0.001, 1.0])
# stanSLS = SLS([1.0, 1.0, 1.0])

# G
# fig, ax = subplots()

# fraczen_result = fraczen.G(t, fraczen.parameters)
# fracSLS_result = fracSLS.G(t, fracSLS.parameters)
# fracSpec_result = fracSpec.G(t, fracSpec.parameters)
# stanSLS_result = stanSLS.G(t, stanSLS.parameters)

# ax[:plot](t, fraczen_result, "-", color="black", linewidth=3, zorder=2)
# ax[:plot](t, fracSLS_result, "--", color="red", linewidth=3, zorder=3)
# ax[:plot](t, fracSpec_result, "-.", color="green", linewidth=3, zorder=3)
# ax[:plot](t, stanSLS_result, "o", alpha=0.1, color="blue", linewidth=3, zorder=1)

# ax[:tick_params]("both", labelsize=24)
# ax[:grid]("on", alpha=0.4)
# show()

# fraczen = FractionalZener([2.0, 1.0, 1.0, 0.0, 0.5, 0.0])
# fracSLS = FractionalSLS([2.0, 1.0, 1.0, 0.5])
# fracSpec = FractionalSpecial([2.0, 1.0, 0.0, 0.5])
# stanSLS = SLS([2.0, 1.0, 0.5])

# # # J
# fig, ax = subplots()

# fraczen_result = fraczen.J(t, fraczen.parameters)
# fracSLS_result = fracSLS.J(t, fracSLS.parameters)
# fracSpec_result = fracSpec.J(t, fracSpec.parameters)
# stanSLS_result = stanSLS.J(t, stanSLS.parameters)

# ax[:plot](t, fraczen_result, "-", color="pink", linewidth=3, zorder=2)
# ax[:plot](t, fracSLS_result, "--", color="red", linewidth=3, zorder=3)
# ax[:plot](t, fracSpec_result, "-.", color="green", linewidth=3, zorder=3)
# ax[:plot](t, stanSLS_result, "o", alpha=0.1, color="blue", linewidth=3, zorder=1)

# ax[:tick_params]("both", labelsize=24)
# ax[:grid]("on", alpha=0.4)
# show()