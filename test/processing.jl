function _resample_strainonly()
    t0 = collect(0.0:0.01:1.0)
    ϵ0 = t0.^2

    t1 = collect(0.0:0.1:1.0)
    ϵ1 = t1.^2

    data0 = RheoTimeData(t = t0, ϵ = ϵ0)
    
    dataout = resample(data0, -10)

    dataout.t==t1 && dataout.ϵ==ϵ1 
end
@test _resample_strainonly()

function _resample_stressonly()
    t0 = collect(0.0:0.01:1.0)
    σ0 = t0.^2

    t1 = collect(0.0:0.1:1.0)
    σ1 = t1.^2

    data0 = RheoTimeData(t = t0, σ = σ0)
    
    dataout = resample(data0, -10)

    dataout.t==t1 && dataout.σ==σ1 
end
@test _resample_stressonly()

function _resample_stressandstrain()
    t0 = collect(0.0:0.01:1.0)
    ϵ0 = t0.^2
    σ0 = t0.^3

    t1 = collect(0.0:0.1:1.0)
    ϵ1 = t1.^2
    σ1 = t1.^3

    data0 = RheoTimeData(t = t0, σ = σ0, ϵ = ϵ0)
    
    dataout = resample(data0, -10)

    dataout.t==t1 && dataout.ϵ==ϵ1 && dataout.σ==σ1 
end
@test _resample_stressandstrain()

function _resample_stressandstrain_multiplesections()
    t0 = collect(-10.0:0.5:10.0)
    ϵ0 = 2*t0
    σ0 = 2*t0

    # data for comparison composed
    # of multiple sections
    t1a = collect(-10.0:0.1:-5.1)
    t1b = collect(-5.0:1.0:0.0)
    t1c = collect(0.05:0.05:4.95)
    t1d = collect(5.0:2.0:10.0)
    t1e = [t0[end]]

    t1 = vcat(t1a, t1b, t1c, t1d, t1e)
    ϵ1 = 2*t1
    σ1 = 2*t1

    data0 = RheoTimeData(t = t0, σ = σ0, ϵ = ϵ0)
    
    dataout = resample(data0, [5, -2, 10, -4]; time_boundaries = [-10.0, -5.0, 0.0, 5.0, 10.0])

    dataout.t==t1 && dataout.ϵ==ϵ1 && dataout.σ==σ1 
end
@test _resample_stressandstrain_multiplesections()

function _cutting_strainonly()
    t0 = collect(0.0:1.0:5.0)
    ϵ0 = t0.^2

    t1 = collect(1.0:1.0:5.0)
    ϵ1 = t1.^2

    data0 = RheoTimeData(t = t0, ϵ = ϵ0)
    
    dataout = cutting(data0, 1.0, 5.0)

    dataout.t==t1 && dataout.ϵ==ϵ1 
end
@test _cutting_strainonly()

function _cutting_stressonly()
    t0 = collect(0.0:1.0:5.0)
    σ0 = t0.^2

    t1 = collect(1.0:1.0:5.0)
    σ1 = t1.^2

    data0 = RheoTimeData(t = t0, σ = σ0)
    
    dataout = cutting(data0, 1.0, 5.0)

    dataout.t==t1 && dataout.σ==σ1 
end
@test _cutting_stressonly()

function _cutting_stressandstrain()
    t0 = collect(0.0:1.0:5.0)
    ϵ0 = t0.^2
    σ0 = t0.^3

    t1 = collect(1.0:1.0:5.0)
    ϵ1 = t1.^2
    σ1 = t1.^3

    data0 = RheoTimeData(t = t0, σ = σ0, ϵ = ϵ0)
    
    dataout = cutting(data0, 1.0, 5.0)

    dataout.t==t1 && dataout.ϵ==ϵ1 && dataout.σ==σ1 
end
@test _cutting_stressandstrain()

function _smooth_stressonly()
    f = 1.0 # Hz
    ω = 2*π*f # rad/s
    t0 = collect(0.0:0.01:4/f)
    σ0 = sin.(ω*t0)

    data0 = RheoTimeData(t = t0, σ = σ0)

    dataout = smooth(data0, 1/f; pad="circular")

    all(i -> isapprox(0.5*data0.σ[i], dataout.σ[i], atol=0.1), eachindex(data0.σ))
end
@test _smooth_stressonly()

function _smooth_strainonly()
    f = 1.0 # Hz
    ω = 2*π*f # rad/s
    t0 = collect(0.0:0.01:4/f)
    ϵ0 = sin.(ω*t0)

    data0 = RheoTimeData(t = t0, ϵ = ϵ0)

    dataout = smooth(data0, 1/f; pad="circular")

    all(i -> isapprox(0.5*data0.ϵ[i], dataout.ϵ[i], atol=0.1), eachindex(data0.ϵ))
end
@test _smooth_strainonly()

function _smooth_stressandstrain()
    f = 1.0 # Hz
    ω = 2*π*f # rad/s
    t0 = collect(0.0:0.01:4/f)
    ϵ0 = sin.(ω*t0)
    σ0 = sin.(ω*t0)

    data0 = RheoTimeData(t = t0, ϵ = ϵ0, σ = σ0)

    dataout = smooth(data0, 1/f; pad="circular")

    test1 = all(i -> isapprox(0.5*data0.ϵ[i], dataout.ϵ[i], atol=0.1), eachindex(data0.ϵ))
    test2 = all(i -> isapprox(0.5*data0.σ[i], dataout.σ[i], atol=0.1), eachindex(data0.σ))

    test1 && test2
end
@test _smooth_stressandstrain()

function _extract_timefromtimedata()
    f = 1.0 # Hz
    ω = 2*π*f # rad/s
    t0 = collect(0.0:0.01:4/f)
    ϵ0 = sin.(ω*t0)
    σ0 = sin.(ω*t0)

    data0 = RheoTimeData(t = t0, ϵ = ϵ0, σ = σ0)

    dataout = extract(data0, time_only)

    data0.t==dataout.t && dataout.σ==[] && dataout.ϵ==[]
end
@test _extract_timefromtimedata()

function _extract_timestressfromtimedata()
    f = 1.0 # Hz
    ω = 2*π*f # rad/s
    t0 = collect(0.0:0.01:4/f)
    ϵ0 = sin.(ω*t0)
    σ0 = sin.(ω*t0)

    data0 = RheoTimeData(t = t0, ϵ = ϵ0, σ = σ0)

    dataout = extract(data0, stress_only)

    data0.t==dataout.t && data0.σ==dataout.σ && dataout.ϵ==[]
end
@test _extract_timestressfromtimedata()

function _extract_timestrainfromtimedata()
    f = 1.0 # Hz
    ω = 2*π*f # rad/s
    t0 = collect(0.0:0.01:4/f)
    ϵ0 = sin.(ω*t0)
    σ0 = sin.(ω*t0)

    data0 = RheoTimeData(t = t0, ϵ = ϵ0, σ = σ0)

    dataout = extract(data0, strain_only)

    data0.t==dataout.t && dataout.σ==[] && data0.ϵ==dataout.ϵ
end
@test _extract_timestrainfromtimedata()

function _extract_freqfromfreqdata()
    ω0 = collect(0.0:0.01:10.0)
    Gp0 = 2*ω0
    Gpp0 = 3*ω0

    data0 = RheoFreqData(ω = ω0, Gp = Gp0, Gpp = Gpp0)

    dataout = extract(data0, freq_only)

    data0.ω==dataout.ω && dataout.Gp==[] && dataout.Gpp==[]
end
@test _extract_freqfromfreqdata()

function _modelfit_const_ramp(tol)
    dt = 0.01
    t = Vector{RHEOS.RheoFloat}(0.0:dt:20.0)
    exact_response = 1.0 .- exp.(-t)
    ramp_loading = t
    ramp_loading_derivative = RHEOS.derivBD(ramp_loading, t)

    modulus = quote α*exp(-t/β) end
    model = RheoModelClass(name = "testmodel", p = [:α, :β], J = modulus, info="none")

    data0 = RheoTimeData(t = t, ϵ = exact_response, σ = ramp_loading)

    init_params = (α=1.0, β=1.0)
    modelout = modelfit(data0, model, stress_imposed, p0=init_params, lo=(α=0.9, β=0.9), hi=(α=1.1, β=1.1))

    found_params = modelout.params
    isapprox(collect(values(found_params)), collect(values(init_params)), atol = tol)
end
@test _modelfit_const_ramp(tol)

function _modelfit_const_ramp_nobounds(tol)
    dt = 0.01
    t = Vector{RHEOS.RheoFloat}(0.0:dt:20.0)
    exact_response = 1.0 .- exp.(-t)
    ramp_loading = t
    ramp_loading_derivative = RHEOS.derivBD(ramp_loading, t)

    modulus = quote α*exp(-t/β) end
    model = RheoModelClass(name = "testmodel", p = [:α, :β], J = modulus, info="none")

    data0 = RheoTimeData(t = t, ϵ = exact_response, σ = ramp_loading)

    init_params = (α=1.0, β=1.0)
    modelout = modelfit(data0, model, stress_imposed, p0=init_params)

    found_params = modelout.params
    isapprox(collect(values(found_params)), collect(values(init_params)), atol = tol)
end
@test _modelfit_const_ramp_nobounds(tol)

function _modelfit_const_ramp_nobounds_singleparam(tol)
    dt = 0.01
    t = Vector{RHEOS.RheoFloat}(0.0:dt:20.0)
    exact_response = 1.0 .- exp.(-t)
    ramp_loading = t
    ramp_loading_derivative = RHEOS.derivBD(ramp_loading, t)

    modulus = quote α*exp(-t) end
    model = RheoModelClass(name = "testmodel", p = [:α], G = modulus, info="none")

    data0 = RheoTimeData(t = t, ϵ = ramp_loading, σ = exact_response)

    init_params = (α=1.0,)
    modelout = modelfit(data0, model, strain_imposed, p0=init_params)

    found_params = modelout.params
    isapprox(collect(values(found_params)), collect(values(init_params)), atol = tol)
end
@test _modelfit_const_ramp_nobounds_singleparam(tol)

function _modelfit_const_ramp_nobounds(tol)
    dt = 0.01
    t = Vector{RHEOS.RheoFloat}(0.0:dt:20.0)
    exact_response = 1.0 .- exp.(-t)
    ramp_loading = t
    ramp_loading_derivative = RHEOS.derivBD(ramp_loading, t)

    modulus = quote α*exp(-t/β) end
    model = RheoModelClass(name = "testmodel", p = [:α, :β], J = modulus, info="none")

    data0 = RheoTimeData(t = t, ϵ = exact_response, σ = ramp_loading)

    init_params = (α=1.0, β=1.0)
    modelout = modelfit(data0, model, stress_imposed, p0=init_params)

    found_params = modelout.params
    isapprox(collect(values(found_params)), collect(values(init_params)), atol = tol)
end
@test _modelfit_const_ramp_nobounds(tol)

function _modelfit_var_ramp_nonsing(tol)
    dt = 0.01
    t = [Vector{RHEOS.RheoFloat}(0.0:dt:(15.0-dt)); Vector{RHEOS.RheoFloat}(15.0:10*dt:20.0)]
    exact_response = 1.0 .- exp.(-t)
    ramp_loading = t

    modulus = quote α*exp.(-t/β) end
    model = RheoModelClass(name = "testmodel", p = [:α, :β], J = modulus, info="none")

    data0 = RheoTimeData(t = t, ϵ = exact_response, σ = ramp_loading)

    init_params = (α=1.3, β=0.7)
    modelout = modelfit(data0, model, stress_imposed, p0=init_params, lo=(α=0.5, β=0.5), hi=(α=1.5, β=1.5))

    found_params = modelout.params
    actual_params = (α=1.0, β=1.0)

    isapprox(collect(values(found_params)), collect(values(actual_params)), atol = tol)
end
@test _modelfit_var_ramp_nonsing(tol)

function _modelfit_const_ramp_sing(tol)
    dt = 0.01
    t = Vector{RHEOS.RheoFloat}(0.0:dt:20.0)
    actual_params = (α=1.0, β=0.5)
    exact_response = actual_params.α*t.^(1.0 - actual_params.β) / (1.0 - actual_params.β)
    ramp_loading = t

    modulus = quote α*t.^(-β) end
    model = RheoModelClass(name = "testmodel", p = [:α, :β], J = modulus, info="none")

    data0 = RheoTimeData(t = t, ϵ = exact_response, σ = ramp_loading)

    init_params = (α=1.3, β=0.7)
    modelout = modelfit(data0, model, stress_imposed, p0=init_params, lo=(α=0.5, β=0.2), hi=(α=1.5, β=1.5))

    found_params = modelout.params

    isapprox(collect(values(found_params)), collect(values(actual_params)), atol = 5*tol)
end
@test _modelfit_const_ramp_sing(tol)

function _modelfit_var_ramp_sing(tol)
    dt = 0.01
    t = [Vector{RHEOS.RheoFloat}(0.0:dt:(15.0-dt)); Vector{RHEOS.RheoFloat}(15.0:10*dt:20.0)]
    actual_params = (α=1.0, β=0.5)
    exact_response = actual_params.α*t.^(1.0 - actual_params.β) / (1.0 - actual_params.β)
    ramp_loading = t

    modulus = quote α*t.^(-β) end
    model = RheoModelClass(name = "testmodel", p = [:α, :β], J = modulus, info="none")

    data0 = RheoTimeData(t = t, ϵ = exact_response, σ = ramp_loading)

    init_params = (α=1.3, β=0.7)
    modelout = modelfit(data0, model, stress_imposed, p0=init_params, lo=(α=0.5, β=0.2), hi=(α=1.5, β=1.5))

    found_params = modelout.params

    isapprox(collect(values(found_params)), collect(values(actual_params)), atol = 5*tol)
end
@test _modelfit_var_ramp_sing(tol)

function _modelpredict_nonsing_const(tol)
    dt = 0.01
    t = Vector{RHEOS.RheoFloat}(0.0:dt:20.0)
    exact_response = 1 .- exp.(-t)
    loading = t
    
    modulus = quote α*exp.(-t/β) end
    modelclass = RheoModelClass(name = "testmodel", p = [:α, :β], G = modulus, info="none")
    model = RheoModel(modelclass, α=1.0, β=1.0)
    data0 = RheoTimeData(t = t, ϵ = loading)

    computed_response = modelpredict(data0, model)

    all(i -> isapprox(exact_response[i], computed_response.σ[i], atol=tol), eachindex(exact_response))
end
@test _modelpredict_nonsing_const(tol)

function _modelpredict_nonsing_var(tol)
    dt = 0.01
    t = [Vector{RHEOS.RheoFloat}(0.0:dt:(15.0-dt)); Vector{RHEOS.RheoFloat}(15.0:10*dt:20.0)]
    exact_response = 1 .- exp.(-t)
    loading = t
    
    modulus = quote α*exp.(-t/β) end
    modelclass = RheoModelClass(name = "testmodel", p = [:α, :β], G = modulus, info="none")
    model = RheoModel(modelclass, α=1.0, β=1.0)
    data0 = RheoTimeData(t = t, ϵ = loading)

    computed_response = modelpredict(data0, model)

    all(i -> isapprox(exact_response[i], computed_response.σ[i], atol=tol), eachindex(exact_response))
end
@test _modelpredict_nonsing_var(tol)

function _modelpredict_sing_const(tol)
    dt = 0.001
    t = Vector{RHEOS.RheoFloat}(0.0:dt:20.0)
    exact_response = t.^(1.0 - 0.5) / (1.0 - 0.5)
    loading = t

    modulus = quote α*t.^(-β) end
    modelclass = RheoModelClass(name = "testmodel", p = [:α, :β], G = modulus, info="none")
    model = RheoModel(modelclass, α=1.0, β=0.5)
    data0 = RheoTimeData(t = t, ϵ = loading)

    computed_response = modelpredict(data0, model)
    
    all(i -> isapprox(exact_response[i], computed_response.σ[i], atol=5*tol), eachindex(exact_response))
end
@test _modelpredict_sing_const(tol)

function _modelpredict_sing_var(tol)
    dt = 0.01
    t = [Vector{RHEOS.RheoFloat}(0.0:dt:(15.0-dt)); Vector{RHEOS.RheoFloat}(15.0:10*dt:20.0)]
    exact_response = t.^(1.0 - 0.5) / (1.0 - 0.5)
    loading = t

    modulus = quote α*t.^(-β) end
    modelclass = RheoModelClass(name = "testmodel", p = [:α, :β], G = modulus, info="none")
    model = RheoModel(modelclass, α=1.0, β=0.5)
    data0 = RheoTimeData(t = t, ϵ = loading)
    
    computed_response = modelpredict(data0, model)
    plot(t, computed_response.σ)
    plot(t, exact_response)
    println(maximum(abs.(computed_response.σ - exact_response)))
    all(i -> isapprox(exact_response[i], computed_response.σ[i], atol=5*tol), eachindex(exact_response))
end
@test _modelpredict_sing_var(tol)
