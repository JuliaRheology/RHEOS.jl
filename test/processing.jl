println("===============================================")
println("Testing processing.jl")
println("===============================================")


function _resample_timeonly()
    d = timeline(t_start=0, t_end=1, step=1//4)
    d = resample(d,dt=0.1)
    (length(d.t)==11) && (d.t[2]==0.1)
end
@test _resample_timeonly()

function _resample_stressonly()
    d = RheoTimeData(σ=[1,2,4,7], t=[0,1,3,6])
    d = resample(d) # make timestep uniform without changing number of points
    (length(d.t)==4) && (d.t[2]==2) && (d.σ[2]≈3)

end
@test _resample_stressonly()

function _resample_strainonly()
    d = RheoTimeData(ϵ=[1,2,4,7], t=[0,1,3,6])
    d = resample(d,t=1:1//4:2)
    (length(d.t)==5) && (d.t[2]==5//4) && (d.ϵ[2]≈9//4)
end
@test _resample_strainonly()

function _resample_stressstrain()
    d = RheoTimeData(ϵ=[1,2,4,7], t=[0,1,3,6], σ=[7,4,2,1])
    d = resample(d,scale=2)
    (length(d.t)==7) && (d.t[4]≈1.875) && (d.σ[4]≈2.6328125) && (d.ϵ[4]≈2.875)
end
@test _resample_stressstrain()



# function _resample_strainonly()
#     t0 = collect(0.0:0.01:1.0)
#     ϵ0 = t0.^2
#
#     t1 = collect(0.0:0.1:1.0)
#     ϵ1 = t1.^2
#
#     data0 = RheoTimeData(t = t0, ϵ = ϵ0)
#
#     dataout = resample(data0, -10)
#
#     dataout.t==t1 && dataout.ϵ==ϵ1
# end
# @test _resample_strainonly()
#
# function _resample_stressonly()
#     t0 = collect(0.0:0.01:1.0)
#     σ0 = t0.^2
#
#     t1 = collect(0.0:0.1:1.0)
#     σ1 = t1.^2
#
#     data0 = RheoTimeData(t = t0, σ = σ0)
#
#     dataout = resample(data0, -10)
#
#     dataout.t==t1 && dataout.σ==σ1
# end
# @test _resample_stressonly()
#
# function _resample_stressandstrain()
#     t0 = collect(0.0:0.01:1.0)
#     ϵ0 = t0.^2
#     σ0 = t0.^3
#
#     t1 = collect(0.0:0.1:1.0)
#     ϵ1 = t1.^2
#     σ1 = t1.^3
#
#     data0 = RheoTimeData(t = t0, σ = σ0, ϵ = ϵ0)
#
#     dataout = resample(data0, -10)
#
#     dataout.t==t1 && dataout.ϵ==ϵ1 && dataout.σ==σ1
# end
# @test _resample_stressandstrain()
#
# function _resample_stressandstrain_multiplesections()
#     t0 = collect(-10.0:0.5:10.0)
#     ϵ0 = 2*t0
#     σ0 = 2*t0
#
#     # data for comparison composed
#     # of multiple sections
#     t1a = collect(-10.0:0.1:-5.1)
#     t1b = collect(-5.0:1.0:0.0)
#     t1c = collect(0.05:0.05:4.95)
#     t1d = collect(5.0:2.0:10.0)
#     t1e = [t0[end]]
#
#     t1 = vcat(t1a, t1b, t1c, t1d, t1e)
#     ϵ1 = 2*t1
#     σ1 = 2*t1
#
#     data0 = RheoTimeData(t = t0, σ = σ0, ϵ = ϵ0)
#
#     dataout = resample(data0, [5, -2, 10, -4]; time_boundaries = [-10.0, -5.0, 0.0, 5.0, 10.0])
#
#     dataout.t==t1 && dataout.ϵ==ϵ1 && dataout.σ==σ1
# end
# @test _resample_stressandstrain_multiplesections()

function _indexweight_tworegiondownsample_includelastel()
    timedata = timeline(t_start=0.0, t_end=10.0, step=1.0)

    indices = indexweight(timedata; elperiods = [-2, -3], time_boundaries=[0.0, 5.0, 10.0])

    indices == [1, 3, 5, 6, 9, 11]
end
@test _indexweight_tworegiondownsample_includelastel()

function _indexweight_tworegiondownsample_nolastel()
    timedata = timeline(t_start=0.0, t_end=10.0, step=1.0)

    indices = indexweight(timedata; elperiods = [-2, -3], time_boundaries=[0.0, 5.0, 10.0], includelast=false)

    indices == [1, 3, 5, 6, 9]
end
@test _indexweight_tworegiondownsample_nolastel()

function _indexweight_tworegionupsample_includelastel()
    timedata = timeline(t_start=0.0, t_end=5.0, step=1.0)

    indices = indexweight(timedata; elperiods = [2, 3], time_boundaries=[0.0, 3.0, 4.0])

    indices == [1, 1, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5]
end
@test _indexweight_tworegionupsample_includelastel()

function _indexweight_tworegionupsample_nolastel()
    timedata = timeline(t_start=0.0, t_end=5.0, step=1.0)

    indices = indexweight(timedata; elperiods = [2, 3], time_boundaries=[0.0, 3.0, 4.0], includelast=false)

    indices == [1, 1, 2, 2, 3, 3, 4, 4, 4]
end
@test _indexweight_tworegionupsample_nolastel()

function _indexweight_mixed_includelastel_finishup()
    timedata = timeline(t_start=0.0, t_end=10.0, step=1.0)

    indices = indexweight(timedata; elperiods = [-2, 3, -3, 4], time_boundaries=[0.0, 3.0, 5.0, 9.0, 10.0])

    indices == [1, 3, 4, 4, 4, 5, 5, 5, 6, 9, 10, 10, 10, 10, 11, 11, 11, 11]
end
@test _indexweight_mixed_includelastel_finishup()

function _indexweight_mixed_includelastel_finishdown()
    timedata = timeline(t_start=0.0, t_end=10.0, step=1.0)

    indices = indexweight(timedata; elperiods = [2, -3, 2, -2], time_boundaries=[0.0, 3.0, 5.0, 7.0, 10.0])

    indices == [1, 1, 2, 2, 3, 3, 4, 6, 6, 7, 7, 8, 10, 11]
end
@test _indexweight_mixed_includelastel_finishdown()

function _indexweight_mixed_nolastel_finishup()
    timedata = timeline(t_start=0.0, t_end=10.0, step=1.0)

    indices = indexweight(timedata; elperiods = [-2, 3, -3, 4], time_boundaries=[0.0, 3.0, 5.0, 9.0, 10.0], includelast=false)

    indices == [1, 3, 4, 4, 4, 5, 5, 5, 6, 9, 10, 10, 10, 10]
end
@test _indexweight_mixed_nolastel_finishup()

function _indexweight_mixed_nolastel_finishdown()
    timedata = timeline(t_start=0.0, t_end=10.0, step=1.0)

    indices = indexweight(timedata; elperiods = [2, -3, 2, -2], time_boundaries=[0.0, 3.0, 5.0, 7.0, 10.0], includelast=false)

    indices == [1, 1, 2, 2, 3, 3, 4, 6, 6, 7, 7, 8, 10]
end
@test _indexweight_mixed_nolastel_finishdown()

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

function _fill_init_params_nothing()
    RHEOS.fill_init_params(Fract_Zener, nothing)==[0.5, 0.8, 0.5, 0.5, 0.5, 0.5]
end
@test _fill_init_params_nothing()

function _fill_init_params_fullyspecified()
    RHEOS.fill_init_params(Fract_Zener, (cₐ=1.0, a=0.51, cᵦ=1.5, β=0.4, cᵧ=2.0, γ=0.7))==[1.0, 0.51, 1.5, 0.4, 2.0, 0.7]
end
@test _fill_init_params_fullyspecified()

function _fill_init_params_semispecified1()
    RHEOS.fill_init_params(Fract_Zener, (a=0.51, cᵦ=1.5, β=0.4, cᵧ=2.0, γ=0.7))==[0.5, 0.51, 1.5, 0.4, 2.0, 0.7]
end
@test _fill_init_params_semispecified1()

function _fill_init_params_semispecified2()
    RHEOS.fill_init_params(Fract_Zener, (cᵦ=1.5, β=0.51, cᵧ=2.0, γ=0.7))==[0.5, 0.999, 1.5, 0.51, 2.0, 0.7]
end
@test _fill_init_params_semispecified2()

function _fill_init_params_semispecified3()
    RHEOS.fill_init_params(Fract_Zener, (cᵦ=1.5, β=0.4, cᵧ=2.0, γ=0.7))==[0.5, 0.8, 1.5, 0.4, 2.0, 0.7]
end
@test _fill_init_params_semispecified3()

function _fill_init_params_semispecified4()
    RHEOS.fill_init_params(Fract_Zener, (cₐ=1.0, a=0.4, cᵦ=1.5, cᵧ=2.0, γ=0.7))==[1.0, 0.4, 1.5, 0.2, 2.0, 0.7]
end
@test _fill_init_params_semispecified4()

function _fill_init_params_semispecified5()
    RHEOS.fill_init_params(Fract_Zener, (cₐ=1.0, cᵦ=1.5, cᵧ=2.0, γ=0.7))==[1.0, 0.8, 1.5, 0.5, 2.0, 0.7]
end
@test _fill_init_params_semispecified5()

function _fill_lower_bounds_nothing()
    isnothing(RHEOS.fill_lower_bounds(Fract_Zener, nothing))
end
@test _fill_lower_bounds_nothing()

function _fill_lower_bounds_fullyspecified()
    RHEOS.fill_lower_bounds(Fract_Zener, (cₐ=1.0, a=0.4, cᵦ=1.5, β=0.51, cᵧ=2.0, γ=0.7))==[1.0, 0.4, 1.5, 0.51, 2.0, 0.7]
end
@test _fill_lower_bounds_fullyspecified()

function _fill_lower_bounds_semispecified()
    RHEOS.fill_lower_bounds(Fract_Zener, (a=0.4, β=0.51, cᵧ=2.0, γ=0.7))==[0.0, 0.4, 0.0, 0.51, 2.0, 0.7]
end
@test _fill_lower_bounds_semispecified()

function _fill_upper_bounds_nothing()
    isnothing(RHEOS.fill_upper_bounds(Fract_Zener, nothing))
end
@test _fill_upper_bounds_nothing()

function _fill_upper_bounds_fullyspecified()
    RHEOS.fill_upper_bounds(Fract_Zener, (cₐ=1.0, a=0.4, cᵦ=1.5, β=0.51, cᵧ=2.0, γ=0.7))==[1.0, 0.4, 1.5, 0.51, 2.0, 0.7]
end
@test _fill_upper_bounds_fullyspecified()

function _fill_upper_bounds_semispecified()
    RHEOS.fill_upper_bounds(Fract_Zener, (a=0.4, β=0.51, cᵧ=2.0, γ=0.7))==[Inf, 0.4, Inf, 0.51, 2.0, 0.7]
end
@test _fill_upper_bounds_semispecified()

function _modelfit_const_ramp_relax(tol)
    dt = 0.01
    t = Vector{RheoFloat}(0.0:dt:20.0)
    exact_response = 1.0 .- exp.(-t)
    ramp_loading = t

    modulus = quote α*exp(-t/β) end
    model = RheoModelClass(name = "testmodel", p = [:α, :β], G = modulus, info="none")

    data0 = RheoTimeData(t = t, ϵ = ramp_loading, σ = exact_response)

    init_params = (α=1.0, β=1.0)
    modelout = modelfit(data0, model, strain_imposed, p0=init_params, lo=(α=0.9, β=0.9), hi=(α=1.1, β=1.1))

    found_params = modelout.params
    isapprox(collect(values(found_params)), collect(values(init_params)), atol = tol)
end
@test _modelfit_const_ramp_relax(tol)

function _modelfit_const_ramp_creep(tol)
    dt = 0.01
    t = Vector{RheoFloat}(0.0:dt:20.0)
    exact_response = 1.0 .- exp.(-t)
    ramp_loading = t

    modulus = quote α*exp(-t/β) end
    model = RheoModelClass(name = "testmodel", p = [:α, :β], J = modulus, info="none")

    data0 = RheoTimeData(t = t, ϵ = exact_response, σ = ramp_loading)

    init_params = (α=1.0, β=1.0)
    modelout = modelfit(data0, model, stress_imposed, p0=init_params, lo=(α=0.9, β=0.9), hi=(α=1.1, β=1.1))

    found_params = modelout.params
    isapprox(collect(values(found_params)), collect(values(init_params)), atol = tol)
end
@test _modelfit_const_ramp_creep(tol)

function _modelfit_const_ramp_nobounds_relax(tol)
    dt = 0.01
    t = Vector{RheoFloat}(0.0:dt:20.0)
    exact_response = 1.0 .- exp.(-t)
    ramp_loading = t

    modulus = quote α*exp(-t/β) end
    model = RheoModelClass(name = "testmodel", p = [:α, :β], G = modulus, info="none")

    data0 = RheoTimeData(t = t, ϵ = ramp_loading, σ = exact_response)

    init_params = (α=1.0, β=1.0)
    modelout = modelfit(data0, model, strain_imposed, p0=init_params)

    found_params = modelout.params
    isapprox(collect(values(found_params)), collect(values(init_params)), atol = tol)
end
@test _modelfit_const_ramp_nobounds_relax(tol)

function _modelfit_const_ramp_nobounds_creep(tol)
    dt = 0.01
    t = Vector{RheoFloat}(0.0:dt:20.0)
    exact_response = 1.0 .- exp.(-t)
    ramp_loading = t

    modulus = quote α*exp(-t/β) end
    model = RheoModelClass(name = "testmodel", p = [:α, :β], J = modulus, info="none")

    data0 = RheoTimeData(t = t, ϵ = exact_response, σ = ramp_loading)

    init_params = (α=1.0, β=1.0)
    modelout = modelfit(data0, model, stress_imposed, p0=init_params)

    found_params = modelout.params
    isapprox(collect(values(found_params)), collect(values(init_params)), atol = tol)
end
@test _modelfit_const_ramp_nobounds_creep(tol)

function _modelfit_const_ramp_nobounds_singleparam_relax(tol)
    dt = 0.01
    t = Vector{RheoFloat}(0.0:dt:20.0)
    exact_response = 1.0 .- exp.(-t)
    ramp_loading = t

    modulus = quote α*exp(-t) end
    model = RheoModelClass(name = "testmodel", p = [:α], G = modulus, info="none")

    data0 = RheoTimeData(t = t, ϵ = ramp_loading, σ = exact_response)

    init_params = (α=1.0,)
    modelout = modelfit(data0, model, strain_imposed, p0=init_params)

    found_params = modelout.params
    isapprox(collect(values(found_params)), collect(values(init_params)), atol = tol)
end
@test _modelfit_const_ramp_nobounds_singleparam_relax(tol)

function _modelfit_const_ramp_nobounds_singleparam_creep(tol)
    dt = 0.01
    t = Vector{RheoFloat}(0.0:dt:20.0)
    exact_response = 1.0 .- exp.(-t)
    ramp_loading = t

    modulus = quote α*exp(-t) end
    model = RheoModelClass(name = "testmodel", p = [:α], J = modulus, info="none")

    data0 = RheoTimeData(t = t, ϵ = exact_response, σ = ramp_loading)

    init_params = (α=1.0,)
    modelout = modelfit(data0, model, stress_imposed, p0=init_params)

    found_params = modelout.params
    isapprox(collect(values(found_params)), collect(values(init_params)), atol = tol)
end
@test _modelfit_const_ramp_nobounds_singleparam_creep(tol)

function _modelfit_const_ramp_nobounds_relax(tol)
    dt = 0.01
    t = Vector{RheoFloat}(0.0:dt:20.0)
    exact_response = 1.0 .- exp.(-t)
    ramp_loading = t

    modulus = quote α*exp(-t/β) end
    model = RheoModelClass(name = "testmodel", p = [:α, :β], G = modulus, info="none")

    data0 = RheoTimeData(t = t, ϵ = ramp_loading, σ = exact_response)

    init_params = (α=1.0, β=1.0)
    modelout = modelfit(data0, model, strain_imposed, p0=init_params)

    found_params = modelout.params
    isapprox(collect(values(found_params)), collect(values(init_params)), atol = tol)
end
@test _modelfit_const_ramp_nobounds_relax(tol)

function _modelfit_const_ramp_nobounds_creep(tol)
    dt = 0.01
    t = Vector{RheoFloat}(0.0:dt:20.0)
    exact_response = 1.0 .- exp.(-t)
    ramp_loading = t

    modulus = quote α*exp(-t/β) end
    model = RheoModelClass(name = "testmodel", p = [:α, :β], J = modulus, info="none")

    data0 = RheoTimeData(t = t, ϵ = exact_response, σ = ramp_loading)

    init_params = (α=1.0, β=1.0)
    modelout = modelfit(data0, model, stress_imposed, p0=init_params)

    found_params = modelout.params
    isapprox(collect(values(found_params)), collect(values(init_params)), atol = tol)
end
@test _modelfit_const_ramp_nobounds_creep(tol)

function _modelfit_var_ramp_nonsing_relax(tol)
    dt = 0.01
    t = [Vector{RheoFloat}(0.0:dt:(15.0-dt)); Vector{RheoFloat}(15.0:10*dt:20.0)]
    exact_response = 1.0 .- exp.(-t)
    ramp_loading = t

    modulus = quote α*exp.(-t/β) end
    model = RheoModelClass(name = "testmodel", p = [:α, :β], G = modulus, info="none")

    data0 = RheoTimeData(t = t, ϵ = ramp_loading, σ = exact_response)

    init_params = (α=1.3, β=0.7)
    modelout = modelfit(data0, model, strain_imposed, p0=init_params, lo=(α=0.5, β=0.5), hi=(α=1.5, β=1.5))

    found_params = modelout.params
    actual_params = (α=1.0, β=1.0)

    isapprox(collect(values(found_params)), collect(values(actual_params)), atol = tol)
end
@test _modelfit_var_ramp_nonsing_relax(tol)

function _modelfit_var_ramp_nonsing_creep(tol)
    dt = 0.01
    t = [Vector{RheoFloat}(0.0:dt:(15.0-dt)); Vector{RheoFloat}(15.0:10*dt:20.0)]
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
@test _modelfit_var_ramp_nonsing_creep(tol)

function _modelfit_const_ramp_sing_relax(tol)
    dt = 0.01
    t = Vector{RheoFloat}(0.0:dt:20.0)
    actual_params = (α=1.0, β=0.5)
    exact_response = actual_params.α*t.^(1.0 - actual_params.β) / (1.0 - actual_params.β)
    ramp_loading = t

    modulus = quote α*t.^(-β) end
    model = RheoModelClass(name = "testmodel", p = [:α, :β], G = modulus, info="none")

    data0 = RheoTimeData(t = t, ϵ = ramp_loading, σ = exact_response)

    init_params = (α=1.3, β=0.7)
    modelout = modelfit(data0, model, strain_imposed, p0=init_params, lo=(α=0.5, β=0.2), hi=(α=1.5, β=1.5))

    found_params = modelout.params

    isapprox(collect(values(found_params)), collect(values(actual_params)), atol = 5*tol)
end
@test _modelfit_const_ramp_sing_relax(tol)

function _modelfit_const_ramp_sing_creep(tol)
    dt = 0.01
    t = Vector{RheoFloat}(0.0:dt:20.0)
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
@test _modelfit_const_ramp_sing_creep(tol)

function _modelfit_var_ramp_sing_relax(tol)
    dt = 0.01
    t = [Vector{RheoFloat}(0.0:dt:(15.0-dt)); Vector{RheoFloat}(15.0:10*dt:20.0)]
    actual_params = (α=1.0, β=0.5)
    exact_response = actual_params.α*t.^(1.0 - actual_params.β) / (1.0 - actual_params.β)
    ramp_loading = t

    modulus = quote α*t.^(-β) end
    model = RheoModelClass(name = "testmodel", p = [:α, :β], G = modulus, info="none")

    data0 = RheoTimeData(t = t, ϵ = ramp_loading, σ = exact_response)

    init_params = (α=1.3, β=0.7)
    modelout = modelfit(data0, model, strain_imposed, p0=init_params, lo=(α=0.5, β=0.2), hi=(α=1.5, β=1.5))

    found_params = modelout.params

    isapprox(collect(values(found_params)), collect(values(actual_params)), atol = 5*tol)
end
@test _modelfit_var_ramp_sing_relax(tol)

function _modelfit_const_ramp_creep_weighted_selftest(tol)
    dt = 0.01
    t = Vector{RheoFloat}(0.0:dt:20.0)
    exact_response = 1.0 .- exp.(-t)
    ramp_loading = t

    modulus = quote α*exp(-t/β) end
    model = RheoModelClass(name = "testmodel", p = [:α, :β], J = modulus, info="none")

    data0 = RheoTimeData(t = t, ϵ = exact_response, σ = ramp_loading)

    actual_params = (α=1.0, β=1.0)
    init_params = (α=1.2, β=0.8)
    modelout = modelfit(data0, model, stress_imposed, p0=init_params, lo=(α=0.7, β=0.7), hi=(α=1.5, β=1.5), weights=collect(Integer, 1:length(t)))

    found_params = modelout.params
    isapprox(collect(values(found_params)), collect(values(actual_params)), atol = tol)
end
@test _modelfit_const_ramp_creep_weighted_selftest(tol)

function _modelfit_const_ramp_sing_creep_weighted_selftest(tol)
    dt = 0.01
    t = Vector{RheoFloat}(0.0:dt:20.0)
    actual_params = (α=1.0, β=0.5)
    exact_response = actual_params.α*t.^(1.0 - actual_params.β) / (1.0 - actual_params.β)
    ramp_loading = t

    modulus = quote α*t.^(-β) end
    model = RheoModelClass(name = "testmodel", p = [:α, :β], J = modulus, info="none")

    data0 = RheoTimeData(t = t, ϵ = exact_response, σ = ramp_loading)

    init_params = (α=1.3, β=0.7)
    modelout = modelfit(data0, model, stress_imposed, p0=init_params, lo=(α=0.5, β=0.2), hi=(α=1.5, β=1.5), weights=collect(Integer, 1:length(t)))

    found_params = modelout.params

    isapprox(collect(values(found_params)), collect(values(actual_params)), atol = 5*tol)
end
@test _modelfit_const_ramp_sing_creep_weighted_selftest(tol)

function _modelfit_const_ramp_creep_weighted_downsampled(tol)
    dt = 0.005
    t = Vector{RheoFloat}(0.0:dt:20.0)
    exact_response = 1.0 .- exp.(-t)
    ramp_loading = t

    modulus = quote α*exp(-t/β) end
    model = RheoModelClass(name = "testmodel", p = [:α, :β], J = modulus, info="none")

    data0 = RheoTimeData(t = t, ϵ = exact_response, σ = ramp_loading)

    indices = indexweight(data0; elperiods = [-2, 2], time_boundaries = [0.0, 5.0, 20.0])

    init_params = (α=1.2, β=0.8)
    modelout_weighted = modelfit(data0, model, stress_imposed, p0=init_params, lo=(α=0.7, β=0.7), hi=(α=1.5, β=1.5), weights=indices)
    modelout = modelfit(data0, model, stress_imposed, p0=init_params, lo=(α=0.7, β=0.7), hi=(α=1.5, β=1.5))

    found_params_weighted = modelout_weighted.params
    found_params = modelout.params
    isapprox(collect(values(found_params_weighted)), collect(values(found_params)), atol = tol)
end
@test _modelfit_const_ramp_creep_weighted_downsampled(tol)

function _modelfit_const_ramp_sing_creep_weighted_downsampled(tol)
    dt = 0.01
    t = Vector{RheoFloat}(0.0:dt:20.0)
    actual_params = (α=1.0, β=0.5)
    exact_response = actual_params.α*t.^(1.0 - actual_params.β) / (1.0 - actual_params.β)
    ramp_loading = t

    modulus = quote α*t.^(-β) end
    model = RheoModelClass(name = "testmodel", p = [:α, :β], J = modulus, info="none")

    data0 = RheoTimeData(t = t, ϵ = exact_response, σ = ramp_loading)

    indices = indexweight(data0; elperiods = [-2, 2], time_boundaries = [0.0, 5.0, 20.0])

    init_params = (α=1.3, β=0.7)
    modelout_weighted = modelfit(data0, model, stress_imposed, p0=init_params, lo=(α=0.5, β=0.2), hi=(α=1.5, β=1.5), weights=indices)
    modelout = modelfit(data0, model, stress_imposed, p0=init_params, lo=(α=0.5, β=0.2), hi=(α=1.5, β=1.5))

    found_params_weighted = modelout_weighted.params
    found_params = modelout.params
    isapprox(collect(values(found_params_weighted)), collect(values(found_params)), atol = tol)
end
@test _modelfit_const_ramp_sing_creep_weighted_downsampled(tol)

function _modelfit_var_ramp_sing_creep(tol)
    dt = 0.01
    t = [Vector{RheoFloat}(0.0:dt:(15.0-dt)); Vector{RheoFloat}(15.0:10*dt:20.0)]
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
@test _modelfit_var_ramp_sing_creep(tol)

function _modelpredict_nonsing_const_relax(tol)
    dt = 0.01
    t = Vector{RheoFloat}(0.0:dt:20.0)
    exact_response = 1 .- exp.(-t)
    loading = t

    modulus = quote α*exp.(-t/β) end
    modelclass = RheoModelClass(name = "testmodel", p = [:α, :β], G = modulus, info="none")
    model = RheoModel(modelclass, α=1.0, β=1.0)
    data0 = RheoTimeData(t = t, ϵ = loading)

    computed_response = modelpredict(data0, model)

    all(i -> isapprox(exact_response[i], computed_response.σ[i], atol=tol), eachindex(exact_response))
end
@test _modelpredict_nonsing_const_relax(tol)

function _modelpredict_nonsing_const_creep(tol)
    dt = 0.01
    t = Vector{RheoFloat}(0.0:dt:20.0)
    exact_response = 1 .- exp.(-t)
    loading = t

    modulus = quote α*exp.(-t/β) end
    modelclass = RheoModelClass(name = "testmodel", p = [:α, :β], J = modulus, info="none")
    model = RheoModel(modelclass, α=1.0, β=1.0)
    data0 = RheoTimeData(t = t, σ = loading)

    computed_response = modelpredict(data0, model)

    all(i -> isapprox(exact_response[i], computed_response.ϵ[i], atol=tol), eachindex(exact_response))
end
@test _modelpredict_nonsing_const_creep(tol)

function _modelpredict_nonsing_var_relax(tol)
    dt = 0.01
    t = [Vector{RheoFloat}(0.0:dt:(15.0-dt)); Vector{RheoFloat}(15.0:10*dt:20.0)]
    exact_response = 1 .- exp.(-t)
    loading = t

    modulus = quote α*exp.(-t/β) end
    modelclass = RheoModelClass(name = "testmodel", p = [:α, :β], G = modulus, info="none")
    model = RheoModel(modelclass, α=1.0, β=1.0)
    data0 = RheoTimeData(t = t, ϵ = loading)

    computed_response = modelpredict(data0, model)

    all(i -> isapprox(exact_response[i], computed_response.σ[i], atol=tol), eachindex(exact_response))
end
@test _modelpredict_nonsing_var_relax(tol)

function _modelpredict_nonsing_var_creep(tol)
    dt = 0.01
    t = [Vector{RheoFloat}(0.0:dt:(15.0-dt)); Vector{RheoFloat}(15.0:10*dt:20.0)]
    exact_response = 1 .- exp.(-t)
    loading = t

    modulus = quote α*exp.(-t/β) end
    modelclass = RheoModelClass(name = "testmodel", p = [:α, :β], J = modulus, info="none")
    model = RheoModel(modelclass, α=1.0, β=1.0)
    data0 = RheoTimeData(t = t, σ = loading)

    computed_response = modelpredict(data0, model)

    all(i -> isapprox(exact_response[i], computed_response.ϵ[i], atol=tol), eachindex(exact_response))
end
@test _modelpredict_nonsing_var_creep(tol)

function _modelpredict_sing_const_relax(tol)
    dt = 0.001
    t = Vector{RheoFloat}(0.0:dt:20.0)
    exact_response = t.^(1.0 - 0.5) / (1.0 - 0.5)
    loading = t

    modulus = quote α*t.^(-β) end
    modelclass = RheoModelClass(name = "testmodel", p = [:α, :β], G = modulus, info="none")
    model = RheoModel(modelclass, α=1.0, β=0.5)
    data0 = RheoTimeData(t = t, ϵ = loading)

    computed_response = modelpredict(data0, model)

    all(i -> isapprox(exact_response[i], computed_response.σ[i], atol=5*tol), eachindex(exact_response))
end
@test _modelpredict_sing_const_relax(tol)

function _modelpredict_sing_const_creep(tol)
    dt = 0.001
    t = Vector{RheoFloat}(0.0:dt:20.0)
    exact_response = t.^(1.0 - 0.5) / (1.0 - 0.5)
    loading = t

    modulus = quote α*t.^(-β) end
    modelclass = RheoModelClass(name = "testmodel", p = [:α, :β], J = modulus, info="none")
    model = RheoModel(modelclass, α=1.0, β=0.5)
    data0 = RheoTimeData(t = t, σ = loading)

    computed_response = modelpredict(data0, model)

    all(i -> isapprox(exact_response[i], computed_response.ϵ[i], atol=5*tol), eachindex(exact_response))
end
@test _modelpredict_sing_const_creep(tol)

function _modelpredict_sing_var_relax(tol)
    dt = 0.01
    t = [Vector{RheoFloat}(0.0:dt:(15.0-dt)); Vector{RheoFloat}(15.0:10*dt:20.0)]
    exact_response = t.^(1.0 - 0.5) / (1.0 - 0.5)
    loading = t

    modulus = quote α*t.^(-β) end
    modelclass = RheoModelClass(name = "testmodel", p = [:α, :β], G = modulus, info="none")
    model = RheoModel(modelclass, α=1.0, β=0.5)
    data0 = RheoTimeData(t = t, ϵ = loading)

    computed_response = modelpredict(data0, model)

    all(i -> isapprox(exact_response[i], computed_response.σ[i], atol=5*tol), eachindex(exact_response))
end
@test _modelpredict_sing_var_relax(tol)

function _modelpredict_sing_var_creep(tol)
    dt = 0.01
    t = [Vector{RheoFloat}(0.0:dt:(15.0-dt)); Vector{RheoFloat}(15.0:10*dt:20.0)]
    exact_response = t.^(1.0 - 0.5) / (1.0 - 0.5)
    loading = t

    modulus = quote α*t.^(-β) end
    modelclass = RheoModelClass(name = "testmodel", p = [:α, :β], J = modulus, info="none")
    model = RheoModel(modelclass, α=1.0, β=0.5)
    data0 = RheoTimeData(t = t, σ = loading)

    computed_response = modelpredict(data0, model)

    all(i -> isapprox(exact_response[i], computed_response.ϵ[i], atol=5*tol), eachindex(exact_response))
end
@test _modelpredict_sing_var_creep(tol)

function _modelstepfit_nonsing_relax(tol)
    dt = 0.01
    t = Vector{RheoFloat}(0.0:dt:20.0)
    actual_params = (α=2.0, β=0.5)
    loading = ones(length(t))
    exact_response = loading[1]*actual_params[1]*exp.(-t/actual_params[2])

    modulus = quote α*exp.(-t/β) end
    model = RheoModelClass(name = "testmodel", p = [:α, :β], G = modulus, info="none")

    data0 = RheoTimeData(t = t, ϵ = loading, σ = exact_response)

    init_params = (α=1.3, β=0.7)
    modelout = modelstepfit(data0, model, strain_imposed, p0=init_params, lo=(α=0.2, β=0.2), hi=(α=3.5, β=1.5))

    found_params = modelout.params

    isapprox(collect(values(actual_params)), collect(values(found_params)), atol=tol)
end
@test _modelstepfit_nonsing_relax(tol)

function _modelstepfit_nonsing_creep(tol)
    dt = 0.01
    t = Vector{RheoFloat}(0.0:dt:20.0)
    actual_params = (α=2.0, β=0.5)
    loading = ones(length(t))
    exact_response = loading[1]*actual_params[1]*exp.(-t/actual_params[2])

    modulus = quote α*exp.(-t/β) end
    model = RheoModelClass(name = "testmodel", p = [:α, :β], J = modulus, info="none")

    data0 = RheoTimeData(t = t, ϵ = exact_response, σ = loading)

    init_params = (α=1.3, β=0.7)
    modelout = modelstepfit(data0, model, stress_imposed, p0=init_params, lo=(α=0.2, β=0.2), hi=(α=3.5, β=1.5))

    found_params = modelout.params

    isapprox(collect(values(actual_params)), collect(values(found_params)), atol=tol)
end
@test _modelstepfit_nonsing_creep(tol)

function _modelstepfit_nonsing_noinitparams_relax(tol)
    dt = 0.01
    t = Vector{RheoFloat}(0.0:dt:20.0)
    actual_params = (α=2.0, β=0.5)
    loading = ones(length(t))
    exact_response = loading[1]*actual_params[1]*exp.(-t/actual_params[2])

    modulus = quote α*exp.(-t/β) end
    model = RheoModelClass(name = "testmodel", p = [:α, :β], G = modulus, info="none")

    data0 = RheoTimeData(t = t, ϵ = loading, σ = exact_response)

    modelout = modelstepfit(data0, model, strain_imposed, lo=(α=0.2, β=0.2), hi=(α=3.5, β=1.5))

    found_params = modelout.params

    isapprox(collect(values(actual_params)), collect(values(found_params)), atol=tol)
end
@test _modelstepfit_nonsing_noinitparams_relax(tol)

function _modelstepfit_nonsing_noinitparams_creep(tol)
    dt = 0.01
    t = Vector{RheoFloat}(0.0:dt:20.0)
    actual_params = (α=2.0, β=0.5)
    loading = ones(length(t))
    exact_response = loading[1]*actual_params[1]*exp.(-t/actual_params[2])

    modulus = quote α*exp.(-t/β) end
    model = RheoModelClass(name = "testmodel", p = [:α, :β], J = modulus, info="none")

    data0 = RheoTimeData(t = t, ϵ = exact_response, σ = loading)

    modelout = modelstepfit(data0, model, stress_imposed, lo=(α=0.2, β=0.2), hi=(α=3.5, β=1.5))

    found_params = modelout.params

    isapprox(collect(values(actual_params)), collect(values(found_params)), atol=tol)
end
@test _modelstepfit_nonsing_noinitparams_creep(tol)

function _modelstepfit_nonsing_nobounds_relax(tol)
    dt = 0.01
    t = Vector{RheoFloat}(0.0:dt:20.0)
    actual_params = (α=2.0, β=0.5)
    loading = ones(length(t))
    exact_response = loading[1]*actual_params[1]*exp.(-t/actual_params[2])

    modulus = quote α*exp.(-t/β) end
    model = RheoModelClass(name = "testmodel", p = [:α, :β], G = modulus, info="none")

    data0 = RheoTimeData(t = t, ϵ = loading, σ = exact_response)

    init_params = (α=1.3, β=0.7)
    modelout = modelstepfit(data0, model, strain_imposed, p0=init_params)

    found_params = modelout.params

    isapprox(collect(values(actual_params)), collect(values(found_params)), atol=tol)
end
@test _modelstepfit_nonsing_nobounds_relax(tol)

function _modelstepfit_nonsing_nobounds_creep(tol)
    dt = 0.01
    t = Vector{RheoFloat}(0.0:dt:20.0)
    actual_params = (α=2.0, β=0.5)
    loading = ones(length(t))
    exact_response = loading[1]*actual_params[1]*exp.(-t/actual_params[2])

    modulus = quote α*exp.(-t/β) end
    model = RheoModelClass(name = "testmodel", p = [:α, :β], J = modulus, info="none")

    data0 = RheoTimeData(t = t, ϵ = exact_response, σ = loading)

    init_params = (α=1.3, β=0.7)
    modelout = modelstepfit(data0, model, stress_imposed, p0=init_params)

    found_params = modelout.params

    isapprox(collect(values(actual_params)), collect(values(found_params)), atol=tol)
end
@test _modelstepfit_nonsing_nobounds_creep(tol)

function _modelstepfit_sing_relax(tol)
    dt = 0.01
    t = Vector{RheoFloat}(0.0:dt:20.0)
    actual_params = (α=2.0, β=0.5)
    loading = ones(length(t))
    exact_response = loading[1]*actual_params[1]*t.^(-actual_params[2])

    modulus = quote α*t.^(-β) end
    model = RheoModelClass(name = "testmodel", p = [:α, :β], G = modulus, info="none")

    data0 = RheoTimeData(t = t, ϵ = loading, σ = exact_response)

    init_params = (α=1.3, β=0.7)
    modelout = modelstepfit(data0, model, strain_imposed, p0=init_params, lo=(α=0.2, β=0.2), hi=(α=3.5, β=1.5))

    found_params = modelout.params

    isapprox(collect(values(actual_params)), collect(values(found_params)), atol=tol)
end
@test _modelstepfit_sing_relax(tol)

function _modelstepfit_sing_creep(tol)
    dt = 0.01
    t = Vector{RheoFloat}(0.0:dt:20.0)
    actual_params = (α=2.0, β=0.5)
    loading = ones(length(t))
    exact_response = loading[1]*actual_params[1]*t.^(-actual_params[2])

    modulus = quote α*t.^(-β) end
    model = RheoModelClass(name = "testmodel", p = [:α, :β], J = modulus, info="none")

    data0 = RheoTimeData(t = t, ϵ = exact_response, σ = loading)

    init_params = (α=1.3, β=0.7)
    modelout = modelstepfit(data0, model, stress_imposed, p0=init_params, lo=(α=0.2, β=0.2), hi=(α=3.5, β=1.5))

    found_params = modelout.params

    isapprox(collect(values(actual_params)), collect(values(found_params)), atol=tol)
end
@test _modelstepfit_sing_creep(tol)

function _modelstepfit_nonsing_nobounds_relax_weighted_selftest(tol)
    dt = 0.01
    t = Vector{RheoFloat}(0.0:dt:20.0)
    actual_params = (α=2.0, β=0.5)
    loading = ones(length(t))
    exact_response = loading[1]*actual_params[1]*exp.(-t/actual_params[2])

    modulus = quote α*exp.(-t/β) end
    model = RheoModelClass(name = "testmodel", p = [:α, :β], G = modulus, info="none")

    data0 = RheoTimeData(t = t, ϵ = loading, σ = exact_response)

    init_params = (α=1.3, β=0.7)
    modelout = modelstepfit(data0, model, strain_imposed; p0=init_params, weights=collect(Integer, 1:length(t)))

    found_params = modelout.params

    isapprox(collect(values(actual_params)), collect(values(found_params)), atol=tol)
end
@test _modelstepfit_nonsing_nobounds_relax_weighted_selftest(tol)

function _modelstepfit_sing_relax_weighted_selftest(tol)
    dt = 0.01
    t = Vector{RheoFloat}(0.0:dt:20.0)
    actual_params = (α=2.0, β=0.5)
    loading = ones(length(t))
    exact_response = loading[1]*actual_params[1]*t.^(-actual_params[2])

    modulus = quote α*t.^(-β) end
    model = RheoModelClass(name = "testmodel", p = [:α, :β], G = modulus, info="none")

    data0 = RheoTimeData(t = t, ϵ = loading, σ = exact_response)

    init_params = (α=1.3, β=0.7)
    modelout = modelstepfit(data0, model, strain_imposed, p0=init_params, lo=(α=0.2, β=0.2), hi=(α=3.5, β=1.5), weights=collect(Integer, 1:length(t)))

    found_params = modelout.params

    isapprox(collect(values(actual_params)), collect(values(found_params)), atol=tol)
end
@test _modelstepfit_sing_relax_weighted_selftest(tol)

function _modelstepfit_nonsing_nobounds_relax_weighted_closeness(tol)
    dt = 0.01
    t = Vector{RheoFloat}(0.0:dt:20.0)
    actual_params = (α=2.0, β=0.5)
    loading = ones(length(t))
    exact_response = loading[1]*actual_params[1]*exp.(-t/actual_params[2])

    modulus = quote α*exp.(-t/β) end
    model = RheoModelClass(name = "testmodel", p = [:α, :β], G = modulus, info="none")

    data0 = RheoTimeData(t = t, ϵ = loading, σ = exact_response)
    indices = indexweight(data0; elperiods = [-2, 2], time_boundaries = [0.0, 5.0, 20.0])

    init_params = (α=1.3, β=0.7)
    modelout_weighted = modelstepfit(data0, model, strain_imposed; p0=init_params, weights=indices)
    modelout = modelstepfit(data0, model, strain_imposed; p0=init_params)

    found_params_weighted = modelout_weighted.params
    found_params = modelout.params

    isapprox(collect(values(found_params_weighted)), collect(values(found_params)), atol=tol)
end
@test _modelstepfit_nonsing_nobounds_relax_weighted_closeness(tol)

function _modelstepfit_sing_relax_weighted_closeness(tol)
    dt = 0.01
    t = Vector{RheoFloat}(0.0:dt:20.0)
    actual_params = (α=2.0, β=0.5)
    loading = ones(length(t))
    exact_response = loading[1]*actual_params[1]*t.^(-actual_params[2])

    modulus = quote α*t.^(-β) end
    model = RheoModelClass(name = "testmodel", p = [:α, :β], G = modulus, info="none")

    data0 = RheoTimeData(t = t, ϵ = loading, σ = exact_response)
    indices = indexweight(data0; elperiods = [-2, 2], time_boundaries = [0.0, 5.0, 20.0])

    init_params = (α=1.3, β=0.7)
    modelout_weighted = modelstepfit(data0, model, strain_imposed, p0=init_params, lo=(α=0.2, β=0.2), hi=(α=3.5, β=1.5), weights=indices)
    modelout = modelstepfit(data0, model, strain_imposed, p0=init_params, lo=(α=0.2, β=0.2), hi=(α=3.5, β=1.5))

    found_params_weighted = modelout_weighted.params
    found_params = modelout.params

    isapprox(collect(values(found_params_weighted)), collect(values(found_params)), atol=tol)
end
@test _modelstepfit_sing_relax_weighted_closeness(tol)

function _modelsteppredict_nonsing_relax(tol)
    dt = 0.01
    t = Vector{RheoFloat}(0.0:dt:20.0)
    exact_response = 1.0*exp.(-t/1.0)
    loading = ones(RheoFloat, length(t))

    modulus = quote α*exp.(-t/β) end
    modelclass = RheoModelClass(name = "testmodel", p = [:α, :β], G = modulus, info="none")
    model = RheoModel(modelclass, α=1.0, β=1.0)

    data0 = RheoTimeData(t = t, ϵ = loading)

    computed_response = modelsteppredict(data0, model)

    all(i -> isapprox(exact_response[i], computed_response.σ[i], atol=tol), eachindex(exact_response))
end
@test _modelsteppredict_nonsing_relax(tol)

function _modelsteppredict_nonsing_creep(tol)
    dt = 0.01
    t = Vector{RheoFloat}(0.0:dt:20.0)
    exact_response = 1.0*exp.(-t/1.0)
    loading = ones(RheoFloat, length(t))

    modulus = quote α*exp.(-t/β) end
    modelclass = RheoModelClass(name = "testmodel", p = [:α, :β], J = modulus, info="none")
    model = RheoModel(modelclass, α=1.0, β=1.0)

    data0 = RheoTimeData(t = t, σ = loading)

    computed_response = modelsteppredict(data0, model)

    all(i -> isapprox(exact_response[i], computed_response.ϵ[i], atol=tol), eachindex(exact_response))
end
@test _modelsteppredict_nonsing_creep(tol)

function _modelsteppredict_nonsing_shifted_relax(tol)
    dt = 0.01
    t = Vector{RheoFloat}(0.0:dt:20.0)
    exact_response = [i>=5.0 ? 1.0*exp(-(i-5.0)) : 0.0 for i in t]
    loading = ones(RheoFloat, length(t))

    modulus = quote α*exp.(-t/β) end
    modelclass = RheoModelClass(name = "testmodel", p = [:α, :β], G = modulus, info="none")
    model = RheoModel(modelclass, α=1.0, β=1.0)

    data0 = RheoTimeData(t = t, ϵ = loading)

    computed_response = modelsteppredict(data0, model; step_on=5.0)

    all(i -> isapprox(exact_response[i], computed_response.σ[i], atol=tol), eachindex(exact_response))
end
@test _modelsteppredict_nonsing_shifted_relax(tol)

function _modelsteppredict_nonsing_shifted_creep(tol)
    dt = 0.01
    t = Vector{RheoFloat}(0.0:dt:20.0)
    exact_response = [i>=5.0 ? 1.0*exp(-(i-5.0)) : 0.0 for i in t]
    loading = ones(RheoFloat, length(t))

    modulus = quote α*exp.(-t/β) end
    modelclass = RheoModelClass(name = "testmodel", p = [:α, :β], J = modulus, info="none")
    model = RheoModel(modelclass, α=1.0, β=1.0)

    data0 = RheoTimeData(t = t, σ = loading)

    computed_response = modelsteppredict(data0, model; step_on=5.0)

    all(i -> isapprox(exact_response[i], computed_response.ϵ[i], atol=tol), eachindex(exact_response))
end
@test _modelsteppredict_nonsing_shifted_creep(tol)

function _modelsteppredict_sing_relax(tol)
    dt = 0.01
    t = Vector{RheoFloat}(0.0:dt:20.0)
    exact_response = t.^(-0.5)
    loading = ones(RheoFloat, length(t))

    modulus = quote α*t.^(-β) end
    modelclass = RheoModelClass(name = "testmodel", p = [:α, :β], G = modulus, info="none")
    model = RheoModel(modelclass, α=1.0, β=0.5)
    data0 = RheoTimeData(t = t, ϵ = loading)

    computed_response = modelsteppredict(data0, model)

    all(i -> isapprox(exact_response[i], computed_response.σ[i], atol=tol), 1:length(t))
end
@test _modelsteppredict_sing_relax(tol)

function _modelsteppredict_sing_creep(tol)
    dt = 0.01
    t = Vector{RheoFloat}(0.0:dt:20.0)
    exact_response = t.^(-0.5)
    loading = ones(RheoFloat, length(t))

    modulus = quote α*t.^(-β) end
    modelclass = RheoModelClass(name = "testmodel", p = [:α, :β], J = modulus, info="none")
    model = RheoModel(modelclass, α=1.0, β=0.5)
    data0 = RheoTimeData(t = t, σ = loading)

    computed_response = modelsteppredict(data0, model)

    all(i -> isapprox(exact_response[i], computed_response.ϵ[i], atol=tol), 1:length(t))
end
@test _modelsteppredict_sing_creep(tol)

function _modelsteppredict_sing_shifted_relax(tol)
    dt = 0.01
    t = Vector{RheoFloat}(0.0:dt:20.0)
    exact_response = [i>=5.0 ? (i - 5.0).^(-0.5) : 0.0 for i in t]
    loading = ones(RheoFloat, length(t))

    modulus = quote α*t.^(-β) end
    modelclass = RheoModelClass(name = "testmodel", p = [:α, :β], G = modulus, info="none")
    model = RheoModel(modelclass, α=1.0, β=0.5)
    data0 = RheoTimeData(t = t, ϵ = loading)

    computed_response = modelsteppredict(data0, model, step_on=5.0)

    stepon_el = RHEOS.closestindex(t, 5.0)

    all(i -> isapprox(exact_response[i], computed_response.σ[i], atol=tol), 1:length(t))
end
@test _modelsteppredict_sing_shifted_relax(tol)

function _modelsteppredict_sing_shifted_creep(tol)
    dt = 0.01
    t = Vector{RheoFloat}(0.0:dt:20.0)
    exact_response = [i>=5.0 ? (i - 5.0).^(-0.5) : 0.0 for i in t]
    loading = ones(RheoFloat, length(t))

    modulus = quote α*t.^(-β) end
    modelclass = RheoModelClass(name = "testmodel", p = [:α, :β], J = modulus, info="none")
    model = RheoModel(modelclass, α=1.0, β=0.5)
    data0 = RheoTimeData(t = t, σ = loading)

    computed_response = modelsteppredict(data0, model, step_on=5.0)

    stepon_el = RHEOS.closestindex(t, 5.0)

    all(i -> isapprox(exact_response[i], computed_response.ϵ[i], atol=tol), 1:length(t))
end
@test _modelsteppredict_sing_shifted_creep(tol)

function _obj_dynamic(tol)
    params = [1.0, 0.1]
    ω = Vector{RheoFloat}(0.0:0.1:100.0)

    dataGp = 2*ω
    dataGpp = 3*ω

    moduliGp(ω, params) = params[1]*2*ω .+ params[2]
    moduliGpp(ω, params) = params[1]*3*ω .+ params[2]

    cost = RHEOS.obj_dynamic(params, nothing, ω, dataGp, dataGpp, moduliGp, moduliGpp)

    costGp = sum((dataGp - moduliGp(ω, params)).^2)
    costGpp = sum((dataGpp - moduliGpp(ω, params)).^2)

    isapprox(cost, (costGp + costGpp))
end
@test _obj_dynamic(tol)

function _obj_dynamic_mean(tol)
    params = [1.0, 0.1]
    ω = Vector{RheoFloat}(0.0:0.1:100.0)

    dataGp = 2*ω
    dataGpp = 3*ω

    meanGp = sum(dataGp)/length(dataGp)
    meanGpp = sum(dataGpp)/length(dataGpp)

    moduliGp(ω, params) = params[1]*2*ω .+ params[2]
    moduliGpp(ω, params) = params[1]*3*ω .+ params[2]

    cost = RHEOS.obj_dynamic_mean(params, nothing, ω, dataGp, dataGpp, moduliGp, moduliGpp, meanGp, meanGpp)

    costGp = sum(((dataGp - moduliGp(ω, params))./meanGp).^2)
    costGpp = sum(((dataGpp - moduliGpp(ω, params))./meanGpp).^2)

    isapprox(cost, (costGp + costGpp))
end
@test _obj_dynamic_mean(tol)

function _obj_dynamic_log(tol)
    params = [1.0, 0.2]
    ω = Vector{RheoFloat}(0.0:0.1:100.0)

    dataGp = 2*ω .+ 0.1
    dataGpp = 3*ω .+ 0.1

    moduliGp(ω, params) = params[1]*2*ω .+ params[2]
    moduliGpp(ω, params) = params[1]*3*ω .+ params[2]

    cost = RHEOS.obj_dynamic_log(params, nothing, ω, dataGp, dataGpp, moduliGp, moduliGpp)

    costGp = sum(log.(dataGp./moduliGp(ω, params)).^2)
    costGpp = sum(log.(dataGpp./moduliGpp(ω, params)).^2)

    isapprox(cost, (costGp + costGpp))
end
@test _obj_dynamic_log(tol)

function _obj_dynamic_local(tol)
    params = [1.0, 0.2]
    ω = Vector{RheoFloat}(0.0:0.1:100.0)

    dataGp = 2*ω .+ 0.1
    dataGpp = 3*ω .+ 0.1

    moduliGp(ω, params) = params[1]*2*ω .+ params[2]
    moduliGpp(ω, params) = params[1]*3*ω .+ params[2]

    cost = RHEOS.obj_dynamic_local(params, nothing, ω, dataGp, dataGpp, moduliGp, moduliGpp)

    costGp = sum((1.0 .- moduliGp(ω, params)./dataGp).^2)
    costGpp = sum((1.0 .- moduliGpp(ω, params)./dataGpp).^2)

    isapprox(cost, (costGp + costGpp))
end
@test _obj_dynamic_local(tol)

function _obj_dynamic_manual(tol)
    params = [1.0, 0.2]
    ω = Vector{RheoFloat}(0.0:0.1:100.0)

    dataGp = 2*ω .+ 0.1
    dataGpp = 3*ω .+ 0.1

    moduliGp(ω, params) = params[1]*2*ω .+ params[2]
    moduliGpp(ω, params) = params[1]*3*ω .+ params[2]

    weights = [10.125, 21.2532602]

    cost = RHEOS.obj_dynamic_manual(params, nothing, ω, dataGp, dataGpp, moduliGp, moduliGpp, weights)

    costGp = sum((dataGp - moduliGp(ω, params)).^2)
    costGpp = sum((dataGpp - moduliGpp(ω, params)).^2)

    isapprox(cost, (costGp*weights[1] + costGpp*weights[2]))
end
@test _obj_dynamic_manual(tol)

function _dynamicmodelfit(tol)
    ω = Vector{RheoFloat}(0.0:0.01:20.0)
    actual_model = RheoModel(SLS_Zener, η=5.0, kᵦ=2.5, kᵧ=7.5)

    dataGp = storagemod(actual_model,ω)
    dataGpp = lossmod(actual_model,ω)
    data0 = RheoFreqData(ω = ω, Gp = dataGp, Gpp = dataGpp)

    init_params = (η=2.0, kᵦ=10.5, kᵧ=3.0)
    modelout = dynamicmodelfit(data0, SLS_Zener, p0=init_params, lo=(η=0.1, kᵦ=0.2, kᵧ=0.3), hi=(η=12.0, kᵦ=Inf, kᵧ=27.0), weights="none", rel_tol=1e-6)

    found_params = modelout.params
    actual_params = actual_model.params

    isapprox(collect(values(found_params)), collect(values(actual_params)), atol = tol)
end
@test _dynamicmodelfit(tol)

function _dynamicmodelfit_noinit(tol)
    ω = Vector{RheoFloat}(0.0:0.01:20.0)
    actual_model = RheoModel(SLS_Zener, η=5.0, kᵦ=2.5, kᵧ=7.5)

    dataGp = storagemod(actual_model,ω)
    dataGpp = lossmod(actual_model,ω)
    data0 = RheoFreqData(ω = ω, Gp = dataGp, Gpp = dataGpp)

    modelout = dynamicmodelfit(data0, SLS_Zener, lo=(η=0.1, kᵦ=0.2, kᵧ=0.3), hi=(η=12.0, kᵦ=Inf, kᵧ=27.0), weights="none", rel_tol=1e-6)

    found_params = modelout.params
    actual_params = actual_model.params

    isapprox(collect(values(found_params)), collect(values(actual_params)), atol = tol)
end
@test _dynamicmodelfit_noinit(tol)

function _dynamicmodelfit_nobounds(tol)
    ω = Vector{RheoFloat}(0.0:0.01:20.0)
    actual_model = RheoModel(SLS_Zener, η=5.0, kᵦ=2.5, kᵧ=7.5)

    dataGp = storagemod(actual_model,ω)
    dataGpp = lossmod(actual_model,ω)
    data0 = RheoFreqData(ω = ω, Gp = dataGp, Gpp = dataGpp)

    init_params = (η=2.0, kᵦ=10.5, kᵧ=3.0)
    modelout = dynamicmodelfit(data0, SLS_Zener, p0=init_params, weights="none", rel_tol=1e-6)

    found_params = modelout.params
    actual_params = actual_model.params

    isapprox(collect(values(found_params)), collect(values(actual_params)), atol = tol)
end
@test _dynamicmodelfit_nobounds(tol)

function _dynamicmodelfit_mean(tol)
    ω = Vector{RheoFloat}(0.0:0.01:20.0)
    actual_model = RheoModel(SLS_Zener, η=5.0, kᵦ=2.5, kᵧ=7.5)

    dataGp = storagemod(actual_model,ω)
    dataGpp = lossmod(actual_model,ω)
    data0 = RheoFreqData(ω = ω, Gp = dataGp, Gpp = dataGpp)

    init_params = (η=2.0, kᵦ=10.5, kᵧ=3.0)
    modelout = dynamicmodelfit(data0, SLS_Zener, p0=init_params, lo=(η=0.1, kᵦ=0.2, kᵧ=0.3), hi=(η=12.0, kᵦ=Inf, kᵧ=27.0), weights="mean", rel_tol=1e-6)

    found_params = modelout.params
    actual_params = actual_model.params

    isapprox(collect(values(found_params)), collect(values(actual_params)), atol = tol)
end
@test _dynamicmodelfit_mean(tol)

function _dynamicmodelfit_log(tol)
    ω = Vector{RheoFloat}(0.01:0.01:20.0)
    actual_model = RheoModel(SLS_Zener, η=5.0, kᵦ=2.5, kᵧ=7.5)

    dataGp = storagemod(actual_model,ω)
    dataGpp = lossmod(actual_model,ω)
    data0 = RheoFreqData(ω = ω, Gp = dataGp, Gpp = dataGpp)

    init_params = (η=2.0, kᵦ=10.5, kᵧ=3.0)
    modelout = dynamicmodelfit(data0, SLS_Zener, p0=init_params, lo=(η=0.1, kᵦ=0.2, kᵧ=0.3), hi=(η=12.0, kᵦ=Inf, kᵧ=27.0), weights="log", rel_tol=1e-6)

    found_params = modelout.params
    actual_params = actual_model.params

    isapprox(collect(values(found_params)), collect(values(actual_params)), atol = tol)
end
@test _dynamicmodelfit_log(tol)

function _dynamicmodelfit_local(tol)
    ω = Vector{RheoFloat}(0.01:0.01:20.0)
    actual_model = RheoModel(SLS_Zener, η=5.0, kᵦ=2.5, kᵧ=7.5)

    dataGp = storagemod(actual_model,ω)
    dataGpp = lossmod(actual_model,ω)
    data0 = RheoFreqData(ω = ω, Gp = dataGp, Gpp = dataGpp)

    init_params = (η=2.0, kᵦ=10.5, kᵧ=3.0)
    modelout = dynamicmodelfit(data0, SLS_Zener, p0=init_params, lo=(η=0.1, kᵦ=0.2, kᵧ=0.3), hi=(η=12.0, kᵦ=Inf, kᵧ=27.0), weights="local", rel_tol=1e-6)

    found_params = modelout.params
    actual_params = actual_model.params

    isapprox(collect(values(found_params)), collect(values(actual_params)), atol = tol)
end
@test _dynamicmodelfit_local(tol)

function _dynamicmodelfit_manualweights(tol)
    ω = Vector{RheoFloat}(0.01:0.01:20.0)
    actual_model = RheoModel(SLS_Zener, η=5.0, kᵦ=2.5, kᵧ=7.5)

    dataGp = storagemod(actual_model,ω)
    dataGpp = lossmod(actual_model,ω)
    data0 = RheoFreqData(ω = ω, Gp = dataGp, Gpp = dataGpp)

    init_params = (η=2.0, kᵦ=10.5, kᵧ=3.0)
    modelout = dynamicmodelfit(data0, SLS_Zener, p0=init_params, lo=(η=0.1, kᵦ=0.2, kᵧ=0.3), hi=(η=12.0, kᵦ=Inf, kᵧ=27.0), weights=[0.5, 1.5], rel_tol=1e-6)

    found_params = modelout.params
    actual_params = actual_model.params

    isapprox(collect(values(found_params)), collect(values(actual_params)), atol = tol)
end
@test _dynamicmodelfit_manualweights(tol)

function _dynamicmodelpredict(tol)
    ω = Vector{RheoFloat}(0.0:0.01:20.0)
    actual_model = RheoModel(SLS_Zener, η=5.0, kᵦ=2.5, kᵧ=7.5)

    dataGp = storagemod(actual_model,ω)
    dataGpp = lossmod(actual_model,ω)
    data0 = RheoFreqData(ω = ω, Gp = dataGp, Gpp = dataGpp)

    datapredicted = dynamicmodelpredict(data0, actual_model)

    test1 = isapprox(dataGp, datapredicted.Gp)
    test2 = isapprox(dataGpp, datapredicted.Gpp)
    test3 = isapprox(ω, datapredicted.ω)

    test1 && test2 && test3
end
@test _dynamicmodelpredict(tol)
