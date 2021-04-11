println("===============================================")
println("Testing definition.jl")
println("===============================================")

function _RheoTimeData_explicit_nolog()
    a = Vector{RheoFloat}(1.0:1.0:3.0)
    b = Vector{RheoFloat}(4.0:1.0:6.0)
    c = Vector{RheoFloat}(7.0:1.0:9.0)

    data = RheoTimeData(a, b, c, nothing)

    data.σ==a && data.ϵ==b && data.t==c && isnothing(data.log)
end
@test _RheoTimeData_explicit_nolog()

function _RheoTimeData_const_nolog()
    eps = Vector{RheoFloat}(1.0:1.0:3.0)
    sig = Vector{RheoFloat}(4.0:1.0:6.0)
    t = Vector{RheoFloat}(7.0:1.0:9.0)

    data = RheoTimeData(ϵ=eps, σ=sig, t=t, savelog=false)
    showlog(data)

    data.σ==sig && data.ϵ==eps && data.t==t && isnothing(data.log)
end
@test _RheoTimeData_const_nolog()

function _RheoTimeData_const_nostrain()
    sig = Vector{RheoFloat}(4.0:1.0:6.0)
    t = Vector{RheoFloat}(7.0:1.0:9.0)

    data = RheoTimeData(σ=sig, t=t)
    showlog(data)

    data.σ==sig && data.ϵ==RheoFloat[] && data.t==t && data.log[1].info.type == stress_only::TimeDataType
end
@test _RheoTimeData_const_nostrain()

function _operators_logs()
    d1=strainfunction(timeline(0:0.5:10),t->exp(-t))
    d2=strainfunction(timeline(0:0.5:10),t->1-exp(-t))
    d=2*d1 - (-d2) + d2

    d3=rheologrun(d.log)

    (d3.ϵ == d.ϵ) && all([ abs(e-2.)<=eps(RheoFloat) for e in d.ϵ ])
end
@test _operators_logs()




function _rheoconv()
    vi64=Int64(1)
    ai64=[vi64,vi64]
    arf=Vector{RheoFloat}([1,2,3])
    rheoconv(RheoFloat(1.0))===rheoconv(vi64) && typeof(rheoconv(ai64)) == Vector{RheoFloat} &&
        arf === rheoconv(arf)  &&  ai64 !== rheoconv(ai64)
end
@test _rheoconv()


function _union_strain1stress2()
    time_instance = timeline(t_start=0.0, t_end=15.0, step=0.2)
    imposed_strain = strainfunction(time_instance, sawtooth(offset=5.0, amp=2, period=5))
    imposed_stress = stressfunction(time_instance, stairs(offset=3.0, amp=0.5, width=4.0))

    combined = imposed_strain|imposed_stress

    (combined.σ==imposed_stress.σ) && (combined.ϵ==imposed_strain.ϵ) && (combined.t==imposed_stress.t)
end
@test _union_strain1stress2()

function _union_stress1strain2()
    time_instance = timeline(t_start=0.0, t_end=15.0, step=0.2)
    imposed_strain = strainfunction(time_instance, sawtooth(offset=5.0, amp=2, period=5))
    imposed_stress = stressfunction(time_instance, stairs(offset=3.0, amp=0.5, width=4.0))

    combined = imposed_stress|imposed_strain

    (combined.σ==imposed_stress.σ) && (combined.ϵ==imposed_strain.ϵ) && (combined.t==imposed_stress.t)
end
@test _union_stress1strain2()

function _freeze_params()
    SLS2_mod = freeze_params( SLS2, G₀=2, η₂=3.5)

    relaxmod(SLS2, 1, [2,1,2,3,3.5]) == relaxmod(SLS2_mod, 1, [1,2,3])
end
@test _freeze_params()




#
#   The moduli functions on arrays are extensively tested as part of the model tests
#


function _scalar_moduli()
    m=RheoModel(Spring, k=2)
    (relaxmod(m))(1) == relaxmod(m, 1) == relaxmod(Spring, k=2, 1) == 2. &&
        (creepcomp(m))(1) == creepcomp(m, 1) == creepcomp(Spring, k=2, 1) == 0.5 &&
        (storagemod(m))(1) == storagemod(m, 1) == storagemod(Spring, k=2, 1) == 2. &&
        (lossmod(m))(1) == lossmod(m, 1) == lossmod(Spring, k=2, 1) == 0.0 &&
        (dynamicmod(m))(1) == dynamicmod(m, 1) == dynamicmod(Spring, k=2, 1) == 2.0 + 0.0*im 
end
@test _scalar_moduli()
