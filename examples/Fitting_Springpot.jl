#!/usr/bin/env julia
using PyPlot
using Revise
using RHEOS

#------------------------------GENERATE DATA AND FIT----------------------------

sp = RheoModel(SpringPot, (cᵦ = 1e3, β = 0.2))
time = time_line(t_end = 200.0)
step_strain = strainfunction(time,hstep())
sp_pred = modelpredict(step_strain,sp)

# p0 = (cᵦ=1e3, β=0.99)
# lb = (cᵦ=0.0, β=0.01)
# ub = (cᵦ=Inf, β=0.999)
fit_sp = modelfit(sp_pred,SpringPot,strain_imposed; verbose = true)
predict_sp = modelpredict(step_strain,fit_sp)
#---------------------------------------------------

fs = RheoModel(FractionalSpecial, (η=1e4, cᵦ=1e3, β=0.2, k=800))
time = time_line(t_end = 200.0)
step_strain = strainfunction(time,hstep())
fs_pred = modelpredict(step_strain,fs)

p0 = (η=1e4, cᵦ=1e3, β=0.5, k=800)
lb = (η=0.0, cᵦ=0.0, β=0.1, k=0.0)
ub = (η=Inf, cᵦ=Inf, β=0.9, k=Inf)
fit_fs = modelfit(fs_pred,FractionalSpecial,strain_imposed; p0 =p0, lo = lb, hi = ub, verbose = true)
predict_fs = modelpredict(step_strain,fit_fs)


#--------------------------------------------------

fm = RheoModel(FractionalMaxwell, (cₐ=10, a=1., cᵦ=5, β=0.5))
time = time_line(t_end = 200.0)
step_strain = strainfunction(time,hstep(amp = 0.3))
fm_pred = modelpredict(step_strain,fm)

p0 = (cₐ=5, a=0.52, cᵦ=2, β=0.50)
lb = (cₐ=0.0, a=0.01, cᵦ=0.0, β=0.01)
ub = (cₐ=1e3, a=0.99, cᵦ=1e3, β=0.99)
fit_fs = modelfit(fm_pred,FractionalMaxwell,strain_imposed; p0 = p0, lo = lb, hi = ub,verbose = true)
predict_fs = modelpredict(step_strain,fit_fs)


plot(strain_tot.t, strain_tot.ϵ)
using NLopt, Test

count = 0 # keep track of # function evaluations

function myfunc(x::Vector, grad::Vector)
    if length(grad) > 0
        grad[1] = 0
        grad[2] = 0.5/sqrt(x[2])
    end

    global count
    count::Int += 1
    println("f_$count($x)")

    sqrt(x[2])
end

function mycon(x::Vector,i)
    if i == 1
        c = (2*x[1] + 0)^3 - x[2]
    elseif i == 2
        c = (-1*x[1] + 1)^3 - x[2]
    end
end


opt = Opt(:GN_ISRES, 2)
opt.lower_bounds = [-Inf, 0.]
opt.xtol_rel = 1e-5
opt.min_objective = myfunc


for i = 1:1:2
    myconst(x,g,a) = mycon(x,i)
    print(myconst)
    inequality_constraint!(opt,(x,g) ->myconst(x,g,[0]))
end
lo = [0.0,0.0]
ub = [10,10]

lower_bounds!(opt, lo)
upper_bounds!(opt,ub)

(minf,minx,ret) = optimize(opt, [1.234, 5.678])
println("got $minf at $minx after $count iterations (returned $ret)")



function mycon(r, x::Vector,grad)
        r[1] = (2*x[1] + 0)^3 - x[2]
        # r[2] = (-1*x[1] + 1)^3 - x[2]
        # r[3] = -x[2]+0.5
end
lo = [0.0,0.0]
ub = [10,10]


opt = Opt(:GN_ISRES, 2)
opt.lower_bounds = [-Inf, 0.]
opt.xtol_rel = 1e-4
opt.min_objective = myfunc

lower_bounds!(opt, lo)
upper_bounds!(opt,ub)

inequality_constraint!(opt,mycon, 0.0)


(minf,minx,ret) = optimize(opt, [1.234, 5.678])
println("got $minf at $minx after $count iterations (returned $ret)")
