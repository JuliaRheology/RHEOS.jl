println("===============================================")
println("Testing poynting-thomson.jl")
println("===============================================")

function Fract_PT_G_reduce_SLS()
    dt = 0.001
    # NaN at 0.0 due to InverseLaplace so start at dt
    t = collect(dt:dt:10.0)

    fracpt = relaxmod(Fract_PT, t, [2.0, 1.0, 1.0, 0.0, 0.5, 0.0])
    stanSLS = relaxmod(SLS_PT, t, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracpt[i], stanSLS[i]), eachindex(t))
end
@test Fract_PT_G_reduce_SLS()

function FractSLS_PT_G_reduce_SLS()
    dt = 0.001
    t = collect(0.0:dt:10.0)

    fracSLS = relaxmod(FractSLS_PT, t, [2.0, 1.0, 1.0, 0.5])
    stanSLS = relaxmod(SLS_PT, t, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracSLS[i], stanSLS[i]), eachindex(t))
end
@test FractSLS_PT_G_reduce_SLS()



function Fract_PT_J_reduce_SLS(tol)
    dt = 0.001
    t = collect(0.0:dt:10.0)

    fracpt = creepcomp(Fract_PT, t, [2.0, 0.999, 1.0, 0.0, 0.5, 0.0])
    stanSLS = creepcomp(SLS_PT, t, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracpt[i], stanSLS[i], atol=tol), eachindex(t))
end
@test Fract_PT_J_reduce_SLS(tol)

function FractSLS_PT_J_reduce_SLS(tol)
    dt = 0.001
    t = collect(0.0:dt:10.0)

    fracSLS = creepcomp(FractSLS_PT, t, [2.0, 0.999, 1.0, 0.5])
    stanSLS = creepcomp(SLS_PT, t, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracSLS[i], stanSLS[i], atol=tol), eachindex(t))
end
@test FractSLS_PT_J_reduce_SLS(tol)



function Fract_PT_Gp_reduce_SLS()
    ω_step = 0.01
    chirp = collect(0.0:ω_step:10000.0)

    fracpt = storagemod(Fract_PT, chirp, [2.0, 1.0, 1.0, 0.0, 0.5, 0.0])
    stanSLS = storagemod(SLS_PT, chirp, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracpt[i], stanSLS[i]), eachindex(chirp))
end
@test Fract_PT_Gp_reduce_SLS()

function FractSLS_PT_Gp_reduce_SLS()
    ω_step = 0.01
    chirp = collect(0.0:ω_step:10000.0)

    fracSLS = storagemod(FractSLS_PT, chirp, [2.0, 1.0, 1.0, 0.5])
    stanSLS = storagemod(SLS_PT, chirp, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracSLS[i], stanSLS[i]), eachindex(chirp))
end
@test FractSLS_PT_Gp_reduce_SLS()



function Fract_PT_Gpp_reduce_SLS()
    ω_step = 0.01
    chirp = collect(0.0:ω_step:10000.0)

    fracpt = lossmod(Fract_PT, chirp, [2.0, 1.0, 1.0, 0.0, 0.5, 0.0])
    stanSLS = lossmod(SLS_PT, chirp, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracpt[i], stanSLS[i]), eachindex(chirp))
end
@test Fract_PT_Gpp_reduce_SLS()

function FractSLS_PT_Gpp_reduce_SLS()
    ω_step = 0.01
    chirp = collect(0.0:ω_step:10000.0)

    fracSLS = lossmod(FractSLS_PT, chirp, [2.0, 1.0, 1.0, 0.5])
    stanSLS = lossmod(SLS_PT, chirp, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracSLS[i], stanSLS[i]), eachindex(chirp))
end
@test FractSLS_PT_Gpp_reduce_SLS()





function Fract_PT_G_reduce_Jeffreys(tol)
    dt = 0.001
    # NaN at 0.0 due to Gamma function so start at dt
    t = collect(dt:dt:10.0)

    fracpt = relaxmod(Fract_PT, t, [2.0, 1.0, 1.0, 0.0, 0.5, 1.0])
    stanJef = relaxmod(Jeffreys_PT, t, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracpt[i], stanJef[i], atol=tol), eachindex(t))
end
@test Fract_PT_G_reduce_Jeffreys(tol)

function FractJeffreys_PT_G_reduce_Jeffreys()
    dt = 0.001
    t = collect(0.0:dt:10.0)

    fracJef = relaxmod(FractJeffreys_PT, t, [2.0, 1.0, 0.0, 0.5])
    stanJef = relaxmod(Jeffreys_PT, t, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracJef[i], stanJef[i]), eachindex(t))
end
@test FractJeffreys_PT_G_reduce_Jeffreys()



function Fract_PT_J_reduce_Jeffreys(tol)
    dt = 0.001
    t = collect(0.0:dt:10.0)

    fracpt = creepcomp(Fract_PT, t, [2.0, 0.999, 1.0, 0.0, 0.5, 1.0])
    stanJef = creepcomp(Jeffreys_PT, t, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracpt[i], stanJef[i], atol=tol), eachindex(t))
end
@test Fract_PT_J_reduce_Jeffreys(tol)

function FractJeffreys_PT_J_reduce_Jeffreys(tol)
    dt = 0.001
    t = collect(0.0:dt:10.0)

    fracJef = creepcomp(FractJeffreys_PT, t, [2.0, 1.0, 0.001, 0.5])
    stanJef = creepcomp(Jeffreys_PT, t, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracJef[i], stanJef[i], atol=tol), eachindex(t))
end
@test FractJeffreys_PT_J_reduce_Jeffreys(tol)



function Fract_PT_Gp_reduce_Jeffreys()
    ω_step = 0.01
    chirp = collect(0.0:ω_step:10000.0)

    fracpt = storagemod(Fract_PT, chirp, [2.0, 1.0, 1.0, 0.0, 0.5, 1.0])
    stanJef = storagemod(Jeffreys_PT, chirp, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracpt[i], stanJef[i]), eachindex(chirp))
end
@test Fract_PT_Gp_reduce_Jeffreys()

function FractJeffreys_PT_Gp_reduce_Jeffreys()
    ω_step = 0.01
    chirp = collect(0.0:ω_step:10000.0)

    fracJef = storagemod(FractJeffreys_PT, chirp, [2.0, 1.0, 0.0, 0.5])
    stanJef = storagemod(Jeffreys_PT, chirp, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracJef[i], stanJef[i]), eachindex(chirp))
end
@test FractJeffreys_PT_Gp_reduce_Jeffreys()



function Fract_PT_Gpp_reduce_Jeffreys()
    ω_step = 0.01
    chirp = collect(0.0:ω_step:10000.0)

    fracpt = lossmod(Fract_PT, chirp, [2.0, 1.0, 1.0, 0.0, 0.5, 1.0])
    stanJef = lossmod(Jeffreys_PT, chirp, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracpt[i], stanJef[i]), eachindex(chirp))
end
@test Fract_PT_Gpp_reduce_Jeffreys()

function FractJeffreys_PT_Gpp_reduce_Jeffreys()
    ω_step = 0.01
    chirp = collect(0.0:ω_step:10000.0)

    fracJef = lossmod(FractJeffreys_PT, chirp, [2.0, 1.0, 0.0, 0.5])
    stanJef = lossmod(Jeffreys_PT, chirp, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracJef[i], stanJef[i]), eachindex(chirp))
end
@test FractJeffreys_PT_Gpp_reduce_Jeffreys()
