println("===============================================")
println("Testing zener.jl")
println("===============================================")

function Fract_Zener_G_reduce_SLS()
    dt = 0.001
    t = collect(0.0:dt:10.0)

    fraczen = relaxmod(Fract_Zener, t, [2.0, 1.0, 1.0, 0.0, 0.5, 0.0])
    stanSLS = relaxmod(SLS_Zener, t, [2.0, 1.0, 0.5])

    all(i -> isapprox(fraczen[i], stanSLS[i]), eachindex(t))
end
@test Fract_Zener_G_reduce_SLS()

function FractSLS_Zener_G_reduce_SLS()
    dt = 0.001
    t = collect(0.0:dt:10.0)

    fracSLS = relaxmod(FractSLS_Zener, t, [2.0, 1.0, 1.0, 0.5])
    stanSLS = relaxmod(SLS_Zener, t, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracSLS[i], stanSLS[i]), eachindex(t))
end
@test FractSLS_Zener_G_reduce_SLS()

function FractSolid_G_reduce_SLS()
    dt = 0.001
    t = collect(0.0:dt:10.0)

    fracSpec = relaxmod(FractSolid, t, [2.0, 1.0, 0.0, 0.5])
    stanSLS = relaxmod(SLS_Zener, t, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracSpec[i], stanSLS[i]), eachindex(t))
end
@test FractSolid_G_reduce_SLS()




function Fract_Zener_J_reduce_SLS()
    dt = 0.001
    # NaN at 0.0 due to InverseLaplace so start at dt
    t = collect(dt:dt:10.0)

    fraczen = creepcomp(Fract_Zener, t, [2.0, 1.0, 1.0, 0.0, 0.5, 0.0])
    stanSLS = creepcomp(SLS_Zener, t, [2.0, 1.0, 0.5])

    all(i -> isapprox(fraczen[i], stanSLS[i]), eachindex(t))
end
@test Fract_Zener_J_reduce_SLS()

function FractSLS_Zener_J_reduce_SLS()
    dt = 0.001
    # NaN at 0.0 due to InverseLaplace so start at dt
    t = collect(dt:dt:10.0)

    fracSLS = creepcomp(FractSLS_Zener, t, [2.0, 1.0, 1.0, 0.5])
    stanSLS = creepcomp(SLS_Zener, t, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracSLS[i], stanSLS[i]), eachindex(t))
end
@test FractSLS_Zener_J_reduce_SLS()

function FractSolid_J_reduce_SLS()
    dt = 0.001
    # NaN at 0.0 due to InverseLaplace so start at dt
    t = collect(dt:dt:10.0)

    fracSpec = creepcomp(FractSolid, t, [2.0, 1.0, 0.0, 0.5])
    stanSLS = creepcomp(SLS_Zener, t, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracSpec[i], stanSLS[i]), eachindex(t))
end
@test FractSolid_J_reduce_SLS()




function Fract_Zener_Gp_reduce_SLS()
    ω_step = 0.01
    chirp = collect(0.0:ω_step:10000.0)

    fraczen = storagemod(Fract_Zener, chirp, [2.0, 1.0, 1.0, 0.0, 0.5, 0.0])
    stanSLS = storagemod(SLS_Zener, chirp, [2.0, 1.0, 0.5])

    all(i -> isapprox(fraczen[i], stanSLS[i]), eachindex(chirp))
end
@test Fract_Zener_Gp_reduce_SLS()

function FractSLS_Zener_Gp_reduce_SLS()
    ω_step = 0.01
    chirp = collect(0.0:ω_step:10000.0)

    fracSLS = storagemod(FractSLS_Zener, chirp, [2.0, 1.0, 1.0, 0.5])
    stanSLS = storagemod(SLS_Zener, chirp, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracSLS[i], stanSLS[i]), eachindex(chirp))
end
@test FractSLS_Zener_Gp_reduce_SLS()

function FractSolid_Gp_reduce_SLS()
    ω_step = 0.01
    chirp = collect(0.0:ω_step:10000.0)

    fracSpec = storagemod(FractSolid, chirp, [2.0, 1.0, 0.0, 0.5])
    stanSLS = storagemod(SLS_Zener, chirp, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracSpec[i], stanSLS[i]), eachindex(chirp))
end
@test FractSolid_Gp_reduce_SLS()




function Fract_Zener_Gpp_reduce_SLS()
    ω_step = 0.01
    chirp = collect(0.0:ω_step:10000.0)

    fraczen = lossmod(Fract_Zener, chirp, [2.0, 1.0, 1.0, 0.0, 0.5, 0.0])
    stanSLS = lossmod(SLS_Zener, chirp, [2.0, 1.0, 0.5])

    all(i -> isapprox(fraczen[i], stanSLS[i]), eachindex(chirp))
end
@test Fract_Zener_Gpp_reduce_SLS()

function FractSLS_Zener_Gpp_reduce_SLS()
    ω_step = 0.01
    chirp = collect(0.0:ω_step:10000.0)

    fracSLS = lossmod(FractSLS_Zener, chirp, [2.0, 1.0, 1.0, 0.5])
    stanSLS = lossmod(SLS_Zener, chirp, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracSLS[i], stanSLS[i]), eachindex(chirp))
end
@test FractSLS_Zener_Gpp_reduce_SLS()

function FractSolid_Gpp_reduce_SLS()
    ω_step = 0.01
    chirp = collect(0.0:ω_step:10000.0)

    fracSpec = lossmod(FractSolid, chirp, [2.0, 1.0, 0.0, 0.5])
    stanSLS = lossmod(SLS_Zener, chirp, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracSpec[i], stanSLS[i]), eachindex(chirp))
end
@test FractSolid_Gpp_reduce_SLS()




function Fract_Zener_G_reduce_Jeffreys()
    dt = 0.001
    # NaN at 0.0 due to Gamma function so start at dt
    t = collect(dt:dt:10.0)

    fraczen = relaxmod(Fract_Zener, t, [2.0, 1.0, 1.0, 0.0, 0.5, 1.0])
    stanJef = relaxmod(Jeffreys_Zener, t, [2.0, 1.0, 0.5])

    all(i -> isapprox(fraczen[i], stanJef[i]), eachindex(t))
end
@test Fract_Zener_G_reduce_Jeffreys()

function FractJeffreys_Zener_G_reduce_Jeffreys()
    dt = 0.001
    # NaN at 0.0 due to Gamma function so start at dt
    t = collect(dt:dt:10.0)

    fracJef = relaxmod(FractJeffreys_Zener, t, [2.0, 1.0, 0.0, 0.5])
    stanJef = relaxmod(Jeffreys_Zener, t, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracJef[i], stanJef[i]), eachindex(t))
end
@test FractJeffreys_Zener_G_reduce_Jeffreys()



function Fract_Zener_J_reduce_Jeffreys()
    dt = 0.001
    # NaN at 0.0 due to InverseLaplace function so start at dt
    t = collect(dt:dt:10.0)

    fraczen = creepcomp(Fract_Zener, t, [2.0, 1.0, 1.0, 0.0, 0.5, 1.0])
    stanJef = creepcomp(Jeffreys_Zener, t, [2.0, 1.0, 0.5])

    all(i -> isapprox(fraczen[i], stanJef[i]), eachindex(t))
end
@test Fract_Zener_J_reduce_Jeffreys()

function FractJeffreys_Zener_J_reduce_Jeffreys()
    dt = 0.001
    # NaN at 0.0 due to InverseLaplace function so start at dt
    t = collect(dt:dt:10.0)

    fracJef = creepcomp(FractJeffreys_Zener, t, [2.0, 1.0, 0.0, 0.5])
    stanJef = creepcomp(Jeffreys_Zener, t, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracJef[i], stanJef[i]), eachindex(t))
end
@test FractJeffreys_Zener_J_reduce_Jeffreys()



function Fract_Zener_Gp_reduce_Jeffreys()
    ω_step = 0.01
    chirp = collect(0.0:ω_step:10000.0)

    fraczen = storagemod(Fract_Zener, chirp, [2.0, 1.0, 1.0, 0.0, 0.5, 1.0])
    stanJef = storagemod(Jeffreys_Zener, chirp, [2.0, 1.0, 0.5])

    all(i -> isapprox(fraczen[i], stanJef[i]), eachindex(chirp))
end
@test Fract_Zener_Gp_reduce_Jeffreys()

function FractJeffreys_Zener_Gp_reduce_Jeffreys()
    ω_step = 0.01
    chirp = collect(0.0:ω_step:10000.0)

    fracJef = storagemod(FractJeffreys_Zener, chirp, [2.0, 1.0, 0.0, 0.5])
    stanJef = storagemod(Jeffreys_Zener, chirp, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracJef[i], stanJef[i]), eachindex(chirp))
end
@test FractJeffreys_Zener_Gp_reduce_Jeffreys()



function Fract_Zener_Gpp_reduce_Jeffreys()
    ω_step = 0.01
    chirp = collect(0.0:ω_step:10000.0)

    fraczen = lossmod(Fract_Zener, chirp, [2.0, 1.0, 1.0, 0.0, 0.5, 1.0])
    stanJef = lossmod(Jeffreys_Zener, chirp, [2.0, 1.0, 0.5])

    all(i -> isapprox(fraczen[i], stanJef[i]), eachindex(chirp))
end
@test Fract_Zener_Gpp_reduce_Jeffreys()

function FractJeffreys_Zener_Gpp_reduce_Jeffreys()
    ω_step = 0.01
    chirp = collect(0.0:ω_step:10000.0)

    fracJef = lossmod(FractJeffreys_Zener, chirp, [2.0, 1.0, 0.0, 0.5])
    stanJef = lossmod(Jeffreys_Zener, chirp, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracJef[i], stanJef[i]), eachindex(chirp))
end
@test FractJeffreys_Zener_Gpp_reduce_Jeffreys()
