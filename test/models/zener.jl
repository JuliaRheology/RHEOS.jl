function Fract_Zener_G_reduce_SLS()
    dt = 0.001
    t = collect(0.0:dt:10.0)

    fraczen = Fract_Zener.Ga(t, [2.0, 1.0, 1.0, 0.0, 0.5, 0.0])
    stanSLS = SLS_Zener.Ga(t, [2.0, 1.0, 0.5])

    all(i -> isapprox(fraczen[i], stanSLS[i]), eachindex(t))
end
@test Fract_Zener_G_reduce_SLS()

function FractSLS_Zener_G_reduce_SLS()
    dt = 0.001
    t = collect(0.0:dt:10.0)

    fracSLS = FractSLS_Zener.Ga(t, [2.0, 1.0, 1.0, 0.5])
    stanSLS = SLS_Zener.Ga(t, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracSLS[i], stanSLS[i]), eachindex(t))
end
@test FractSLS_Zener_G_reduce_SLS()

function FractSolid_G_reduce_SLS()
    dt = 0.001
    t = collect(0.0:dt:10.0)

    fracSpec = FractSolid.Ga(t, [2.0, 1.0, 0.0, 0.5])
    stanSLS = SLS_Zener.Ga(t, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracSpec[i], stanSLS[i]), eachindex(t))
end
@test FractSolid_G_reduce_SLS()

function Fract_Zener_J_reduce_SLS()
    dt = 0.001
    # NaN at 0.0 due to InverseLaplace so start at dt
    t = collect(dt:dt:10.0)

    fraczen = Fract_Zener.Ja(t, [2.0, 1.0, 1.0, 0.0, 0.5, 0.0])
    stanSLS = SLS_Zener.Ja(t, [2.0, 1.0, 0.5])

    all(i -> isapprox(fraczen[i], stanSLS[i]), eachindex(t))
end
@test Fract_Zener_J_reduce_SLS()

function FractSLS_Zener_J_reduce_SLS()
    dt = 0.001
    # NaN at 0.0 due to InverseLaplace so start at dt
    t = collect(dt:dt:10.0)

    fracSLS = FractSLS_Zener.Ja(t, [2.0, 1.0, 1.0, 0.5])
    stanSLS = SLS_Zener.Ja(t, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracSLS[i], stanSLS[i]), eachindex(t))
end
@test FractSLS_Zener_J_reduce_SLS()

function FractSolid_J_reduce_SLS()
    dt = 0.001
    # NaN at 0.0 due to InverseLaplace so start at dt
    t = collect(dt:dt:10.0)

    fracSpec = FractSolid.Ja(t, [2.0, 1.0, 0.0, 0.5])
    stanSLS = SLS_Zener.Ja(t, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracSpec[i], stanSLS[i]), eachindex(t))
end
@test FractSolid_J_reduce_SLS()

function Fract_Zener_Gp_reduce_SLS()
    ω_step = 0.01
    chirp = collect(0.0:ω_step:10000.0)

    fraczen = Fract_Zener.Gpa(chirp, [2.0, 1.0, 1.0, 0.0, 0.5, 0.0])
    stanSLS = SLS_Zener.Gpa(chirp, [2.0, 1.0, 0.5])

    all(i -> isapprox(fraczen[i], stanSLS[i]), eachindex(chirp))
end
@test Fract_Zener_Gp_reduce_SLS()

function FractSLS_Zener_Gp_reduce_SLS()
    ω_step = 0.01
    chirp = collect(0.0:ω_step:10000.0)

    fracSLS = FractSLS_Zener.Gpa(chirp, [2.0, 1.0, 1.0, 0.5])
    stanSLS = SLS_Zener.Gpa(chirp, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracSLS[i], stanSLS[i]), eachindex(chirp))
end
@test FractSLS_Zener_Gp_reduce_SLS()

function FractSolid_Gp_reduce_SLS()
    ω_step = 0.01
    chirp = collect(0.0:ω_step:10000.0)

    fracSpec = FractSolid.Gpa(chirp, [2.0, 1.0, 0.0, 0.5])
    stanSLS = SLS_Zener.Gpa(chirp, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracSpec[i], stanSLS[i]), eachindex(chirp))
end
@test FractSolid_Gp_reduce_SLS()

function Fract_Zener_Gpp_reduce_SLS()
    ω_step = 0.01
    chirp = collect(0.0:ω_step:10000.0)

    fraczen = Fract_Zener.Gppa(chirp, [2.0, 1.0, 1.0, 0.0, 0.5, 0.0])
    stanSLS = SLS_Zener.Gppa(chirp, [2.0, 1.0, 0.5])

    all(i -> isapprox(fraczen[i], stanSLS[i]), eachindex(chirp))
end
@test Fract_Zener_Gpp_reduce_SLS()

function FractSLS_Zener_Gpp_reduce_SLS()
    ω_step = 0.01
    chirp = collect(0.0:ω_step:10000.0)

    fracSLS = FractSLS_Zener.Gppa(chirp, [2.0, 1.0, 1.0, 0.5])
    stanSLS = SLS_Zener.Gppa(chirp, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracSLS[i], stanSLS[i]), eachindex(chirp))
end
@test FractSLS_Zener_Gpp_reduce_SLS()

function FractSolid_Gpp_reduce_SLS()
    ω_step = 0.01
    chirp = collect(0.0:ω_step:10000.0)

    fracSpec = FractSolid.Gppa(chirp, [2.0, 1.0, 0.0, 0.5])
    stanSLS = SLS_Zener.Gppa(chirp, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracSpec[i], stanSLS[i]), eachindex(chirp))
end
@test FractSolid_Gpp_reduce_SLS()

function Fract_Zener_G_reduce_Jeffreys()
    dt = 0.001
    # NaN at 0.0 due to Gamma function so start at dt
    t = collect(dt:dt:10.0)

    fraczen = Fract_Zener.Ga(t, [2.0, 1.0, 1.0, 0.0, 0.5, 1.0])
    stanJef = Jeffreys_Zener.Ga(t, [2.0, 1.0, 0.5])

    all(i -> isapprox(fraczen[i], stanJef[i]), eachindex(t))
end
@test Fract_Zener_G_reduce_Jeffreys()

function FractJeffreys_Zener_G_reduce_Jeffreys()
    dt = 0.001
    # NaN at 0.0 due to Gamma function so start at dt
    t = collect(dt:dt:10.0)

    fracJef = FractJeffreys_Zener.Ga(t, [2.0, 1.0, 0.0, 0.5])
    stanJef = Jeffreys_Zener.Ga(t, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracJef[i], stanJef[i]), eachindex(t))
end
@test FractJeffreys_Zener_G_reduce_Jeffreys()

function Fract_Zener_J_reduce_Jeffreys()
    dt = 0.001
    # NaN at 0.0 due to InverseLaplace function so start at dt
    t = collect(dt:dt:10.0)

    fraczen = Fract_Zener.Ja(t, [2.0, 1.0, 1.0, 0.0, 0.5, 1.0])
    stanJef = Jeffreys_Zener.Ja(t, [2.0, 1.0, 0.5])

    all(i -> isapprox(fraczen[i], stanJef[i]), eachindex(t))
end
@test Fract_Zener_J_reduce_Jeffreys()

function FractJeffreys_Zener_J_reduce_Jeffreys()
    dt = 0.001
    # NaN at 0.0 due to InverseLaplace function so start at dt
    t = collect(dt:dt:10.0)

    fracJef = FractJeffreys_Zener.Ja(t, [2.0, 1.0, 0.0, 0.5])
    stanJef = Jeffreys_Zener.Ja(t, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracJef[i], stanJef[i]), eachindex(t))
end
@test FractJeffreys_Zener_J_reduce_Jeffreys()

function Fract_Zener_Gp_reduce_Jeffreys()
    ω_step = 0.01
    chirp = collect(0.0:ω_step:10000.0)

    fraczen = Fract_Zener.Gpa(chirp, [2.0, 1.0, 1.0, 0.0, 0.5, 1.0])
    stanJef = Jeffreys_Zener.Gpa(chirp, [2.0, 1.0, 0.5])

    all(i -> isapprox(fraczen[i], stanJef[i]), eachindex(chirp))
end
@test Fract_Zener_Gp_reduce_Jeffreys()

function FractJeffreys_Zener_Gp_reduce_Jeffreys()
    ω_step = 0.01
    chirp = collect(0.0:ω_step:10000.0)

    fracJef = FractJeffreys_Zener.Gpa(chirp, [2.0, 1.0, 0.0, 0.5])
    stanJef = Jeffreys_Zener.Gpa(chirp, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracJef[i], stanJef[i]), eachindex(chirp))
end
@test FractJeffreys_Zener_Gp_reduce_Jeffreys()

function Fract_Zener_Gpp_reduce_Jeffreys()
    ω_step = 0.01
    chirp = collect(0.0:ω_step:10000.0)

    fraczen = Fract_Zener.Gppa(chirp, [2.0, 1.0, 1.0, 0.0, 0.5, 1.0])
    stanJef = Jeffreys_Zener.Gppa(chirp, [2.0, 1.0, 0.5])

    all(i -> isapprox(fraczen[i], stanJef[i]), eachindex(chirp))
end
@test Fract_Zener_Gpp_reduce_Jeffreys()

function FractJeffreys_Zener_Gpp_reduce_Jeffreys()
    ω_step = 0.01
    chirp = collect(0.0:ω_step:10000.0)

    fracJef = FractJeffreys_Zener.Gppa(chirp, [2.0, 1.0, 0.0, 0.5])
    stanJef = Jeffreys_Zener.Gppa(chirp, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracJef[i], stanJef[i]), eachindex(chirp))
end
@test FractJeffreys_Zener_Gpp_reduce_Jeffreys()