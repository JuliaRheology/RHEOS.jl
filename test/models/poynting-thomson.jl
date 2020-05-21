function Fract_PT_G_reduce_SLS()
    dt = 0.001
    # NaN at 0.0 due to InverseLaplace so start at dt
    t = collect(dt:dt:10.0)

    fracpt = Fract_PT.Ga(t, [2.0, 1.0, 1.0, 0.0, 0.5, 0.0])
    stanSLS = SLS_PT.Ga(t, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracpt[i], stanSLS[i]), eachindex(t))
end
@test Fract_PT_G_reduce_SLS()

function FractSLS_PT_G_reduce_SLS()
    dt = 0.001
    t = collect(0.0:dt:10.0)

    fracSLS = Ga(FractSLS_PT)(t, [2.0, 1.0, 1.0, 0.5])
    stanSLS = SLS_PT.Ga(t, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracSLS[i], stanSLS[i]), eachindex(t))
end
@test FractSLS_PT_G_reduce_SLS()

function Fract_PT_J_reduce_SLS(tol)
    dt = 0.001
    t = collect(0.0:dt:10.0)

    fracpt = Ja(Fract_PT)(t, [2.0, 0.999, 1.0, 0.0, 0.5, 0.0])
    stanSLS = SLS_PT.Ja(t, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracpt[i], stanSLS[i], atol=tol), eachindex(t))
end
@test Fract_PT_J_reduce_SLS(tol)

function FractSLS_PT_J_reduce_SLS(tol)
    dt = 0.001
    t = collect(0.0:dt:10.0)

    fracSLS = Ja(FractSLS_PT)(t, [2.0, 0.999, 1.0, 0.5])
    stanSLS = SLS_PT.Ja(t, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracSLS[i], stanSLS[i], atol=tol), eachindex(t))
end
@test FractSLS_PT_J_reduce_SLS(tol)

function Fract_PT_Gp_reduce_SLS()
    ω_step = 0.01
    chirp = collect(0.0:ω_step:10000.0)

    fracpt = Fract_PT.Gpa(chirp, [2.0, 1.0, 1.0, 0.0, 0.5, 0.0])
    stanSLS = SLS_PT.Gpa(chirp, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracpt[i], stanSLS[i]), eachindex(chirp))
end
@test Fract_PT_Gp_reduce_SLS()

function FractSLS_PT_Gp_reduce_SLS()
    ω_step = 0.01
    chirp = collect(0.0:ω_step:10000.0)

    fracSLS = FractSLS_PT.Gpa(chirp, [2.0, 1.0, 1.0, 0.5])
    stanSLS = SLS_PT.Gpa(chirp, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracSLS[i], stanSLS[i]), eachindex(chirp))
end
@test FractSLS_PT_Gp_reduce_SLS()

function Fract_PT_Gpp_reduce_SLS()
    ω_step = 0.01
    chirp = collect(0.0:ω_step:10000.0)

    fracpt = Fract_PT.Gppa(chirp, [2.0, 1.0, 1.0, 0.0, 0.5, 0.0])
    stanSLS = SLS_PT.Gppa(chirp, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracpt[i], stanSLS[i]), eachindex(chirp))
end
@test Fract_PT_Gpp_reduce_SLS()

function FractSLS_PT_Gpp_reduce_SLS()
    ω_step = 0.01
    chirp = collect(0.0:ω_step:10000.0)

    fracSLS = FractSLS_PT.Gppa(chirp, [2.0, 1.0, 1.0, 0.5])
    stanSLS = SLS_PT.Gppa(chirp, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracSLS[i], stanSLS[i]), eachindex(chirp))
end
@test FractSLS_PT_Gpp_reduce_SLS()

function Fract_PT_G_reduce_Jeffreys(tol)
    dt = 0.001
    # NaN at 0.0 due to Gamma function so start at dt
    t = collect(dt:dt:10.0)

    fracpt = Fract_PT.Ga(t, [2.0, 1.0, 1.0, 0.0, 0.5, 1.0])
    stanJef = Jeffreys_PT.Ga(t, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracpt[i], stanJef[i], atol=tol), eachindex(t))
end
@test Fract_PT_G_reduce_Jeffreys(tol)

function FractJeffreys_PT_G_reduce_Jeffreys()
    dt = 0.001
    t = collect(0.0:dt:10.0)

    fracJef = Ga(FractJeffreys_PT)(t, [2.0, 1.0, 0.0, 0.5])
    stanJef = Jeffreys_PT.Ga(t, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracJef[i], stanJef[i]), eachindex(t))
end
@test FractJeffreys_PT_G_reduce_Jeffreys()

function Fract_PT_J_reduce_Jeffreys(tol)
    dt = 0.001
    t = collect(0.0:dt:10.0)

    fracpt = Ja(Fract_PT)(t, [2.0, 0.999, 1.0, 0.0, 0.5, 1.0])
    stanJef = Jeffreys_PT.Ja(t, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracpt[i], stanJef[i], atol=tol), eachindex(t))
end
@test Fract_PT_J_reduce_Jeffreys(tol)

function FractJeffreys_PT_J_reduce_Jeffreys(tol)
    dt = 0.001
    t = collect(0.0:dt:10.0)

    fracJef = Ja(FractJeffreys_PT)(t, [2.0, 1.0, 0.001, 0.5])
    stanJef = Jeffreys_PT.Ja(t, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracJef[i], stanJef[i], atol=tol), eachindex(t))
end
@test FractJeffreys_PT_J_reduce_Jeffreys(tol)

function Fract_PT_Gp_reduce_Jeffreys()
    ω_step = 0.01
    chirp = collect(0.0:ω_step:10000.0)

    fracpt = Fract_PT.Gpa(chirp, [2.0, 1.0, 1.0, 0.0, 0.5, 1.0])
    stanJef = Jeffreys_PT.Gpa(chirp, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracpt[i], stanJef[i]), eachindex(chirp))
end
@test Fract_PT_Gp_reduce_Jeffreys()

function FractJeffreys_PT_Gp_reduce_Jeffreys()
    ω_step = 0.01
    chirp = collect(0.0:ω_step:10000.0)

    fracJef = FractJeffreys_PT.Gpa(chirp, [2.0, 1.0, 0.0, 0.5])
    stanJef = Jeffreys_PT.Gpa(chirp, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracJef[i], stanJef[i]), eachindex(chirp))
end
@test FractJeffreys_PT_Gp_reduce_Jeffreys()

function Fract_PT_Gpp_reduce_Jeffreys()
    ω_step = 0.01
    chirp = collect(0.0:ω_step:10000.0)

    fracpt = Fract_PT.Gppa(chirp, [2.0, 1.0, 1.0, 0.0, 0.5, 1.0])
    stanJef = Jeffreys_PT.Gppa(chirp, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracpt[i], stanJef[i]), eachindex(chirp))
end
@test Fract_PT_Gpp_reduce_Jeffreys()

function FractJeffreys_PT_Gpp_reduce_Jeffreys()
    ω_step = 0.01
    chirp = collect(0.0:ω_step:10000.0)

    fracJef = FractJeffreys_PT.Gppa(chirp, [2.0, 1.0, 0.0, 0.5])
    stanJef = Jeffreys_PT.Gppa(chirp, [2.0, 1.0, 0.5])

    all(i -> isapprox(fracJef[i], stanJef[i]), eachindex(chirp))
end
@test FractJeffreys_PT_Gpp_reduce_Jeffreys()