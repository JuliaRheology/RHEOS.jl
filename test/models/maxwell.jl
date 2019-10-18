function Fract_Maxwell_G_reduce()
    dt = 0.001
    t = collect(0.0:dt:10.0)

    fm_model = Fract_Maxwell.Ga(t, [1.0, 1.0, 1.0, 0.0])
    maxwell = Maxwell.Ga(t, [1.0, 1.0])

    all(i -> isapprox(fm_model[i], maxwell[i]), eachindex(t))
end
@test Fract_Maxwell_G_reduce()

function FractS_Maxwell_G_reduce()
    dt = 0.001
    t = collect(0.0:dt:10.0)

    fms_model = FractS_Maxwell.Ga(t, [1.0, 1.0, 1.0])
    maxwell = Maxwell.Ga(t, [1.0, 1.0])

    all(i -> isapprox(fms_model[i], maxwell[i]), eachindex(t))
end
@test FractS_Maxwell_G_reduce()

function FractD_Maxwell_G_reduce()
    dt = 0.001
    t = collect(0.0:dt:10.0)

    fmd_model = FractD_Maxwell.Ga(t, [1.0, 1.0, 0.0])
    maxwell = Maxwell.Ga(t, [1.0, 1.0])

    all(i -> isapprox(fmd_model[i], maxwell[i]), eachindex(t))
end
@test FractD_Maxwell_G_reduce()

function Fract_Maxwell_J_reduce()
    dt = 0.001
    t = collect(0.0:dt:10.0)

    fm_model = Fract_Maxwell.Ja(t, [1.0, 1.0, 1.0, 0.0])
    maxwell = Maxwell.Ja(t, [1.0, 1.0])

    all(i -> isapprox(fm_model[i], maxwell[i]), eachindex(t))
end
@test Fract_Maxwell_J_reduce()

function FractS_Maxwell_J_reduce()
    dt = 0.001
    t = collect(0.0:dt:10.0)

    fms_model = FractS_Maxwell.Ja(t, [1.0, 1.0, 1.0])
    maxwell = Maxwell.Ja(t, [1.0, 1.0])

    all(i -> isapprox(fms_model[i], maxwell[i]), eachindex(t))
end
@test FractS_Maxwell_J_reduce()

function FractD_Maxwell_J_reduce()
    dt = 0.001
    t = collect(0.0:dt:10.0)

    fmd_model = FractD_Maxwell.Ja(t, [1.0, 1.0, 0.0])
    maxwell = Maxwell.Ja(t, [1.0, 1.0])

    all(i -> isapprox(fmd_model[i], maxwell[i]), eachindex(t))
end
@test FractD_Maxwell_J_reduce()

function Fract_Maxwell_Gp_reduce()
    ω_step = 0.01;
    chirp = collect(0.0:ω_step:10000.0);

    fm_model = Fract_Maxwell.Gpa(t, [1.0, 1.0, 1.0, 0.0])
    maxwell = Maxwell.Gpa(t, [1.0, 1.0])

    all(i -> isapprox(fm_model[i], maxwell[i]), eachindex(t))
end
@test Fract_Maxwell_Gp_reduce()

function FractS_Maxwell_Gp_reduce()
    ω_step = 0.01;
    chirp = collect(0.0:ω_step:10000.0);

    fms_model = FractS_Maxwell.Gpa(t, [1.0, 1.0, 1.0])
    maxwell = Maxwell.Gpa(t, [1.0, 1.0])

    all(i -> isapprox(fms_model[i], maxwell[i]), eachindex(t))
end
@test FractS_Maxwell_Gp_reduce()

function FractD_Maxwell_Gp_reduce()
    ω_step = 0.01;
    chirp = collect(0.0:ω_step:10000.0);

    fmd_model = FractD_Maxwell.Gpa(t, [1.0, 1.0, 0.0])
    maxwell = Maxwell.Gpa(t, [1.0, 1.0])

    all(i -> isapprox(fmd_model[i], maxwell[i]), eachindex(t))
end
@test FractD_Maxwell_Gp_reduce()

function Fract_Maxwell_Gp_reduce()
    ω_step = 0.01;
    chirp = collect(0.0:ω_step:10000.0);

    fm_model = Fract_Maxwell.Gppa(t, [1.0, 1.0, 1.0, 0.0])
    maxwell = Maxwell.Gppa(t, [1.0, 1.0])

    all(i -> isapprox(fm_model[i], maxwell[i]), eachindex(t))
end
@test Fract_Maxwell_Gp_reduce()

function FractS_Maxwell_Gp_reduce()
    ω_step = 0.01;
    chirp = collect(0.0:ω_step:10000.0);

    fms_model = FractS_Maxwell.Gppa(t, [1.0, 1.0, 1.0])
    maxwell = Maxwell.Gppa(t, [1.0, 1.0])

    all(i -> isapprox(fms_model[i], maxwell[i]), eachindex(t))
end
@test FractS_Maxwell_Gp_reduce()

function FractD_Maxwell_Gp_reduce()
    ω_step = 0.01;
    chirp = collect(0.0:ω_step:10000.0);

    fmd_model = FractD_Maxwell.Gppa(t, [1.0, 1.0, 0.0])
    maxwell = Maxwell.Gppa(t, [1.0, 1.0])

    all(i -> isapprox(fmd_model[i], maxwell[i]), eachindex(t))
end
@test FractD_Maxwell_Gp_reduce()
