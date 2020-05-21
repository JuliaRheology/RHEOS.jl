function Fract_KelvinVoigt_G_reduce()
    dt = 0.001
    # fractional modulus is NaN not inf at 0.0 so start at dt
    t = collect(dt:dt:10.0)

    fm_model = Fract_KelvinVoigt.Ga(t, [1.0, 1.0, 1.0, 0.0])
    kelvinvoigt = KelvinVoigt.Ga(t, [1.0, 1.0])

    all(i -> isapprox(fm_model[i], kelvinvoigt[i]), eachindex(t))
end
@test Fract_KelvinVoigt_G_reduce()

function FractS_KelvinVoigt_G_reduce()
    dt = 0.001
    # fractional modulus is NaN not inf at 0.0 so start at dt
    t = collect(dt:dt:10.0)

    fms_model = FractS_KelvinVoigt.Ga(t, [1.0, 1.0, 1.0])
    kelvinvoigt = KelvinVoigt.Ga(t, [1.0, 1.0])

    all(i -> isapprox(fms_model[i], kelvinvoigt[i]), eachindex(t))
end
@test FractS_KelvinVoigt_G_reduce()

function FractD_KelvinVoigt_G_reduce()
    dt = 0.001
    t = collect(0.0:dt:10.0)

    fmd_model = FractD_KelvinVoigt.Ga(t, [1.0, 1.0, 0.0])
    kelvinvoigt = KelvinVoigt.Ga(t, [1.0, 1.0])

    all(i -> isapprox(fmd_model[i], kelvinvoigt[i]), eachindex(t))
end
@test FractD_KelvinVoigt_G_reduce()

function Fract_KelvinVoigt_J_reduce(tol)
    dt = 0.001
    t = collect(0.0:dt:10.0)

    fm_model = Ja(Fract_KelvinVoigt)(t, [1.0, 0.99999, 1.0, 0.00001])
    kelvinvoigt = KelvinVoigt.Ja(t, [1.0, 1.0])

    all(i -> isapprox(fm_model[i], kelvinvoigt[i], atol=tol), eachindex(t))
end
@test Fract_KelvinVoigt_J_reduce(tol)

function FractS_KelvinVoigt_J_reduce(tol)
    dt = 0.001
    t = collect(0.0:dt:10.0)

    fms_model = Ja(FractS_KelvinVoigt)(t, [1.0, 0.99999, 1.0])
    kelvinvoigt = KelvinVoigt.Ja(t, [1.0, 1.0])

    all(i -> isapprox(fms_model[i], kelvinvoigt[i], atol=tol), eachindex(t))
end
@test FractS_KelvinVoigt_J_reduce(tol)

function FractD_KelvinVoigt_J_reduce(tol)
    dt = 0.001
    t = collect(0.0:dt:10.0)

    fmd_model = Ja(FractD_KelvinVoigt)(t, [1.0, 1.0, 0.00001])
    kelvinvoigt = KelvinVoigt.Ja(t, [1.0, 1.0])

    all(i -> isapprox(fmd_model[i], kelvinvoigt[i], atol=tol), eachindex(t))
end
@test FractD_KelvinVoigt_J_reduce(tol)

function Fract_KelvinVoigt_Gp_reduce()
    ω_step = 0.01;
    chirp = collect(0.0:ω_step:10000.0);

    fm_model = Fract_KelvinVoigt.Gpa(chirp, [1.0, 1.0, 1.0, 0.0])
    kelvinvoigt = KelvinVoigt.Gpa(chirp, [1.0, 1.0])

    all(i -> isapprox(fm_model[i], kelvinvoigt[i]), eachindex(chirp))
end
@test Fract_KelvinVoigt_Gp_reduce()

function FractS_KelvinVoigt_Gp_reduce()
    ω_step = 0.01;
    chirp = collect(0.0:ω_step:10000.0);

    fms_model = FractS_KelvinVoigt.Gpa(chirp, [1.0, 1.0, 1.0])
    kelvinvoigt = KelvinVoigt.Gpa(chirp, [1.0, 1.0])

    all(i -> isapprox(fms_model[i], kelvinvoigt[i]), eachindex(chirp))
end
@test FractS_KelvinVoigt_Gp_reduce()

function FractD_KelvinVoigt_Gp_reduce()
    ω_step = 0.01;
    chirp = collect(0.0:ω_step:10000.0);

    fmd_model = FractD_KelvinVoigt.Gpa(chirp, [1.0, 1.0, 0.0])
    kelvinvoigt = KelvinVoigt.Gpa(chirp, [1.0, 1.0])

    all(i -> isapprox(fmd_model[i], kelvinvoigt[i]), eachindex(chirp))
end
@test FractD_KelvinVoigt_Gp_reduce()

function Fract_KelvinVoigt_Gp_reduce()
    ω_step = 0.01;
    chirp = collect(0.0:ω_step:10000.0);

    fm_model = Fract_KelvinVoigt.Gppa(chirp, [1.0, 1.0, 1.0, 0.0])
    kelvinvoigt = KelvinVoigt.Gppa(chirp, [1.0, 1.0])

    all(i -> isapprox(fm_model[i], kelvinvoigt[i]), eachindex(chirp))
end
@test Fract_KelvinVoigt_Gp_reduce()

function FractS_KelvinVoigt_Gp_reduce()
    ω_step = 0.01;
    chirp = collect(0.0:ω_step:10000.0);

    fms_model = FractS_KelvinVoigt.Gppa(chirp, [1.0, 1.0, 1.0])
    kelvinvoigt = KelvinVoigt.Gppa(chirp, [1.0, 1.0])

    all(i -> isapprox(fms_model[i], kelvinvoigt[i]), eachindex(chirp))
end
@test FractS_KelvinVoigt_Gp_reduce()

function FractD_KelvinVoigt_Gp_reduce()
    ω_step = 0.01;
    chirp = collect(0.0:ω_step:10000.0);

    fmd_model = FractD_KelvinVoigt.Gppa(chirp, [1.0, 1.0, 0.0])
    kelvinvoigt = KelvinVoigt.Gppa(chirp, [1.0, 1.0])

    all(i -> isapprox(fmd_model[i], kelvinvoigt[i]), eachindex(chirp))
end
@test FractD_KelvinVoigt_Gp_reduce()
