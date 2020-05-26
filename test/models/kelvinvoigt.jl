println("===============================================")
println("Testing kelvinvoigt.jl")
println("===============================================")

function Fract_KelvinVoigt_G_reduce()
    dt = 0.001
    # fractional modulus is NaN not inf at 0.0 so start at dt
    t = collect(dt:dt:10.0)

    fm_model = relaxmod(Fract_KelvinVoigt, t, [1.0, 1.0, 1.0, 0.0])
    kelvinvoigt = relaxmod(KelvinVoigt, t, [1.0, 1.0])

    all(i -> isapprox(fm_model[i], kelvinvoigt[i]), eachindex(t))
end
@test Fract_KelvinVoigt_G_reduce()

function FractS_KelvinVoigt_G_reduce()
    dt = 0.001
    # fractional modulus is NaN not inf at 0.0 so start at dt
    t = collect(dt:dt:10.0)

    fms_model = relaxmod(FractS_KelvinVoigt, t, [1.0, 1.0, 1.0])
    kelvinvoigt = relaxmod(KelvinVoigt, t, [1.0, 1.0])

    all(i -> isapprox(fms_model[i], kelvinvoigt[i]), eachindex(t))
end
@test FractS_KelvinVoigt_G_reduce()

function FractD_KelvinVoigt_G_reduce()
    dt = 0.001
    t = collect(0.0:dt:10.0)

    fmd_model = relaxmod(FractD_KelvinVoigt, t, [1.0, 1.0, 0.0])
    kelvinvoigt = relaxmod(KelvinVoigt, t, [1.0, 1.0])

    all(i -> isapprox(fmd_model[i], kelvinvoigt[i]), eachindex(t))
end
@test FractD_KelvinVoigt_G_reduce()



function Fract_KelvinVoigt_J_reduce(tol)
    dt = 0.001
    t = collect(0.0:dt:10.0)

    fm_model = creepmod(Fract_KelvinVoigt, t, [1.0, 0.99999, 1.0, 0.00001])
    kelvinvoigt = creepmod(KelvinVoigt, t, [1.0, 1.0])

    all(i -> isapprox(fm_model[i], kelvinvoigt[i], atol=tol), eachindex(t))
end
@test Fract_KelvinVoigt_J_reduce(tol)

function FractS_KelvinVoigt_J_reduce(tol)
    dt = 0.001
    t = collect(0.0:dt:10.0)

    fms_model = creepmod(FractS_KelvinVoigt, t, [1.0, 0.99999, 1.0])
    kelvinvoigt = creepmod(KelvinVoigt, t, [1.0, 1.0])

    all(i -> isapprox(fms_model[i], kelvinvoigt[i], atol=tol), eachindex(t))
end
@test FractS_KelvinVoigt_J_reduce(tol)

function FractD_KelvinVoigt_J_reduce(tol)
    dt = 0.001
    t = collect(0.0:dt:10.0)

    fmd_model = creepmod(FractD_KelvinVoigt, t, [1.0, 1.0, 0.00001])
    kelvinvoigt = creepmod(KelvinVoigt, t, [1.0, 1.0])

    all(i -> isapprox(fmd_model[i], kelvinvoigt[i], atol=tol), eachindex(t))
end
@test FractD_KelvinVoigt_J_reduce(tol)



function Fract_KelvinVoigt_Gp_reduce()
    ω_step = 0.01;
    chirp = collect(0.0:ω_step:10000.0);

    fm_model = storagemod(Fract_KelvinVoigt, chirp, [1.0, 1.0, 1.0, 0.0])
    kelvinvoigt = storagemod(KelvinVoigt, chirp, [1.0, 1.0])

    all(i -> isapprox(fm_model[i], kelvinvoigt[i]), eachindex(chirp))
end
@test Fract_KelvinVoigt_Gp_reduce()

function FractS_KelvinVoigt_Gp_reduce()
    ω_step = 0.01;
    chirp = collect(0.0:ω_step:10000.0);

    fms_model = storagemod(FractS_KelvinVoigt, chirp, [1.0, 1.0, 1.0])
    kelvinvoigt = storagemod(KelvinVoigt, chirp, [1.0, 1.0])

    all(i -> isapprox(fms_model[i], kelvinvoigt[i]), eachindex(chirp))
end
@test FractS_KelvinVoigt_Gp_reduce()

function FractD_KelvinVoigt_Gp_reduce()
    ω_step = 0.01;
    chirp = collect(0.0:ω_step:10000.0);

    fmd_model = storagemod(FractD_KelvinVoigt, chirp, [1.0, 1.0, 0.0])
    kelvinvoigt = storagemod(KelvinVoigt, chirp, [1.0, 1.0])

    all(i -> isapprox(fmd_model[i], kelvinvoigt[i]), eachindex(chirp))
end
@test FractD_KelvinVoigt_Gp_reduce()



function Fract_KelvinVoigt_Gp_reduce()
    ω_step = 0.01;
    chirp = collect(0.0:ω_step:10000.0);

    fm_model = lossmod(Fract_KelvinVoigt, chirp, [1.0, 1.0, 1.0, 0.0])
    kelvinvoigt = lossmod(KelvinVoigt, chirp, [1.0, 1.0])

    all(i -> isapprox(fm_model[i], kelvinvoigt[i]), eachindex(chirp))
end
@test Fract_KelvinVoigt_Gp_reduce()

function FractS_KelvinVoigt_Gp_reduce()
    ω_step = 0.01;
    chirp = collect(0.0:ω_step:10000.0);

    fms_model = lossmod(FractS_KelvinVoigt, chirp, [1.0, 1.0, 1.0])
    kelvinvoigt = lossmod(KelvinVoigt, chirp, [1.0, 1.0])

    all(i -> isapprox(fms_model[i], kelvinvoigt[i]), eachindex(chirp))
end
@test FractS_KelvinVoigt_Gp_reduce()

function FractD_KelvinVoigt_Gp_reduce()
    ω_step = 0.01;
    chirp = collect(0.0:ω_step:10000.0);

    fmd_model = lossmod(FractD_KelvinVoigt, chirp, [1.0, 1.0, 0.0])
    kelvinvoigt = lossmod(KelvinVoigt, chirp, [1.0, 1.0])

    all(i -> isapprox(fmd_model[i], kelvinvoigt[i]), eachindex(chirp))
end
@test FractD_KelvinVoigt_Gp_reduce()
