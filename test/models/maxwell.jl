println("===============================================")
println("Testing maxwell.jl")
println("===============================================")

function Fract_Maxwell_G_reduce()
    dt = 0.001
    t = collect(0.0:dt:10.0)

    fm_model = relaxmod(Fract_Maxwell, t, [1.0, 1.0, 1.0, 0.0])
    maxwell = relaxmod(Maxwell, t, [1.0, 1.0])

    all(i -> isapprox(fm_model[i], maxwell[i]), eachindex(t))
end
@test Fract_Maxwell_G_reduce()

function FractS_Maxwell_G_reduce()
    dt = 0.001
    t = collect(0.0:dt:10.0)

    fms_model = relaxmod(FractS_Maxwell, t, [1.0, 1.0, 1.0])
    maxwell = relaxmod(Maxwell, t, [1.0, 1.0])

    all(i -> isapprox(fms_model[i], maxwell[i]), eachindex(t))
end
@test FractS_Maxwell_G_reduce()

function FractD_Maxwell_G_reduce()
    dt = 0.001
    t = collect(0.0:dt:10.0)

    fmd_model = relaxmod(FractD_Maxwell, t, [1.0, 1.0, 0.0])
    maxwell = relaxmod(Maxwell, t, [1.0, 1.0])

    all(i -> isapprox(fmd_model[i], maxwell[i]), eachindex(t))
end
@test FractD_Maxwell_G_reduce()



function Fract_Maxwell_J_reduce()
    dt = 0.001
    t = collect(0.0:dt:10.0)

    fm_model = creepmod(Fract_Maxwell, t, [1.0, 1.0, 1.0, 0.0])
    maxwell = creepmod(Maxwell, t, [1.0, 1.0])

    all(i -> isapprox(fm_model[i], maxwell[i]), eachindex(t))
end
@test Fract_Maxwell_J_reduce()

function FractS_Maxwell_J_reduce()
    dt = 0.001
    t = collect(0.0:dt:10.0)

    fms_model = creepmod(FractS_Maxwell, t, [1.0, 1.0, 1.0])
    maxwell = creepmod(Maxwell, t, [1.0, 1.0])

    all(i -> isapprox(fms_model[i], maxwell[i]), eachindex(t))
end
@test FractS_Maxwell_J_reduce()

function FractD_Maxwell_J_reduce()
    dt = 0.001
    t = collect(0.0:dt:10.0)

    fmd_model = creepmod(FractD_Maxwell, t, [1.0, 1.0, 0.0])
    maxwell = creepmod(Maxwell, t, [1.0, 1.0])

    all(i -> isapprox(fmd_model[i], maxwell[i]), eachindex(t))
end
@test FractD_Maxwell_J_reduce()



function Fract_Maxwell_Gp_reduce()
    ω_step = 0.01;
    chirp = collect(0.0:ω_step:10000.0);

    fm_model = storagemod(Fract_Maxwell, chirp, [1.0, 1.0, 1.0, 0.0])
    maxwell = storagemod(Maxwell, chirp, [1.0, 1.0])

    all(i -> isapprox(fm_model[i], maxwell[i]), eachindex(chirp))
end
@test Fract_Maxwell_Gp_reduce()

function FractS_Maxwell_Gp_reduce()
    ω_step = 0.01;
    chirp = collect(0.0:ω_step:10000.0);

    fms_model = storagemod(FractS_Maxwell, chirp, [1.0, 1.0, 1.0])
    maxwell = storagemod(Maxwell, chirp, [1.0, 1.0])

    all(i -> isapprox(fms_model[i], maxwell[i]), eachindex(chirp))
end
@test FractS_Maxwell_Gp_reduce()

function FractD_Maxwell_Gp_reduce()
    ω_step = 0.01;
    chirp = collect(0.0:ω_step:10000.0);

    fmd_model = storagemod(FractD_Maxwell, chirp, [1.0, 1.0, 0.0])
    maxwell = storagemod(Maxwell, chirp, [1.0, 1.0])

    all(i -> isapprox(fmd_model[i], maxwell[i]), eachindex(chirp))
end
@test FractD_Maxwell_Gp_reduce()



function Fract_Maxwell_Gpp_reduce()
    ω_step = 0.01;
    chirp = collect(0.0:ω_step:10000.0);

    fm_model = lossmod(Fract_Maxwell, chirp, [1.0, 1.0, 1.0, 0.0])
    maxwell = lossmod(Maxwell, chirp, [1.0, 1.0])

    all(i -> isapprox(fm_model[i], maxwell[i]), eachindex(chirp))
end
@test Fract_Maxwell_Gpp_reduce()

function FractS_Maxwell_Gpp_reduce()
    ω_step = 0.01;
    chirp = collect(0.0:ω_step:10000.0);

    fms_model = lossmod(FractS_Maxwell, chirp, [1.0, 1.0, 1.0])
    maxwell = lossmod(Maxwell, chirp, [1.0, 1.0])

    all(i -> isapprox(fms_model[i], maxwell[i]), eachindex(chirp))
end
@test FractS_Maxwell_Gpp_reduce()

function FractD_Maxwell_Gpp_reduce()
    ω_step = 0.01;
    chirp = collect(0.0:ω_step:10000.0);

    fmd_model = lossmod( FractD_Maxwell, chirp, [1.0, 1.0, 0.0])
    maxwell = lossmod(Maxwell, chirp, [1.0, 1.0])

    all(i -> isapprox(fmd_model[i], maxwell[i]), eachindex(chirp))
end
@test FractD_Maxwell_Gpp_reduce()
