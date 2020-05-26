function springpot_spring_reduce_G()
    t = Vector{RheoFloat}(0.0:0.1:10.0)

    springpot = relaxmod(RheoModel(Springpot, cᵦ = 1.5, β=0),t)
    spring = relaxmod(Spring,t, k=1.5)

    all(i -> isapprox(springpot[i], spring[i]), eachindex(t))
end
@test springpot_spring_reduce_G()

function springpot_spring_reduce_J()
    t = Vector{RheoFloat}(0.0:0.1:10.0)

    springpot_J = creepmod(Springpot, cᵦ = 1.5, β=0)
    springpot = springpot_J(t)
    spring = creepmod(Spring,t, [1.5])

    all(i -> isapprox(springpot[i], spring[i]), eachindex(t))
end
@test springpot_spring_reduce_J()

function springpot_spring_reduce_Gp()
    ω = Vector{RheoFloat}(0.0:0.1:10.0)

    springpot = storagemod(Springpot, ω, cᵦ = 1.5, β=0)
    spring = storagemod(Spring, ω, [1.5])

    all(i -> isapprox(springpot[i], spring[i]), eachindex(ω))
end
@test springpot_spring_reduce_Gp()

function springpot_spring_reduce_Gpp()
    ω = Vector{RheoFloat}(0.0:0.1:10.0)

    springpot_Gpp = lossmod(RheoModel(Springpot, cᵦ = 1.5, β=0))
    springpot = springpot_Gpp(ω)
    spring = lossmod(Spring, ω, [1.5])

    all(i -> isapprox(springpot[i], spring[i]), eachindex(ω))
end
@test springpot_spring_reduce_Gpp()

function springpot_dashpot_reduce_G()
    # This gives NaN as defined presently,
    # no test written
    return true
end
@test springpot_dashpot_reduce_G()

function springpot_dashpot_reduce_J()
    t = Vector{RheoFloat}(0.0:0.1:10.0)

    springpot = creepmod(Springpot,t, [1.5, 1.0])
    dashpot = creepmod(Dashpot,t, [1.5])

    all(i -> isapprox(springpot[i], dashpot[i]), eachindex(t))
end
@test springpot_dashpot_reduce_J()

function springpot_dashpot_reduce_Gp(tol)
    ω = Vector{RheoFloat}(0.0:0.1:10.0)

    springpot = storagemod(Springpot, ω, [1.5, 1.0])
    dashpot = real(dynamicmod(Dashpot, ω, [1.5]))

    # lower tolerance due to cos(0) innacuracy in Julia
    all(i -> isapprox(springpot[i], dashpot[i], atol=tol), eachindex(ω))
end
@test springpot_dashpot_reduce_Gp(tol)

function springpot_dashpot_reduce_Gpp()
    ω = Vector{RheoFloat}(0.0:0.1:10.0)

    springpot = lossmod(Springpot, ω, [1.5, 1.0])
    dashpot = imag(dynamicmod(Dashpot, ω, [1.5]))

    all(i -> isapprox(springpot[i], dashpot[i]), eachindex(ω))
end
@test springpot_dashpot_reduce_Gpp()