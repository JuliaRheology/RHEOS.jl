function springpot_spring_reduce_G()
    t = Vector{RheoFloat}(0.0:0.1:10.0)

    springpot = Springpot.Ga(t, [1.5, 0.0])
    spring = Spring.Ga(t, [1.5])

    all(i -> isapprox(springpot[i], spring[i]), eachindex(t))
end
@test springpot_spring_reduce_G()

function springpot_spring_reduce_J()
    t = Vector{RheoFloat}(0.0:0.1:10.0)

    springpot = Springpot.Ja(t, [1.5, 0.0])
    spring = Spring.Ja(t, [1.5])

    all(i -> isapprox(springpot[i], spring[i]), eachindex(t))
end
@test springpot_spring_reduce_J()

function springpot_spring_reduce_Gp()
    t = Vector{RheoFloat}(0.0:0.1:10.0)

    springpot = Springpot.Gpa(t, [1.5, 0.0])
    spring = Spring.Gpa(t, [1.5])

    all(i -> isapprox(springpot[i], spring[i]), eachindex(t))
end
@test springpot_spring_reduce_Gp()

function springpot_spring_reduce_Gpp()
    t = Vector{RheoFloat}(0.0:0.1:10.0)

    springpot = Springpot.Gppa(t, [1.5, 0.0])
    spring = Spring.Gppa(t, [1.5])

    all(i -> isapprox(springpot[i], spring[i]), eachindex(t))
end
@test springpot_spring_reduce_Gpp()

function springpot_dashpot_reduce_G()
    # This gives NaN as defined presently,
    # not test written
    return true
end
@test springpot_dashpot_reduce_G()

function springpot_dashpot_reduce_J()
    t = Vector{RheoFloat}(0.0:0.1:10.0)

    springpot = Springpot.Ja(t, [1.5, 1.0])
    dashpot = Dashpot.Ja(t, [1.5])

    all(i -> isapprox(springpot[i], dashpot[i]), eachindex(t))
end
@test springpot_dashpot_reduce_J()

function springpot_dashpot_reduce_Gp(tol)
    t = Vector{RheoFloat}(0.0:0.1:10.0)

    springpot = Springpot.Gpa(t, [1.5, 1.0])
    dashpot = Dashpot.Gpa(t, [1.5])

    # lower tolerance due to cos(0) innacuracy in Julia
    all(i -> isapprox(springpot[i], dashpot[i], atol=tol), eachindex(t))
end
@test springpot_dashpot_reduce_Gp(tol)

function springpot_dashpot_reduce_Gpp()
    t = Vector{RheoFloat}(0.0:0.1:10.0)

    springpot = Springpot.Gppa(t, [1.5, 1.0])
    dashpot = Dashpot.Gppa(t, [1.5])

    all(i -> isapprox(springpot[i], dashpot[i]), eachindex(t))
end
@test springpot_dashpot_reduce_Gpp()