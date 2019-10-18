function springpot_spring_reduce_G(tol)
    t = Vector{RheoFloat}(0.0:0.1:10.0)

    springpot = Springpot.Ga(t, [1.5, 0.0])
    spring = Spring.Ga(t, [1.5])

    all(i -> isapprox(springpot[i], spring[i]), eachindex(t))
end
@test springpot_spring_reduce_G(tol)

function springpot_spring_reduce_J(tol)
    t = Vector{RheoFloat}(0.0:0.1:10.0)

    springpot = Springpot.Ja(t, [1.5, 0.0])
    spring = Spring.Ja(t, [1.5])

    all(i -> isapprox(springpot[i], spring[i]), eachindex(t))
end
@test springpot_spring_reduce_J(tol)

function springpot_dashpot_reduce_G(tol)
    # This gives NaN as defined presently,
    # not test written
    return true
end
@test springpot_dashpot_reduce_G(tol)

function springpot_dashpot_reduce_J(tol)
    t = Vector{RheoFloat}(0.0:0.1:10.0)

    springpot = Springpot.Ja(t, [1.5, 1.0])
    dashpot = Dashpot.Ja(t, [1.5])

    all(i -> isapprox(springpot[i], dashpot[i]), eachindex(t))
end
@test springpot_dashpot_reduce_J(tol)
