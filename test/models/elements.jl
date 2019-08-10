function springpot_spring_reduce_G(tol)
    t = Vector{RheoFloat}(0.0:0.1:10.0)

    springpot = Springpot.Ga(t, [1.5, 0.0])
    spring = Spring.Ga(t, [1.5])

    all(i -> isapprox(springpot[i], spring[i]), eachindex(t))
end
@test springpot_spring_reduce_G(tol)