fl = ILt(s -> 1/s^3, talbot)
@test isapprox( fl(1.0), 0.5; atol =  1e-9)
@test isapprox( fl.([1.0,2.0]), [0.5, 2.0]; atol =  1e-9)

fl = ILt(s -> 1/s^3, gwr,8)
@test isapprox( fl(1.0), 0.5; atol = 1e-5)

p = ILtPair(Weeks(s -> s/(1+s^2)), cos)

e1 = abserr(p, 10.0)
optimize(p,10.0)
e2 = abserr(p, 10.0)
@test e2 < e1

@test typeof(Weeks(iltpair_power(5))) <: ILtPair

# Default method is talbot
fl1 = ILt( s -> 1/s)
@test fl1.iltfunc == talbot
