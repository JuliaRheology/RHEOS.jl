println("===============================================")
println("Testing interface.jl")
println("===============================================")


function _interface_f()
    AFMitf = AFM(4.0)
    ramp = RheoTimeData(AFMitf, t = Vector(0:0.1:10), d = Vector(0:0.05:5))
    model=RheoModel(Spring, k=2)
    data = modelpredict(ramp,model)

    # check expected force for d = 1.0
    isapprox(data[AFMitf].f[21], 16/3., atol=tol)
end
@test _interface_f()


function _interface_d()
    AFMitf = AFM(2)
    ramp = RheoTimeData(AFMitf, t = Vector(0:0.1:10), f = Vector(0:0.05:5))
    model=RheoModel(Spring, k=0.25)
    data = modelpredict(ramp,model)

    # check expected displacement for f = 1.0
    isapprox(data[AFMitf].d[21], 2*(0.75^(2/3)), atol=tol)
end
@test _interface_d()
