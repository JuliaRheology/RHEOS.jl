tol = 1e-2

function _trapz(tol)
    x = [i for i in -5.0:0.01:5.0]
    y = 3*x.^2

    (RHEOS.trapz(y, x) - 250.0)<tol
end
@test _trapz(tol)

function _derivCD(tol)
    x = collect(0.0:0.001:1.5)
    y = x.^2
    dy = 2*x
    dy_numeric = RHEOS.derivCD(y, x)

    all(i -> (dy[i]-dy_numeric[i])<tol, 1:length(dy_numeric))
end
@test _derivCD(tol)

function _derivBD(tol)
    x = collect(0.0:0.001:1.5)
    y = x.^2
    dy = 2*x
    dy_numeric = RHEOS.derivBD(y, x)

    all(i -> (dy[i]-dy_numeric[i])<tol, 1:length(dy_numeric))
end
@test _derivBD(tol)

@test RHEOS.quasinull([-1.0])==true
@test RHEOS.quasinull([1.0])==false

@test RHEOS.constantcheck([1.0, 2.0, 3.0])==true
@test RHEOS.constantcheck([1.0, 2.0, 4.0])==false

@test RHEOS.sampleratecompare([1.0, 2.0, 3.0], [2.0, 3.0, 4.0])==true
@test RHEOS.sampleratecompare([1.0, 2.0, 3.0], [0.1, 0.2, 0.3])==false

# @test _var_resample()

@test RHEOS.closestindex([1.0, 2.0, 3.0], 1.7)==2

@test RHEOS.closestindices([1.0, 2.0, 3.0], [1.7, 3.4])==[2, 3]

@test RHEOS.mapback([1.0, 1.23, 1.24, 1.5], [1.0, 1.1, 1.2, 1.3, 1.4, 1.5])==[1, 3, 6]
@test RHEOS.mapback([1.0, 1.23, 1.24, 1.5], [1.0, 1.1, 1.2, 1.3, 1.4, 1.5]; return_indices=false)==[1.0, 1.2, 1.5]

# function _downsample()
# end
# @test RHEOS.
