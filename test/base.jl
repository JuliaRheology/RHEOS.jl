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
    dy_numeric = derivCD(y, x)

    all(i -> (dy[i]-dy_numeric[i])<tol, 1:length(dy_numeric))
end
@test _derivCD(tol)

function _derivBD(tol)
    x = collect(0.0:0.001:1.5)
    y = x.^2
    dy = 2*x
    dy_numeric = derivBD(y, x)

    all(i -> (dy[i]-dy_numeric[i])<tol, 1:length(dy_numeric))
end
@test _derivBD(tol)

@test RHEOS.quasinull([-1.0])==true
@test RHEOS.quasinull([1.0])==false

@test RHEOS.constantcheck([1.0, 2.0, 3.0])==true
@test RHEOS.constantcheck([1.0, 2.0, 4.0])==false

@test RHEOS.sampleratecompare([1.0, 2.0, 3.0], [2.0, 3.0, 4.0])==true
@test RHEOS.sampleratecompare([1.0, 2.0, 3.0], [0.1, 0.2, 0.3])==false

# How to test smoothing works as expected?
# x = collect(0.0:0.01:30)
# f = 0.2
# ys = RHEOS.smoothgauss(y, x, 1/f)
# y = sin.(2π*f*x)
# ytest = maximum(ys)*sin.(2π*f*x)
# plot(x, y)
# plot(x, ys, "--")
# plot(x, ytest, "-.")
# function _smoothgauss



