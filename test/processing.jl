# to convert to 'smooth' function
# function _smoothgauss(tol)
#     x = collect(0.0:0.01:30)

#     f = 0.2
#     y = sin.(2π*f*x)

#     ys = RHEOS.smoothgauss(y, x, 1/f; pad="circular")

#     reduced_power = 0.5642720830237535
#     ytest = reduced_power*sin.(2π*f*x)

#     all(i -> (ys[i]-ytest[i])<tol, 1:length(x))
# end
# @test _smoothgauss(tol)