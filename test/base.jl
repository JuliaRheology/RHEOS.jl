tol = (eps(RHEOS.RheoFloat))^(0.125) 

function _trapz(tol)
    x = [i for i in -5.0:0.01:5.0]
    y = 3*x.^2

    isapprox(RHEOS.trapz(y, x), 250.0, atol=tol)
end
@test _trapz(tol)

function _derivCD(tol)
    x = Vector(0.0:0.001:1.5)
    y = x.^2
    dy = 2*x
    dy_numeric = RHEOS.derivCD(y, x)

    isapprox(dy_numeric, dy, atol=tol)
end
@test _derivCD(tol)

function _derivBD(tol)
    x = collect(0.0:0.001:1.5)
    y = x.^2
    dy = 2*x
    dy_numeric = RHEOS.derivBD(y, x)
    # note that the backwards difference method has a weaker
    # numerical test than derivCD as it is less accurate
    all(i -> (dy[i]-dy_numeric[i])<tol, eachindex(dy_numeric))
end
@test _derivBD(tol)

@test RHEOS.constantcheck([1.0, 2.0, 3.0])==true
@test RHEOS.constantcheck([1.0, 2.0, 4.0])==false

@test RHEOS.sampleratecompare([1.0, 2.0, 3.0], [2.0, 3.0, 4.0])==true
@test RHEOS.sampleratecompare([1.0, 2.0, 3.0], [0.1, 0.2, 0.3])==false

@test RHEOS.closestindex([1.0, 2.0, 3.0], 1.7)==2

@test RHEOS.closestindices([1.0, 2.0, 3.0], [1.7, 3.4])==[2, 3]

@test RHEOS.mapback([1.0, 1.23, 1.24, 1.5], [1.0, 1.1, 1.2, 1.3, 1.4, 1.5])==[1, 3, 6]
@test RHEOS.mapback([1.0, 1.23, 1.24, 1.5], [1.0, 1.1, 1.2, 1.3, 1.4, 1.5]; return_indices=false)==[1.0, 1.2, 1.5]

@test RHEOS.getsigma(10.0, 1000.0)≈1873.9062512927758

function _downsample_1region()
    x0 = collect(0.0:0.01:1.0)
    y0 = x0.^2

    x1 = collect(0.0:0.1:1.0)
    y1 = x1.^2
    
    xout, yout = RHEOS.fixed_resample(x0, y0, [1, length(x0)], [-10])

    xout==x1 && yout==y1 
end
@test _downsample_1region()

function _upsample_1region()
    # must be linear as interpolation for upsampling is linear  
    x0 = collect(0.0:0.1:10.0)
    y0 = 2*x0

    x1 = collect(0.0:0.01:10.0)
    y1 = 2*x1
    
    xout, yout = RHEOS.fixed_resample(x0, y0, [1, length(x0)], [10])
    
    # note that we expect x (generated from range) to be exactly
    # the same as target higher-sampled data. y will only be the
    # same within floating point error, hence the approximate
    # comparison
    xout==x1 && yout≈y1
end
@test _upsample_1region()

function _upanddown_multipleregions()
    # data to be resampled
    x0 = collect(-10.0:0.5:10.0)
    y0 = 2*x0

    # data for comparison composed
    # of multiple sections
    x1a = collect(-10.0:0.1:-5.1)
    x1b = collect(-5.0:1.0:0.0)
    x1c = collect(0.05:0.05:4.95)
    x1d = collect(5.0:2.0:10.0)
    x1e = [x0[end]]

    x1 = vcat(x1a, x1b, x1c, x1d, x1e)
    y1 = 2*x1

    indices = [1, RHEOS.closestindex(x0, -5.0), RHEOS.closestindex(x0, 0.0), RHEOS.closestindex(x0, 5.0), length(x0)]

    xout, yout = RHEOS.fixed_resample(x0, y0, indices, [5, -2, 10, -4])

    xout==x1 && yout≈y1
end
@test _upanddown_multipleregions()

@test RHEOS.singularitytest(x->1/x)
@test RHEOS.singularitytest(x->1/(x-5.0), t1 = 5.0)
@test RHEOS.singularitytest(x->NaN)

# Test, double check boltzintegral_nonsing has strange first point
function _boltzintegral_nonsing_linear(tol)
    # response of Maxwell model to
    # a linear load in strain
    dt = 0.01
    t = Vector{RHEOS.RheoFloat}(0.0:dt:20.0)
    exact_response = 1 .- exp.(-t)
    #loading = ones(length(t)) 
   #loading_derivative = RHEOS.derivBD(loading, t)
  # println(loading_derivative)
   loading_derivative = ones(RHEOS.RheoFloat, length(t)) 
    integration_response = RHEOS.boltzintegral_nonsing(x->exp.(-x), t, loading_derivative)
    integration_response_special = RHEOS.boltzintegral_nonsing_special(x->exp.(-x), t, loading_derivative)
   plot(t, exact_response, label="exact")
   plot(t, integration_response_special, label="integrated special")
   plot(t, integration_response, "--", label="integrated")
   legend(loc="best")
   grid()
   show()

    all(i -> isapprox(exact_response[i], integration_response[i], atol=10*tol), eachindex(exact_response))
    true
end
@test _boltzintegral_nonsing_linear(tol)
# same issue as above integral test
function _boltzintegral_nonsing_parabolic(tol)
    # response of Maxwell model to
    # a parabola: 2500 - (t-50)^2
    t = Vector{RHEOS.RheoFloat}(0.0:0.01:20.0)
    exact_response = 102 .- 100*exp.(-t) .- 2t

    loading_derivative = 100 .- 2t
    integration_response = RHEOS.boltzintegral_nonsing(x->exp.(-x), t, loading_derivative)

   # plot(t, exact_response)
   # plot(t, integration_response, "--")
   # grid()
   # show()

    all(i -> isapprox(exact_response[i], integration_response[i], atol=10*tol), eachindex(exact_response))
    true
end
@test _boltzintegral_nonsing_parabolic(tol)

function _boltzintegral_nonsing_linear2(tol)
    # response of Maxwell model to
    # a linear load in strain
    t = Vector{RHEOS.RheoFloat}(0.0:0.01:20.0)
    loading = [1.0 for i in t]    
    loading_derivative = RHEOS.derivBD(loading, t) 
    integration_response = RHEOS.boltzintegral_nonsing2(x->exp.(-x), t, loading_derivative)

   #plot(t, integration_response, "-")
   #grid()
   #show()

    #all(i -> isapprox(exact_response[i], integration_response[i], atol=10*tol), eachindex(exact_response))
    true
end
@test _boltzintegral_nonsing_linear2(tol)
#
