#!/usr/bin/env julia

"""
    contact_hertz(f::Array{Float64,1}, δ::Array{Float64,1}; _param::Float64 = NaN64)

Apply spherical Hertz contact model to find approximate contact point.
"""
function contact_hertz(f::Array{Float64,1}, δ::Array{Float64,1}; _param::Float64 = NaN64)

    # regularise displacement
    δ = -δ - minimum(-δ)

    # radius of sphere
    R = _param
    
    # poisson ratio - assume incompressible
    ν = 0.5

    # get approximate δ0 and YM values
    (init_δ₀, init_YM, tilt_approx, offset_approx) =  hertz_approx(f, δ, R, ν)

    # least squares fit
    opt = Opt(:LN_SBPLX, 4) # 4 parameters to fit
    
    bounds_param = 0.25 # +/- limits allowance from initial contact point guess
    lower_bounds!(opt, [init_δ₀*(1 - bounds_param), 0.0, -5.0, -5.0])
    upper_bounds!(opt, [init_δ₀*(1 + bounds_param), 1e10, 5.0, 5.0])

    xtol_rel!(opt,1e-5)

    min_objective!(opt, (params, grad) -> hertz_cost(f, δ, R, ν, params[1], params[2], params[3], params[4]))

    params_init = [init_δ₀, init_YM, tilt_approx, offset_approx]

    (minf, minx, ret) = optimize(opt, params_init)

    println("min parameters $minx")
    println("min cost: $minf")

    f_fitted = hertz_model(δ, R, ν, minx[1], minx[2], minx[3], minx[4])
    plot(δ, f)
    plot(δ, f_fitted, "--")
    show()
end

"""
    contact_threshold(f::Array{Float64,1}, δ::Array{Float64,1}; _param::Float64 = NaN64)

Applies force threshold to find approximate contact point.
"""
function contact_threshold(f::Array{Float64,1}, δ::Array{Float64,1}; _param::Float64 = NaN64)

    # get index of contact
    cp_index = minimum(find(f .>= _param))

    # return 
    cp_index

end

"""
    contact_none(f::Array{Float64,1}, δ::Array{Float64,1}; _param::Float64 = NaN64)

No contact model applied.
"""
function contact_none(f::Array{Float64,1}, δ::Array{Float64,1}; _param::Float64 = NaN64)

    # return index 1
    1
    
end

###########################################
# Convenience functions for contact_hertz #
###########################################
"""
Generate approximate Hertz fit a, b, E, cp to use as initial conditions in proper fit.
"""
function hertz_approx(f::Array{Float64,1}, δ::Array{Float64,1}, R::Float64, ν::Float64)

    # approximate contact point, assume linear near end of approach section
    # and trace line down to axis
    last_section = round(Int64, 0.98*length(f)):length(f)

    avg_gradient = mean(deriv(f[last_section], δ[last_section]))

    c = f[end] - avg_gradient*δ[end]

    approx_δ₀ = -c/avg_gradient

    # transform to find approx YM
    δᵦ = (heaviside(δ - approx_δ₀).*δ).^(3/2)

    prefactor = (4/3)*sqrt(R)/(1 - ν^2)

    avg_gradient_transformed = mean(deriv(f[last_section], δᵦ[last_section]))

    approx_YM = avg_gradient_transformed/prefactor

    # estimate offset and tilt
    first_section = 1:round(Int64, 0.5*length(f))

    df_dh = deriv(f, δ)

    tilt_approx = mean(df_dh[first_section])

    offset_approx = mean(f[first_section])

    # return
    approx_δ₀, approx_YM, tilt_approx, offset_approx
end

"""
Heaviside function
"""
function heaviside(x)
   0.5*(sign.(x) + 1)
end

"""
Hertz spherical model 
"""
function hertz_model(δ::Array{Float64,1}, R::Float64, ν::Float64, δ₀::Float64, E::Float64, tilt::Float64, offset::Float64)
    # define prefactor for Hertz sphere model
    prefactor = (4/3)*sqrt(R)/(1 - ν^2)

    δ_shift = δ - δ₀

    # force
    F = tilt*δ + offset + prefactor*E*((δ_shift.*heaviside(δ_shift)).^(3/2))
end

"""
Hertz cost function
"""
function hertz_cost(f::Array{Float64,1}, δ::Array{Float64,1}, R::Float64, ν::Float64, δ₀::Float64, E::Float64, tilt::Float64, offset::Float64)
    # difference
    Δ = f - hertz_model(δ, R, ν, δ₀, E, tilt, offset)

    # regularisation parameter for offset and tilt
    λ = 0.1

    # sum of squares
    cost = sum(Δ.^2) + λ*(tilt + offset)^2
end