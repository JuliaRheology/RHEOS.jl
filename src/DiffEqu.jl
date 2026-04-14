using NumFracDiff

#-----------------------------------------------------
function generate_GL_weights(order, length, weights)
    
    weights[1] = 1.0

    alpha_plus_one = order + 1.0
    
    for k in 2:length
        weights[k] = weights[k-1] * (1.0 - alpha_plus_one/(k - 1.0)) 
    end
end


struct SupportVectorsExt
    rhs::Vector{Float64}
    lhs::Vector{Float64}
    weights::Vector{Float64}
    deriv::Vector{Float64}
end

using FFTW
using DSP

"""
    compute_GL_frac_deriv_fftfilt(data, weights, order, dt)
Compute the vector of the fractional derivatives using the FFT-based filtering method and the GL weigths.
This method is generally faster than the FDM, since it has a complexity of O(N log N) instead of O(N^2), but works only on serial.
"""
function compute_GL_frac_deriv_fftfilt(data, weights, order, dt)
    
    deriv_res = fftfilt(weights, data)
    
    inv_dt_pow = 1.0 / (dt^order)
    deriv_res .*= inv_dt_pow
    
    return deriv_res
end
#-------------------------------------------------------
function modelfit(data::RheoTimeData, 
    model::RheoModelClass,
    modloading::LoadingType,
    fittype::Differential = Differential();
    method=RL(),
    p0::Union{NamedTuple,Nothing,Dict} = nothing,
    lo::Union{NamedTuple,Nothing,Dict} = nothing,
    hi::Union{NamedTuple,Nothing,Dict} = nothing,
    verbose::Bool = false,
    rel_tol_f::Union{Real,Nothing} = nothing,
    rel_tol_x::Union{Real,Nothing} = isnothing(rel_tol_f) ? 1e-4 : nothing,
    diff_method="BD",
    weights::Union{Nothing,Vector{T}} = nothing,
    optmethod::Union{Symbol,String}= :LN_SBPLX, 
    opttimeout::Union{Real,Nothing} = nothing,
    optmaxeval::Union{Integer,Nothing} = nothing,
    allowconstraints=true) where T <: Integer

p0a = fill_init_params(model, symbol_to_unicode(p0))
loa = fill_lower_bounds(model, symbol_to_unicode(lo))
hia = fill_upper_bounds(model, symbol_to_unicode(hi))

check = rheotimedatatype(data)
@assert (check == strain_and_stress) "Both stress and strain are required"

# check provided weights are all valid
if !isnothing(weights)
@assert isempty(weights[weights.<1]) "Invalid weighting indices provided"
end

# use correct method for derivative
if diff_method=="BD"
deriv = derivBD
elseif diff_method=="CD"
deriv = derivCD
end

strain_deriv = deriv(data.ϵ, data.t)
stress_deriv = deriv(data.σ, data.t)

equation = model.C

# get time step (only needed for convolution, which requires constant dt so t[2]-t[1] is sufficient)
dt = data.t[2] - data.t[1]

# time must start at 0 for convolution to work properly!
t_zeroed = data.t .- minimum(data.t)

# fit
is_constant = constantcheck(data.t)

# indices weighting only used for constant sample-rate data
if !is_constant && !isnothing(weights)
@warn "Indices weighting not used as variable sample-rate data has been provided"
end

# Perform fitting
# TODO pass the equation instead of the modulus -> probably the entire method needs to be completely rewritten
(minf, minx, ret), timetaken, bytes, gctime, memalloc =
@timed leastsquares_init(   p0a,
                            loa,
                            hia,
                            equation,
                            modloading,
                            t_zeroed,
                            dt,
                            data.ϵ,
                            strain_deriv,
                            data.σ,
                            stress_deriv,
                            model._constraint,
                            fittype;
                            method= method,
                            insight = verbose,
                            constant_sampling = is_constant,
                            singularity = false,
                            rel_tol_x = rel_tol_x,
                            rel_tol_f = rel_tol_f,
                            indweights = weights,
                            optmethod = Symbol(optmethod),
                            opttimeout = opttimeout,
                            optmaxeval = optmaxeval,
                            allowconstraints=allowconstraints)

println("Time: $timetaken s, Why: $ret, Parameters: $minx, Error: $minf")

nt = NamedTuple{Tuple(model.freeparams)}(minx)

if data.log !== nothing
# Preparation of data for log item
info=(comment="Fiting rheological model to data", model_name=model.name, model_params=nt, time_taken=timetaken, stop_reason=ret, error=minf)
params=(model=model, modloading=modloading)
keywords=(p0=p0, lo=lo, hi=hi, rel_tol_x=rel_tol_x, diff_method=diff_method)
# Add data to the log
push!(data.log, RheoLogItem( (type=:analysis, funct=:modelfit, params=params, keywords=keywords), info))
end

return RheoModel(model, nt);
end 



function leastsquares_init(params_init::Vector{RheoFloat},
    low_bounds::RheovecOrNone,
    hi_bounds::RheovecOrNone, 
    equation::DiffEqu{DiffScaFree},
    modloading::LoadingType,
    time_series::Vector{RheoFloat},
    dt::RheoFloat,
    strain::Vector{RheoFloat},
    strain_deriv::Vector{RheoFloat},
    stress::Vector{RheoFloat},
    stress_deriv::Vector{RheoFloat},
    constraint::Union{Vector{FWConstraint},Nothing},
    fittype::Differential;
    method=RL(),
    insight::Bool = false,
    constant_sampling::Bool=true,
    singularity::Bool = false,
    rel_tol_x::Union{Real,Nothing} = nothing,
    rel_tol_f::Union{Real,Nothing} = nothing,
    indweights = nothing,
    optmethod::Symbol = :LN_SBPLX,
    opttimeout::Union{Real,Nothing} = nothing,
    optmaxeval::Union{Integer,Nothing} = nothing,
    allowconstraints= true)


    # initialise NLOpt.Opt object with :LN_SBPLX Subplex algorithm or COBYLA if constrained
    if !(optmethod in [:LN_AUGLAG, :LN_COBYLA]) && constraint ≠ nothing && allowconstraints
        optmethod = :LN_COBYLA
    end
    opt = Opt(optmethod, length(params_init))
    # opt = Opt(:LN_BOBYQA, length(params_init))    # Passing tests
    # opt = Opt(:LN_COBYLA, length(params_init))    # Failing test - not precise enough?

    # set optimiser stopping criteria

    # wall clock timeout
    if !isnothing(opttimeout)
        opttimeout = convert(Float64, opttimeout)
        maxtime!(opt, opttimeout)
    end

    # evaluation cycle ceiling
    if !isnothing(optmaxeval)
        maxeval!(opt, optmaxeval)
    end

    # input parameter change tolerance
    if !isnothing(rel_tol_x)
        rel_tol_x = convert(Float64, rel_tol_x)
        xtol_rel!(opt, rel_tol_x)
    end

    # objective function change tolerance 
    if !isnothing(rel_tol_f)
        rel_tol_f = convert(Float64, rel_tol_f)
        ftol_rel!(opt, rel_tol_f)
    end

    # set lower bounds and upper bounds unless they take null value
    if !isnothing(low_bounds)
        low_bounds = convert(Vector{Float64},low_bounds)
        lower_bounds!(opt, low_bounds)
    end

    if !isnothing(hi_bounds)
        hi_bounds = convert(Vector{Float64}, hi_bounds)
        upper_bounds!(opt, hi_bounds)
    end

    # Convert to float64 to avoid conversion by NLOpt
    params_init = convert(Vector{Float64},params_init)
    time_series = convert(Vector{Float64},time_series)
    strain_deriv = convert(Vector{Float64},strain_deriv)
    stress_deriv = convert(Vector{Float64},stress_deriv)
    dt = convert(Float64, dt)


    data_struct = SupportVectors(
        zeros(length(time_series)),
        zeros(length(time_series)),
    )

    prob = NumDiffProblem(dt=dt,order=0.5,n=length(time_series),method=method)
    ws = init_workspace(prob)
    min_objective!(opt, (params, grad) -> obj_const(params, equation,
                                time_series, dt, 
                                strain, strain_deriv, 
                                stress, stress_deriv,data_struct, prob, ws;
                                _insight = insight))


    if constraint ≠ nothing && allowconstraints
        for c in eachindex(constraint)
            inequality_constraint!(opt,nlopt_constraint_wrapper(constraint[c]), 1e-8)
        end
    end


    # minimise objective func, minx are the parameters resulting in minimum
    (minf, minx, ret) = NLopt.optimize(opt, params_init)

    # return all
    return (convert(RheoFloat,minf), convert(Vector{RheoFloat},minx), ret)

    end

function leastsquares_init_fft(params_init::Vector{RheoFloat},
        low_bounds::RheovecOrNone,
        hi_bounds::RheovecOrNone, 
        equation::DiffEqu{DiffScaFree},
        modloading::LoadingType,
        time_series::Vector{RheoFloat},
        dt::RheoFloat,
        strain::Vector{RheoFloat},
        strain_deriv::Vector{RheoFloat},
        stress::Vector{RheoFloat},
        stress_deriv::Vector{RheoFloat},
        constraint::Union{Vector{FWConstraint},Nothing},
        fittype::Differential;
        method::DiffMethod=GL(),
        insight::Bool = false,
        constant_sampling::Bool=true,
        singularity::Bool = false,
        rel_tol_x::Union{Real,Nothing} = nothing,
        rel_tol_f::Union{Real,Nothing} = nothing,
        indweights = nothing,
        optmethod::Symbol = :LN_SBPLX,
        opttimeout::Union{Real,Nothing} = nothing,
        optmaxeval::Union{Integer,Nothing} = nothing,
        allowconstraints::Bool=true)


    # initialise NLOpt.Opt object with :LN_SBPLX Subplex algorithm
    if !(optmethod in [:LN_AUGLAG, :LN_COBYLA]) && constraint ≠ nothing && allowconstraints
    optmethod = :LN_COBYLA
    end

    opt = Opt(optmethod, length(params_init))
    # opt = Opt(:LN_BOBYQA, length(params_init))    # Passing tests
    # opt = Opt(:LN_COBYLA, length(params_init))    # Failing test - not precise enough?

    # set optimiser stopping criteria

    # wall clock timeout
    if !isnothing(opttimeout)
    opttimeout = convert(Float64, opttimeout)
    maxtime!(opt, opttimeout)
    end

    # evaluation cycle ceiling
    if !isnothing(optmaxeval)
    maxeval!(opt, optmaxeval)
    end

    # input parameter change tolerance
    if !isnothing(rel_tol_x)
    rel_tol_x = convert(Float64, rel_tol_x)
    xtol_rel!(opt, rel_tol_x)
    end

    # objective function change tolerance 
    if !isnothing(rel_tol_f)
    rel_tol_f = convert(Float64, rel_tol_f)
    ftol_rel!(opt, rel_tol_f)
    end

    # set lower bounds and upper bounds unless they take null value
    if !isnothing(low_bounds)
    low_bounds = convert(Vector{Float64},low_bounds)
    lower_bounds!(opt, low_bounds)
    end

    if !isnothing(hi_bounds)
    hi_bounds = convert(Vector{Float64}, hi_bounds)
    upper_bounds!(opt, hi_bounds)
    end

    # Convert to float64 to avoid conversion by NLOpt
    params_init = convert(Vector{Float64},params_init)
    time_series = convert(Vector{Float64},time_series)
    strain_deriv = convert(Vector{Float64},strain_deriv)
    stress_deriv = convert(Vector{Float64},stress_deriv)
    dt = convert(Float64, dt)


    data_struct = SupportVectors(
    zeros(length(time_series)),
    zeros(length(time_series))
    )
    prob = NumDiffProblem(dt=dt,order=0.5,n=length(time_series),method=method)
    ws = init_workspace(prob)

    min_objective!(opt, (params, grad) -> obj_const(params, equation,
                                    time_series, dt, 
                                    strain, strain_deriv, 
                                    stress, stress_deriv,data_struct,prob,ws;
                                    _insight = insight))


    if constraint ≠ nothing && allowconstraints
    for c in eachindex(constraint)
    inequality_constraint!(opt,nlopt_constraint_wrapper(constraint[c]), 1e-8)
    end
    end

    # minimise objective func, minx are the parameters resulting in minimum
    (minf, minx, ret) = NLopt.optimize(opt, params_init)

    return (convert(RheoFloat,minf), convert(Vector{RheoFloat},minx), ret)

end


#=
--------------------------------
Functions used to compute the cost function for parameter fitting
--------------------------------
=#

"""
    obj_const(params, equation, time_series, dt, strain, strain_deriv, stress, stress_deriv,data_struct; SM=true, _insight::Bool = false)

Compute the cost function to minimize during parameter fitting, using the GL method with Short-Memory.
"""
function obj_const(params, equation, time_series, dt, strain, strain_deriv, stress, stress_deriv,data_struct, prob, ws;_insight::Bool = false)
    numerical_coeffs = get_coeffs(equation, params)
    @. data_struct.rhs = 0
    @. data_struct.lhs = 0
    
    L=length(strain)
    # Compute the rhs of the equation
    for term in numerical_coeffs.strain
        order = term[1]
        coeff = term[2]
        if order == 0.0
            data_struct.rhs .+= coeff * strain
        elseif order == 1.0
            data_struct.rhs .+= coeff * strain_deriv
        else
            update_order!(prob,ws,order)
            compute!(prob.method,ws,strain,prob)
            @. data_struct.rhs += coeff * ws.deriv
        end
    end
    
    for term in numerical_coeffs.stress
        order = term[1]
        coeff = term[2]
        if order == 0.0
            data_struct.lhs .+= coeff * stress
        elseif order == 1.0
            data_struct.lhs .+= coeff * stress_deriv
        else
            update_order!(prob,ws,order)
            compute!(prob.method,ws,stress,prob)
            @. data_struct.lhs += coeff * ws.deriv
        end
    end

    cost = sum((data_struct.lhs - data_struct.rhs).^2)
    return cost
end

"""
    obj_const_fft(params, equation, time_series, dt, strain, strain_deriv, stress, stress_deriv, data_struct; _insight::Bool = false)

Compute the cost function to minimize during parameter fitting, using the GL method with FFT.
"""
function obj_const_fft(params, equation, time_series, dt, strain, strain_deriv, stress, stress_deriv, data_struct; _insight::Bool = false)

    numerical_coeffs = get_coeffs(equation, params)
    @. data_struct.rhs = 0
    @. data_struct.lhs = 0
    
    L=length(strain)
    # Compute the rhs of the equation
    for term in numerical_coeffs.strain
        order = term[1]
        coeff = term[2]
        if order == 0.0
            data_struct.rhs .+= coeff * strain
        elseif order == 1.0
            data_struct.rhs .+= coeff * strain_deriv
        else
            generate_GL_weights(order, L, data_struct.weights)
            data_struct.rhs .+= coeff * compute_GL_frac_deriv_fftfilt(strain, data_struct.weights, order, dt)
        end
    end
    
    for term in numerical_coeffs.stress
        order = term[1]
        coeff = term[2]
        if order == 0.0
            data_struct.lhs .+= coeff * stress
        elseif order == 1.0
            data_struct.lhs .+= coeff * stress_deriv
        else
            generate_GL_weights(order, L, data_struct.weights)
            data_struct.lhs .+= coeff * compute_GL_frac_deriv_fftfilt(stress, data_struct.weights, order, dt)
        end
    end

    cost = sum((data_struct.lhs - data_struct.rhs).^2)
    return cost
end

function nlopt_constraint_wrapper(fw::FWConstraint)
    return (x::Vector, grad::Vector) -> begin
        fx = fw(x)
        if length(grad) > 0
            eps = 1e-8
            for i in eachindex(x)
                xh = copy(x)
                xh[i] += eps
                grad[i] = (fw(xh) - fx) / eps
            end
        end
        return fx
    end
end

struct SupportVectors
    rhs::Vector{Float64}
    lhs::Vector{Float64}
end


function get_coeffs(equation::DiffEqu{DiffScaFree}, params::Vector{Float64})

    strain_coeffs = Dict{Float64, Float64}()
    for term in equation.leftde
        order_val = term.order(params)
        coef_val = term.coef(params)
        strain_coeffs[order_val] = get(strain_coeffs, order_val, 0.0) + coef_val
    end

    stress_coeffs = Dict{Float64, Float64}()
    for term in equation.rightde
        order_val = term.order(params)
        coef_val = term.coef(params)
        stress_coeffs[order_val] = get(stress_coeffs, order_val, 0.0) + coef_val
    end

    return (strain = strain_coeffs, stress = stress_coeffs)
end


#=
--------------------------------
Predicting functions
--------------------------------
=#

function _modelpredictGL(data::RheoTimeData, equation ,diff_method)

    # Create the data structure to contain the right and left hand side of the equation, as well as the weights for the GL derivative
    data_struct = SupportVectors(
        zeros(length(data.t)),
        zeros(length(data.t))
    )

    # Define the derivative function and the input data based on the type of data provided
    if diff_method=="BD"
        deriv = derivBD
    elseif diff_method=="CD"
        deriv = derivCD
    end

    check = rheotimedatatype(data)
    if (check == strain_only)
        unknown = equation.rightde
        dependency = equation.leftde
        input = data.ϵ
        t = "strain"
    elseif (check == stress_only)
        unknown = equation.leftde
        dependency = equation.rightde
        input = data.σ
        t="stress"
    end

    n = length(data.t)
    dt = data.t[2] - data.t[1]
    deriv_data = deriv(input, data.t)   # First derivative of the input data

    prob = NumDiffProblem(dt=dt,order=0.5,n=length(data.t),method=GL())
    ws = init_workspace(prob)

    computed = zeros(n)
    denominator = 0.0

    # Compute the right-hand side of the equation based on the known terms
    for c in dependency
        if c.order == 0.0
            data_struct.rhs .+= c.coef * input
        elseif c.order == 1.0
            data_struct.rhs .+= c.coef * deriv_data
        else
            # generate_GL_weights(c.order, n, data_struct.weights)
            # data_struct.rhs .+= c.coef * compute_GL_frac_deriv(input, data_struct.deriv, data_struct.weights, c.order, dt)
            update_order!(prob,ws,c.order)
            compute!(prob.method,ws,input,prob)
            @. data_struct.rhs += c.coef * ws.deriv
        end
    end

    # Compute the denominator and the binomial coefficient weights for non-integer orders
    # Store the binomial coefficients in a dictionary to avoid redundant calculations for repeated orders
    bin_coeffs = Dict{Float64, Vector{Float64}}()
    for c in unknown
        if c.order == 0.0
            denominator += c.coef
        elseif c.order == 1.0
            denominator += c.coef / dt
        elseif round(c.order) != c.order && !haskey(bin_coeffs, c.order)
            denominator += c.coef / (dt^c.order)
            bin_coeffs[c.order] = zeros(n)
            generate_GL_weights(c.order, n, bin_coeffs[c.order])
        else
            denominator += c.coef / (dt^c.order)
        end
    end

    # Compute the first value of the computed array
    inv_denominator = 1.0 / denominator
    computed[1] = data_struct.rhs[1] * inv_denominator

    # Compute the rest of the values by updating the rhs with previous computed values for the unknown terms and dividing by the denominator
    for i in 2:n
        val = data_struct.rhs[i]

        for c in unknown
            if c.order == 1.0
                val += computed[i-1] * c.coef / dt

            elseif round(c.order) != c.order
                weights = bin_coeffs[c.order]
                coef_over_dt = c.coef / (dt^c.order) 

                conv_sum = 0.0
                @inbounds @simd for k in 1:(i-1)
                    conv_sum += weights[k+1] * computed[i-k]
                end

                val -= coef_over_dt * conv_sum

            elseif c.order != 0.0
                println("Order not implemented yet: $c")
            end
        end

        computed[i] = val * inv_denominator
    end

    return(computed, input, t)
end

# function _modelpredictFFT_stress(data::RheoTimeData, equation)

#     # Define the input data based on the type of data provided
#     unknown = equation.rightde
#     dependency = equation.leftde
#     input = data.ϵ

#     n  = length(data.t)
#     dt = data.t[2] - data.t[1]
#     L  = nextpow(2, 2n - 1)

#     # Pre-allocate each buffer
#     input_padded = zeros(Float64, L)
#     weights_padded = zeros(Float64, L)
#     input_padded[1:n] .= input

#     fft_size = L ÷ 2 + 1
#     input_fft = zeros(Complex{Float64}, fft_size)
#     weights_fft = zeros(Complex{Float64}, fft_size)
#     rhs_fft = zeros(Complex{Float64}, fft_size)
#     lhs_fft = zeros(Complex{Float64}, fft_size)
#     ifft_buf = zeros(Float64, L)

#     # Prepare FFTW plans using MEASURE as flag. Since the input is real, we can use rfft and irfft.
#     forward_plan = plan_rfft(input_padded; flags=FFTW.MEASURE)
#     inverse_plan = plan_irfft(rhs_fft, L; flags=FFTW.MEASURE)

#     # Transform input to frequency domain
#     mul!(input_fft, forward_plan, input_padded)

#     # Compute RHS in frequency domain (dependency terms)
#     fill!(rhs_fft, 0.0)
#     for c in dependency
#         inv_dt_pow = 1.0 / (dt^c.order)

#         # Prepare the weights buffer, which will vary for each term, based on the order of the derivative.
#         fill!(weights_padded, 0.0)
#         if c.order == 0.0
#             @inbounds @simd for i in 1:fft_size
#                 rhs_fft[i] += c.coef * input_fft[i]
#             end
#             continue
#         elseif c.order == 1.0
#             weights_padded[1] =  1.0
#             weights_padded[2] = -1.0
#         else
#             generate_GL_weights(c.order, n, weights_padded)
#         end

#         # Trasform weights to frequency domain
#         mul!(weights_fft, forward_plan, weights_padded)

#         # Update RHS in frequency domain using input and weights
#         @inbounds @simd for i in 1:fft_size
#             rhs_fft[i] += (c.coef * inv_dt_pow) * weights_fft[i] * input_fft[i]
#         end
#     end

#     # Compute LHS in frequency domain (unknown terms)
#     fill!(lhs_fft, 0.0)
#     for c in unknown
#         inv_dt_pow = 1.0 / (dt^c.order)

#         # Prepare the weights buffer, which will vary for each term, based on the order of the derivative.
#         fill!(weights_padded, 0.0)
#         if c.order == 0.0
#             @inbounds @simd for i in 1:fft_size
#                 lhs_fft[i] += c.coef
#             end
#             continue
#         elseif c.order == 1.0
#             weights_padded[1] =  1.0
#             weights_padded[2] = -1.0
#         else
#             generate_GL_weights(c.order, n, weights_padded)
#         end

#         # Trasform weights to frequency domain
#         mul!(weights_fft, forward_plan, weights_padded)

#         @inbounds @simd for i in 1:fft_size
#             lhs_fft[i] += (c.coef * inv_dt_pow) * weights_fft[i]
#         end
#     end

#     # Solve in frequency domain: output_fft = rhs_fft / lhs_fft
#     output_fft = zeros(Complex{Float64}, fft_size)
#     @inbounds @simd for i in 1:fft_size
#         output_fft[i] = rhs_fft[i] / lhs_fft[i]
#     end

#     # Return to time domain
#     mul!(ifft_buf, inverse_plan, output_fft)

#     return (ifft_buf[1:n], input, "strain")

# end


# function _modelpredictFFT_strain(data::RheoTimeData, equation)

#     # Define the input data based on the type of data provided
#     unknown = equation.leftde
#     dependency = equation.rightde
#     input = data.σ

#     n  = length(data.t)
#     dt = data.t[2] - data.t[1]
#     L  = nextpow(2, 2n - 1)

#     # Pre-allocate each buffer
#     input_padded = zeros(Float64, L)
#     weights_padded = zeros(Float64, L)
#     input_padded[1:n] .= input

#     input_padded[1] = 0.0
#     input_padded[2:n+1] .= input[1:n]   

#     fft_size = L ÷ 2 + 1
#     input_fft = zeros(Complex{Float64}, fft_size)
#     weights_fft = zeros(Complex{Float64}, fft_size)
#     rhs_fft = zeros(Complex{Float64}, fft_size)
#     lhs_fft = zeros(Complex{Float64}, fft_size)
#     ifft_buf = zeros(Float64, L)

#     min_order = minimum([c.order for c in unknown])

#     # Prepare FFTW plans using MEASURE as flag. Since the input is real, we can use rfft and irfft.
#     forward_plan = plan_rfft(input_padded; flags=FFTW.ESTIMATE)
#     inverse_plan = plan_irfft(rhs_fft, L; flags=FFTW.ESTIMATE)

#     # Transform input to frequency domain
#     mul!(input_fft, forward_plan, input_padded)

#     # Compute RHS in frequency domain (dependency terms)
#     fill!(rhs_fft, 0.0)
#     for c in dependency
#         new_order = c.order - min_order

#         inv_dt_pow = 1.0 / (dt^new_order)

#         # Prepare the weights buffer, which will vary for each term, based on the order of the derivative.
#         fill!(weights_padded, 0.0)
#         if new_order == 0.0
#             @inbounds @simd for i in 1:fft_size
#                 rhs_fft[i] += c.coef * input_fft[i]
#             end
#             continue
#         else
#             generate_GL_weights(new_order, n, weights_padded)
#         end

#         # Trasform weights to frequency domain
#         mul!(weights_fft, forward_plan, weights_padded)

#         # Update RHS in frequency domain using input and weights
#         @inbounds @simd for i in 1:fft_size
#             rhs_fft[i] += (c.coef * inv_dt_pow) * weights_fft[i] * input_fft[i]
#         end
#     end

#     # Compute LHS in frequency domain (unknown terms)
#     fill!(lhs_fft, 0.0)
#     for c in unknown
#         new_order = c.order - min_order

#         inv_dt_pow = 1.0 / (dt^new_order)

#         # Prepare the weights buffer, which will vary for each term, based on the order of the derivative.
#         fill!(weights_padded, 0.0)
#         if new_order == 0.0
#             @inbounds @simd for i in 1:fft_size
#                 lhs_fft[i] += c.coef
#             end
#             continue
#         else
#             generate_GL_weights(new_order, n, weights_padded)
#         end

#         # Trasform weights to frequency domain
#         mul!(weights_fft, forward_plan, weights_padded)

#         @inbounds @simd for i in 1:fft_size
#             lhs_fft[i] += (c.coef * inv_dt_pow) * weights_fft[i]
#         end
#     end

#     # Solve in frequency domain: output_fft = rhs_fft / lhs_fft
#     output_fft = zeros(Complex{Float64}, fft_size)
#     @inbounds for i in 1:fft_size
#         output_fft[i] = rhs_fft[i] / lhs_fft[i]
#     end

#     # Return to time domain
#     mul!(ifft_buf, inverse_plan, output_fft)

#     return (ifft_buf[1:n], input, "stress")
# end


"""
    modelpredictGL(data::RheoTimeData, model::RheoModel)

Given an incomplete data set (only either stress or strain missing) and model with values substituted into
parameters (`RheoModel`), return a new dataset based on the model using the Grunwald-Letnikov algorithm for the fractional derivatives.
A complete `RheoTimeData` of type `strain_and_stress` is returned.
"""
function modelpredict(data::RheoTimeData, model::RheoModel,predtype::Differential = Differential())

    check = rheotimedatatype(data)
    @assert (check == strain_only)||(check == stress_only) "Need either strain only or stress only data. Data provided: " * string(check)
    if check == strain_only
        sigma, epsilon, pred_mod = _modelpredictGL(data, model.C, "BD")
    else check == stress_only
        epsilon, sigma, pred_mod = _modelpredictGL(data, model.C, "BD")
    end
    log = logadd_process(data, :modelpredict, params=(model,), 
                         comment="Predicted data - modulus: $pred_mod, parameters:$(model.fixedparams)" ) 

    return RheoTimeData(sigma, epsilon, data.t, log)

end

"""
    modelpredictFFT(data::RheoTimeData, model::RheoModelDiff)

Given an incomplete data set (only either stress or strain missing) and model with values substituted into
parameters (`RheoModel`), return a new dataset based on the model using the Fast Fourier Transform.
A complete `RheoTimeData` of type `strain_and_stress` is returned.
"""
# function modelpredictFFT(data::RheoTimeData, model::RheoModel)

#     check = rheotimedatatype(data)
#     @assert (check == strain_only)||(check == stress_only) "Need either strain only or stress only data. Data provided: " * string(check)
#     if check == strain_only
#         sigma, epsilon, pred_mod = _modelpredictFFT_stress(data, model.C)
#     else check == stress_only
#         epsilon, sigma, pred_mod = _modelpredictFFT_strain(data, model.C)
#     end
#     log = logadd_process(data, :modelpredict, params=(model,), 
#                          comment="Predicted data - modulus: $pred_mod, parameters:$(model.fixedparams)" ) 

#     return RheoTimeData(sigma, epsilon, data.t, log)

# end