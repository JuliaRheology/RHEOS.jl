using NumFracDiff


function modelfitDiff(data::RheoTimeData, 
    model::RheoModelClass,
    modloading::LoadingType;
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
    optmaxeval::Union{Integer,Nothing} = nothing) where T <: Integer

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
@timed leastsquares_initLHSRHS_parallel(   p0a,
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
                            model._constraint;
                            method= method,
                            insight = verbose,
                            constant_sampling = is_constant,
                            singularity = false,
                            rel_tol_x = rel_tol_x,
                            rel_tol_f = rel_tol_f,
                            indweights = weights,
                            optmethod = Symbol(optmethod),
                            opttimeout = opttimeout,
                            optmaxeval = optmaxeval)

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

return (RheoModel(model, nt), minf);
end 



function leastsquares_initLHSRHS_parallel(params_init::Vector{RheoFloat},
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
    constraint::Union{Vector{FWConstraint},Nothing};
    method=RL(),
    insight::Bool = false,
    constant_sampling::Bool=true,
    singularity::Bool = false,
    rel_tol_x::Union{Real,Nothing} = nothing,
    rel_tol_f::Union{Real,Nothing} = nothing,
    indweights = nothing,
    optmethod::Symbol = :LN_SBPLX,
    opttimeout::Union{Real,Nothing} = nothing,
    optmaxeval::Union{Integer,Nothing} = nothing)


# initialise NLOpt.Opt object with :LN_SBPLX Subplex algorithm
if !(optmethod in [:LN_AUGLAG, :LN_COBYLA]) && constraint ≠ nothing
optmethod = :LN_COBYLA
end
println(optmethod)
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
# ws = NumDiffWorkspace(zeros(RheoFloat,length(time_series)),zeros(RheoFloat,length(time_series)))
ws = init_workspace(prob)
min_objective!(opt, (params, grad) -> obj_const(params, equation,
                            time_series, dt, 
                            strain, strain_deriv, 
                            stress, stress_deriv,data_struct, prob, ws;
                            _insight = insight))


if constraint ≠ nothing
for c in eachindex(constraint)
inequality_constraint!(opt,nlopt_constraint_wrapper(constraint[c]), 1e-8)
end
end


# minimise objective func, minx are the parameters resulting in minimum
(minf, minx, ret) = NLopt.optimize(opt, params_init)

# return all
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
            # if SM
            #     L=calculate_optimal_L(strain,order)
            # end
            # generate_GL_weights(order,L,data_struct.weights)
            # compute_GL_frac_deriv(strain, data_struct.deriv, data_struct.weights, order, dt; L=L)
            # @. data_struct.rhs += coeff * data_struct.deriv
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
            # if SM
            #     L = calculate_optimal_L(stress,order)
            # end
            # generate_GL_weights(order,L,data_struct.weights)
            # compute_GL_frac_deriv(stress, data_struct.deriv, data_struct.weights, order, dt; L=L)
            # @. data_struct.lhs += coeff * data_struct.deriv

            update_order!(prob,ws,order)
            compute!(prob.method,ws,stress,prob)
            @. data_struct.rhs += coeff * ws.deriv
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