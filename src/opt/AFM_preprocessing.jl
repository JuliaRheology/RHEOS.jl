function contact_point(self::RheoTimeData, interface::Interface, cp::String, param::NamedTuple; savelog = true,  sec = :)

    contactmodels = Dict("threshold" => contact_threshold,
                         "hertz" => contact_hertz);
    
    data = self[interface];
    
    cp_index = contactmodels[cp](data[:f], data[:d], param);

    t = self.t[cp_index:end] .- self.t[cp_index];
    d = data[:d][cp_index:end] .- data[:d][cp_index];
    f = data[:f][cp_index:end] .- data[:f][cp_index];

    ϵ,σ = interface.to_ϵσ(d, f);
    
    log = self.log == nothing ? nothing : [self.log; RheoLogItem( (type=:process, funct=:contact_point, params=(cp = cp, param = param ), keywords=(sec = sec, ) ),
    (comment="Contact point detection and data shift",) ) ]

 return RheoTimeData(σ, ϵ, t, log)

end


"""
    contact_threshold(f::Array{Float64,1}, d::Array{Float64,1}; _param::NamedTuple)

Applies force threshold to find approximate contact point.
"""
function contact_threshold(force::Array{RheoFloat}, disp::Array{RheoFloat}, param::NamedTuple)

    # get index of contact
    cp_index = findfirst(force .>= param.threshold);

    # return 
    return cp_index

end


"""
    contact_hertz(f::Array{Float64,1}, δ::Array{Float64,1}; _param::Float64 = NaN64)

    Apply spherical Hertz contact model to find approximate contact point.
"""
function contact_hertz(f::Array{RheoFloat}, δ::Array{RheoFloat}, param::NamedTuple)

    # regularise displacement
    δ = δ .- minimum(δ);

    # radius of sphere
    R = param.R;
    
    # poisson ratio - assume incompressible
    ν = 0.5;

    # get approximate δ0 and YM values
    (init_δ₀, init_YM, tilt_approx, offset_approx) =  hertz_approx(f, δ, R, ν);

    # least squares fit
    opt = Opt(:LN_SBPLX, 4); # 4 parameters to fit
    
    bounds_param = 0.25 # +/- limits allowance from initial contact point guess
    lower_bounds!(opt, [init_δ₀*(1 - bounds_param), 0.0, -5.0, -5.0]);
    upper_bounds!(opt, [init_δ₀*(1 + bounds_param), 1e10, 5.0, 5.0]);

    xtol_rel!(opt,1e-5);

    min_objective!(opt, (params, grad) -> hertz_cost(f, δ, R, ν, params[1], params[2], params[3], params[4]));

    params_init = [init_δ₀, init_YM, tilt_approx, offset_approx];


    (minf, minx, ret) = optimize(opt, params_init);

    # return index of contact point
    δ₀_index = argmin(abs.(δ .- minx[1]));
    
end


###########################################
# Convenience functions for contact_hertz #
###########################################
"""
Generate approximate Hertz fit a, b, E, cp to use as initial conditions in proper fit.
"""

function hertz_approx(f::Array{RheoFloat,1}, δ::Array{RheoFloat,1}, R::RheoFloat, ν::RheoFloat)
    
    deriv = derivCD(derivCD(f, δ),δ);
    value_max, index_max = findmax(deriv);

    plot(δ, deriv)

    #coeff = (1-index_max/length(f))

    # approximate contact point, assume linear near end of approach section
    # and trace line down to axis
    #last_section = round(Int64, coeff*length(f)):length(f)
    if index_max > 5    
        last_section = (index_max-5):length(f);
    else
        last_section = (index_max):length(f);
    end

    avg_gradient = mean(derivBD(f[last_section], δ[last_section]));

    c = f[end] - avg_gradient*δ[end];

    approx_δ₀ = -c/avg_gradient;

    # transform to find approx YM
    δᵦ = (heaviside(δ .- approx_δ₀).*δ).^(3/2);

    prefactor = (4/3)*sqrt(R)/(1 - ν^2);

    avg_gradient_transformed = mean(derivBD(f[last_section], δᵦ[last_section]));

    approx_YM = avg_gradient_transformed/prefactor;

    # estimate offset and tilt
    first_section = 1:round(Int64, 0.5*length(f));

    df_dh = derivBD(f, δ);

    tilt_approx = mean(df_dh[first_section]);

    offset_approx = mean(f[first_section]);

    if approx_δ₀<0
        approx_δ₀ = 0
    end

    if approx_YM<0
        approx_YM = 0
    end

    # return
    approx_δ₀, approx_YM, tilt_approx, offset_approx;
end

"""
Heaviside function
"""
function heaviside(x)
   0.5*(sign.(x) .+ 1);
end

"""
Hertz spherical model 
"""
function hertz_model(δ::Array{RheoFloat,1}, R::RheoFloat, ν::RheoFloat, δ₀::RheoFloat, E::RheoFloat, tilt::RheoFloat, offset::RheoFloat)
    # define prefactor for Hertz sphere model
    prefactor = (4/3)*sqrt(R)/(1 - ν^2);

    δ_shift = δ - δ₀;

    # force
    F = tilt*δ + offset + prefactor*E*((δ_shift.*heaviside(δ_shift)).^(3/2));
end

"""
Hertz cost function
"""
function hertz_cost(f::Array{RheoFloat,1}, δ::Array{RheoFloat,1}, R::RheoFloat, ν::RheoFloat, δ₀::RheoFloat, E::RheoFloat, tilt::RheoFloat, offset::RheoFloat)
    # difference
    Δ = f - hertz_model(δ, R, ν, δ₀, E, tilt, offset);

    # regularisation parameter for offset and tilt
    λ = 0.1;

    # sum of squares
    cost = sum(Δ.^2); # Maybe not needed + λ*(tilt + offset)^2
end

