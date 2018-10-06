paramsDefault = [0.5, 0.5, 1.0, 0.2]

function vectorizedSLS(t::Array{T,1}, params::Array{T,1}) where T<:Real
    cₐ, a, kᵦ, kᵧ = params

    G = kᵦ*mittleff.(a, -(kᵦ/cₐ)*t.^a) .+ kᵧ
end;

function numberSLS(t::T, params::Vector{T}) where T<:Real
    cₐ, a, kᵦ, kᵧ = params

    G = kᵦ*mittleff(a, -(kᵦ/cₐ)*t^a) + kᵧ
end

broadcastSLS(t::Vector{T}, params::Vector{T}) where T<:Real = numberSLS.(t, (params,));
