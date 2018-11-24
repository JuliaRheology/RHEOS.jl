#!/usr/bin/env julia

"""
    RheologyData(σ::Vector{T}, ϵ::Vector{T}, t::Vector{T}, sampling::String, log::Vector{String}) where T<:Real

RheologyData struct contains stress, strain and time data.

If preferred, an instance can be generated manually by just providing the three data
vectors in the right order, sampling type will be checked automatically. If loading
partial data (either stress or strain), fill the other vector as a vector of zeros
of the same length as the others.

# Fields

- σ: stress
- ϵ: strain
- t: time
- sampling: sampling type, either "constant" or "variable"
- log: a log of struct's events, e.g. preprocessing
"""
struct RheologyData

    σ::Vector{T} where T<:Real
    ϵ::Vector{T} where T<:Real
    t::Vector{T} where T<:Real

    sampling::String

    log::Vector{String}

end

function RheologyData(colnames::Vector{String}, data1::Vector{T}, data2::Vector{T}, data3::Vector{T}=zeros(length(data2)); filedir::String="none", log::Vector{String}=Array{String}(undef, 0)) where T<:Real

    # checks
    @assert length(data1)==length(data2) "Data arrays must be same length"
    @assert length(data1)==length(data3) "Data arrays must be same length"

    # get data in correct variables
    data = [data1, data2, data3]
    σ = Vector{T}(undef, 0)
    ϵ = Vector{T}(undef, 0)
    t = Vector{T}(undef, 0)

    for (i, v) in enumerate(colnames)
        if v == "stress"
            σ = data[i]
        elseif v == "strain"
            ϵ = data[i]
        elseif v == "time"
            t = data[i]
        else
            @assert false "Incorrect Column Names"
        end
    end

    # test for NaNs
    newstartingval = 1
    for i in 1:length(σ)
        if !isnan(σ[i]) && !isnan(ϵ[i])
            newstartingval = i
            break
        end
    end

    # adjust starting point accordingly to remove NaNs in σ, ϵ
    σ = σ[newstartingval:end]
    ϵ = ϵ[newstartingval:end]
    t = t[newstartingval:end]

    # readjust time to account for NaN movement and/or negative time values
    # t = t - minimum(t)

    # set up log for three cases, file dir given, derived from other data so filedir
    # in log already, no file name or log given.
    if filedir != "none" || length(log)==0
        if length(colnames)<3
            push!(log, "partial data loaded from:")
        elseif length(colnames)==3
            push!(log, "complete data loaded from:")
        end
        push!(log, filedir)
    end

    sampling = constantcheck(t) ? "constant" : "variable"

    # return class with all fields initialised
    RheologyData(σ, ϵ, t, sampling, log)

end

function RheologyData(σ::Vector{T}, ϵ::Vector{T}, t::Vector{T}; log::Vector{String}=["Manually Created."]) where T<:Real

    sampling = constantcheck(t) ? "constant" : "variable"

    RheologyData(σ, ϵ, t, sampling, log)

end

function RheologyData(data::Vector{T}, t::Vector{T}, log::Vector{String}) where T<:Real
    # used when generating data so always constant
    RheologyData(data, data, t, "constant", log)

end

function +(self1::RheologyData, self2::RheologyData)

    @assert sampleratecompare(self1.t, self2.t) "Step size must be same for both datasets"

    # get time array
    if length(self1.t) >= length(self2.t)
        t = self1.t
    else
        t = self2.t
    end

    # init data array and fill by summing over each argument's indices
    σ = zeros(length(t))
    ϵ = zeros(length(t))
    # sum with last value propagating (hanging)
    for i in 1:length(t)
        if i<=length(self1.t)
            σ[i] += self1.σ[i]
            ϵ[i] += self1.ϵ[i]
        else
            σ[i] += self1.σ[end]
            ϵ[i] += self1.ϵ[end]
        end

        if i<=length(self2.t)
            σ[i] += self2.σ[i]
            ϵ[i] += self2.ϵ[i]
        else
            σ[i] += self2.σ[end]
            ϵ[i] += self2.ϵ[end]
        end
    end

    # log
    log = vcat(self1.log, self2.log, ["previous two logs added"])

    RheologyData(σ, ϵ, t, "constant", log)

end

function -(self1::RheologyData, self2::RheologyData)

    @assert sampleratecompare(self1.t, self2.t) "Step size must be same for both datasets"

    # get time array
    if length(self1.t) >= length(self2.t)
        t = self1.t
    else
        t = self2.t
    end

    # init data array and fill by summing over each argument's indices
    σ = zeros(length(t))
    ϵ = zeros(length(t))
    # sum with last value propagating (hanging)
    for i in 1:length(t)
        if i<=length(self1.t)
            σ[i] += self1.σ[i]
            ϵ[i] += self1.ϵ[i]
        else
            σ[i] += self1.σ[end]
            ϵ[i] += self1.ϵ[end]
        end

        if i<=length(self2.t)
            σ[i] -= self2.σ[i]
            ϵ[i] -= self2.ϵ[i]
        else
            σ[i] -= self2.σ[end]
            ϵ[i] -= self2.ϵ[end]
        end
    end

    # log
    log = vcat(self1.log, self2.log, ["previous two logs added"])

    RheologyData(σ, ϵ, t, "constant", log)

end

function *(self1::RheologyData, self2::RheologyData)

    @assert sampleratecompare(self1.t, self2.t) "Step size must be same for both datasets"

    # get time array
    if length(self1.t) >= length(self2.t)
        t  = self1.t
    else
        t = self2.t
    end

    # init data array and fill by summing over each argument's indices
    σ = zeros(length(t))
    ϵ = zeros(length(t))
    # sum with last value propagating (hanging)
    for i in 1:length(t)
        if i<=length(self1.t) && i<=length(self2.t)
            σ[i] = self1.σ[i]*self2.σ[i]
            ϵ[i] = self1.ϵ[i]*self2.ϵ[i]

        elseif i<=length(self1.t) && i>length(self2.t)
            σ[i] = self1.σ[i]
            ϵ[i] = self1.ϵ[i]

        elseif i>length(self1.t) && i<=length(self2.t)
            σ[i] = self2.σ[i]
            ϵ[i] = self2.ϵ[i]

        end
    end

    # log
    log = vcat(self1.log, self2.log, ["previous two logs added"])

    RheologyData(σ, ϵ, t, "constant", log)

end

function *(self::RheologyData, operand::Real)

    ϵ = self.ϵ*operand
    σ = self.σ*operand

    # log
    log = vcat(self.log, ["multiplied data by $operand"])

    RheologyData(σ, ϵ, self.t, "constant", log)

end

function *(operand::Real, self::RheologyData)

    # multiplication commutes so call function as defined for opposite operand order
    return self*operand

end

"""
    RheologyModel(G::T, J::T, Gp::T, Gpp::T, parameters::Vector{S<:Real} log::Vector{String}) where T<:Function

RheologyModel contains a model's various moduli, parameters, and log of activity.

For incomplete models, an alternative constructor is available where all arguments
are keyword arguments and moduli not provided default to a null modulus which
always returns [-1.0].

# Fields

- G: Relaxation modulus
- J: Creep modulus
- Gp: Storage modulus
- Gpp: Loss modulus
- parameters: Used for predicting and as default starting parameters in fitting
- log: a log of struct's events, e.g. what file it was fitted to
"""
struct RheologyModel

    G::Function

    J::Function

    Gp::Function

    Gpp::Function

    parameters::Vector{T} where T<:Real

    log::Vector{String}

end

function null_modulus(t::Vector{T}, params::Vector{T}) where T<:Real
    return [-1.0]
end

RheologyModel(;G::Function = null_modulus,
               J::Function = null_modulus,
               Gp::Function = null_modulus,
               Gpp::Function = null_modulus,
               params::Vector{T} where T<:Real = [-1.0],
               log::Vector{String} = ["model created by user with parameters $params"]) = RheologyModel(G, J, Gp, Gpp, params, log)

"""
    RheologyDynamic(Gp::Vector{T}, Gpp::Vector{T}, ω::Vector{T}, log::Vector{String}) where T<:Real

RheologyDynamic contains storage modulus, loss modulus and frequency data.

If preferred, an instance can be generated manually by just providing the three data
vectors in the right order.

# Fields

- Gp: storage modulus
- Gpp: loss modulus
- ω: frequency
- log: a log of struct's events, e.g. preprocessing
"""
struct RheologyDynamic

    # original data
    Gp::Vector{T} where T<:Real
    Gpp::Vector{T} where T<:Real
    ω::Vector{T} where T<:Real

    # operations applied, stores history of which functions (including arguments)
    log::Vector{String}

end

function RheologyDynamic(colnames::Vector{String}, data1::Vector{T}, data2::Vector{T}, data3::Vector{T}; filedir::String="none", log::Vector{String}=Array{String}(undef, 0)) where T<:Real

    # checks
    @assert length(data1)==length(data2) "Data arrays must be same length"
    @assert length(data1)==length(data3) "Data arrays must be same length"

    # get data in correct variables
    data = [data1, data2, data3]
    Gp = Array{T}(undef, 0)
    Gpp = Array{T}(undef, 0)
    ω = Array{T}(undef, 0)

    for (i, v) in enumerate(colnames)
        if v == "Gp"
            Gp = data[i]
        elseif v == "Gpp"
            Gpp = data[i]
        elseif v == "frequency"
            ω = data[i]
        else
            @assert false "Incorrect Column Names"
        end
    end

    # set up log for three cases, file dir given, derived from other data so filedir
    # in log already, no file name or log given.
    if filedir != "none" || length(log)==0
        push!(log, string("complete data loaded from:", filedir))
    end



    # return class with all fields initialised
    RheologyDynamic(Gp, Gpp, ω, log)

end

function RheologyDynamic(Gp::Vector{T}, Gpp::Vector{T}, ω::Vector{T}; log::Vector{String}=["Manually Created."]) where T<:Real

    RheologyDynamic(Gp, Gpp, ω, log)

end