#!/usr/bin/env julia

export RheoTimeData, RheoTimeDataType, RheoFreqData, RheoFreqDataType, check_time_data_consistency
export LoadingType, strain_imposed, stress_imposed
export TimeDataType, time_only, strain_only, stress_only, strain_and_stress
export FreqDataType, invalid_freq_data, frec_only, with_modulus









# Empty vector mostly used as default parameters to indicate missing/unspecified data.
empty_rheodata_vector=RheoFloat[]


"""
    RheoTimeData(;σ::Vector{T1}, ϵ::Vector{T2}, t::Vector{T3}, log::Vector{String}) where {T1<:Real, T2<:Real, T3<:Real}

RheoTimeData struct contains stress, strain and time data.

If preferred, an instance can be generated manually by just providing the three data
vectors in the right order, sampling type will be checked automatically. If loading
partial data (either stress or strain), fill the other vector as a vector of zeros
of the same length as the others.

# Fields

- σ: stress
- ϵ: strain
- t: time
- log: a log of struct's events, e.g. preprocessing
"""
struct RheoTimeData

    σ::Vector{RheoFloat}
    ϵ::Vector{RheoFloat}
    t::Vector{RheoFloat}

    log::Vector{String}
    #log::Dict{Any,Any}

end

function RheoTimeData(;ϵ::Vector{T1} = empty_rheodata_vector, σ::Vector{T2} = empty_rheodata_vector, t::Vector{T3} = empty_rheodata_vector, source="User provided data.")  where {T1<:Real, T2<:Real, T3<:Real}
    s_datatype = string("Data type: ",check_time_data_consistency(t,ϵ,σ))
    RheoTimeData(convert(Vector{RheoFloat},σ), convert(Vector{RheoFloat},ϵ), convert(Vector{RheoFloat},t), [source,s_datatype])   #Dict{Any,Any}("source"=> source, "datatype"=>check_time_data_consistency(t,ϵ,σ))
end



@enum TimeDataType invalid_time_data=-1 time_only=0 strain_only=1 stress_only=2 strain_and_stress=3

function check_time_data_consistency(t,e,s)
    @assert (length(t)>0)  "Time data empty"

    sdef=(s != empty_rheodata_vector)
    edef=(e != empty_rheodata_vector)
    if (sdef && edef)
        @assert (length(s)==length(t)) && (length(e)==length(t)) "Time data length inconsistent"
        return strain_and_stress
    end

    if (sdef && (!edef))
        @assert (length(s)==length(t)) "Time data length inconsistent"
        return stress_only
    end

    if (edef && (!sdef))
        @assert (length(e)==length(t)) "Time data length inconsistent"
        return strain_only
    end

    if ((!edef) && (!sdef))
        return time_only
    end

    return invalid_time_data

end

function RheoTimeDataType(d::RheoTimeData)
    return check_time_data_consistency(d.t,d.ϵ,d.σ)
end

@enum LoadingType strain_imposed=1 stress_imposed=2




"""
    RheoFreqData(Gp::Vector{T1}, Gpp::Vector{T2}, ω::Vector{T3}, log::Vector{String}) where {T1<:Real, T2<:Real, T3<:Real}

RheologyDynamic contains storage modulus, loss modulus and frequency data.

If preferred, an instance can be generated manually by just providing the three data
vectors in the right order.

# Fields

- Gp: storage modulus
- Gpp: loss modulus
- ω: frequency
- log: a log of struct's events, e.g. preprocessing
"""
struct RheoFreqData

    # Complex modulus data
    Gp::Vector{RheoFloat}
    Gpp::Vector{RheoFloat}
    ω::Vector{RheoFloat}

    # operations applied, stores history of which functions (including arguments)
    log::Vector{String}

end


function RheoFreqData(;Gp::Vector{T1} = empty_rheodata_vector, Gpp::Vector{T2} = empty_rheodata_vector, ω::Vector{T3} = empty_rheodata_vector, info="User provided data.")  where {T1<:Real, T2<:Real, T3<:Real}
    s_datatype = string("Data type: ",check_freq_data_consistency(ω,Gp,Gpp))
    RheoFreqData(convert(Vector{RheoFloat},Gp), convert(Vector{RheoFloat},Gpp), convert(Vector{RheoFloat},ω), [info,s_datatype])
end


@enum FreqDataType invalid_freq_data=-1 frec_only=0 with_modulus=1

function check_freq_data_consistency(o,gp,gpp)
    @assert (length(o)>0)  "Freq data empty"

    gpdef=(gp != empty_rheodata_vector)
    gppdef=(gpp != empty_rheodata_vector)
    if (gpdef && gppdef)
        @assert (length(gp)==length(o)) && (length(gpp)==length(o)) "Data length inconsistent"
        return with_modulus
    end

    if ((!gpdef) && (!gppdef))
        return frec_only
    end

    return invalid_freq_data

end

function RheoFreqDataType(d::RheoFreqData)
    return check_freq_data_consistency(d.ω,d.Gp,d.Gpp)
end
















function +(self1::RheoTimeData, self2::RheoTimeData)

    type1 = RheoTimeDataType(self1)
    type2 = RheoTimeDataType(self2)
    @assert (type1!=invalid_time_data) "Addition error: first parameter invalid"
    @assert (type2!=invalid_time_data) "Addition error: second parameter invalid"
    @assert (type1==type2) "Addition error: parameters inconsistent"
    @assert (type1!=time_only) "Addition error: time only data cannot be added"
    @assert (self1.t == self2.t) "Addition error: timelines inconsistent"

    # Operation on the logs - not so clear what it should be
    log = vcat(self1.log, self2.log, ["previous two logs added"])

    if (type1==strain_only) && (type2==strain_only)
        return RheoTimeData(empty_rheodata_vector, self1.ϵ+self2.ϵ, self1.t, log)
    end

    if (type1==stress_only) && (type2==stress_only)
        return RheoTimeData(self1.σ+self2.σ, empty_rheodata_vector, self1.t, log)
    end

    if (type1==strain_and_stress) && (type2==strain_and_stress)
        return RheoTimeData(self1.σ+self2.σ, self1.ϵ+self2.ϵ, self1.t, log)
    end

end


function -(self1::RheoTimeData, self2::RheoTimeData)

    type1 = RheoTimeDataType(self1)
    type2 = RheoTimeDataType(self2)
    @assert (type1!=invalid_time_data) "Addition error: first parameter invalid"
    @assert (type2!=invalid_time_data) "Addition error: second parameter invalid"
    @assert (type1==type2) "Addition error: parameters inconsistent"
    @assert (type1!=time_only) "Addition error: time only data cannot be added"
    @assert (self1.t == self2.t) "Addition error: timelines inconsistent"

    # Operation on the logs - not so clear what it should be
    log = vcat(self1.log, self2.log, ["previous two logs added"])

    if (type1==strain_only) && (type2==strain_only)
        return RheoTimeData(empty_rheodata_vector, self1.ϵ-self2.ϵ, self1.t, log)
    end

    if (type1==stress_only) && (type2==stress_only)
        return RheoTimeData(self1.σ-self2.σ, empty_rheodata_vector, self1.t, log)
    end

    if (type1==strain_and_stress) && (type2==strain_and_stress)
        return RheoTimeData(self1.σ-self2.σ, self1.ϵ-self2.ϵ, self1.t, log)
    end

end


function -(self1::RheoTimeData)

    type1 = RheoTimeDataType(self1)
    @assert (type1!=invalid_time_data) "unary - error: parameter invalid"
    @assert (type1!=time_only) "unary - error: time only data cannot be manipulated this way"

    # log
    log = vcat(self1.log, ["multiplied by -1"])

    return RheoTimeData(-self1.σ, -self1.ϵ, self1.t, log)

end



function *(operand::Real, self1::RheoTimeData)

    type1 = RheoTimeDataType(self1)
    @assert (type1!=invalid_time_data) "* error: parameter invalid"
    @assert (type1!=time_only) "* error: time only data cannot be manipulated this way"

    # log
    log = vcat(self1.log, ["multiplied data by $operand"])

    return RheoTimeData(operand*self1.σ, operand*self1.ϵ, self1.t, log)
end

function *(self1::RheoTimeData, operand::Real)
    return operand*self1
end





#  Time shift operator >>
#  shift time by a certain amount, trash the end and pad at the start with 0

#  Union operator |
#  Combine stress_only data and strain_only data into stress_strain













struct RheoModelClass

    name::String
    params::Vector{Symbol}

    G::FunctionWrapper{RheoFloat,Tuple{RheoFloat,Array{RheoFloat,1}}}
    Ga::FunctionWrapper{Array{RheoFloat,1},Tuple{Array{RheoFloat,1},Array{RheoFloat,1}}}
    J::FunctionWrapper{RheoFloat,Tuple{RheoFloat,Array{RheoFloat,1}}}
    Ja::FunctionWrapper{Array{RheoFloat,1},Tuple{Array{RheoFloat,1},Array{RheoFloat,1}}}
    Gp::FunctionWrapper{RheoFloat,Tuple{RheoFloat,Array{RheoFloat,1}}}
    Gpa::FunctionWrapper{Array{RheoFloat,1},Tuple{Array{RheoFloat,1},Array{RheoFloat,1}}}
    Gpp::FunctionWrapper{RheoFloat,Tuple{RheoFloat,Array{RheoFloat,1}}}
    Gppa::FunctionWrapper{Array{RheoFloat,1},Tuple{Array{RheoFloat,1},Array{RheoFloat,1}}}

    constraint::FunctionWrapper{Bool,Tuple{Array{RheoFloat,1}}}

    info::String
    expressions::NamedTuple

end


function Base.show(io::IO, m::RheoModelClass)
    print(io, "\nModel name: $(m.name)")
    ps=join([string(s) for s in m.params], ", ", " and ")
    print(io, "\n\nFree parameters: $ps\n")
    print(io, m.info)
    return
end


export rheoconv,invLaplace

rheoconv(t::Real) = RheoFloat(t)
rheoconv(t::Array{T,1}) where T<:Real = convert(Vector{RheoFloat},t)

invLaplace(f::Function, t::Array{RheoFloat}) = InverseLaplace.talbotarr(f, t)
invLaplace(f::Function, t::RheoFloat) = InverseLaplace.talbot(f, t)

# Define Null Expression, and make it default parameter
# Model name and parameters should be required --> assert
# info should have generic default value, such as "Model with $n params named $s..."

function RheoModelClass(;name::String,
                         p::Array{Symbol}=[],
                         G::Expr = quote NaN end,
                         J::Expr = quote NaN end,
                         Gp::Expr = quote NaN end,
                         Gpp::Expr = quote NaN end,
                         constraint::Expr = quote true end,
                         info)

    # Expression to unpack parameter array into suitably names variables in the moduli expressions
    unpack_expr = Meta.parse(string(join(string.(p), ","), ",=params"))
    expressions = (G=G,J=J,Gp=Gp,Gpp=Gpp,constraint=constraint)

    @eval return(RheoModelClass($name, $p,
        ((t,params) -> begin $unpack_expr; $G; end)                 |> FunctionWrapper{RheoFloat,Tuple{RheoFloat,Array{RheoFloat,1}}},
        ((ta,params) -> begin $unpack_expr; [$G for t in ta]; end)  |> FunctionWrapper{Array{RheoFloat,1},Tuple{Array{RheoFloat,1},Array{RheoFloat,1}}},
        ((t,params) -> begin $unpack_expr; $J; end)                 |> FunctionWrapper{RheoFloat,Tuple{RheoFloat,Array{RheoFloat,1}}},
        ((ta,params) -> begin $unpack_expr; [$J for t in ta]; end)  |> FunctionWrapper{Array{RheoFloat,1},Tuple{Array{RheoFloat,1},Array{RheoFloat,1}}},
        ((t,params) -> begin $unpack_expr; $Gp; end)                |> FunctionWrapper{RheoFloat,Tuple{RheoFloat,Array{RheoFloat,1}}},
        ((ta,params) -> begin $unpack_expr; [$Gp for t in ta]; end) |> FunctionWrapper{Array{RheoFloat,1},Tuple{Array{RheoFloat,1},Array{RheoFloat,1}}},
        ((t,params) -> begin $unpack_expr; $Gpp; end)               |> FunctionWrapper{RheoFloat,Tuple{RheoFloat,Array{RheoFloat,1}}},
        ((ta,params) -> begin $unpack_expr; [$Gpp for t in ta]; end) |> FunctionWrapper{Array{RheoFloat,1},Tuple{Array{RheoFloat,1},Array{RheoFloat,1}}},
        (params -> begin $unpack_expr; $constraint; end)            |> FunctionWrapper{Bool,Tuple{Array{RheoFloat,1}}},
        $info, $expressions) )
end






struct RheoModel

    G::FunctionWrapper{RheoFloat,Tuple{RheoFloat}}
    Ga::FunctionWrapper{Array{RheoFloat,1},Tuple{Array{RheoFloat,1}}}
    J::FunctionWrapper{RheoFloat,Tuple{RheoFloat}}
    Ja::FunctionWrapper{Array{RheoFloat,1},Tuple{Array{RheoFloat,1}}}
    Gp::FunctionWrapper{RheoFloat,Tuple{RheoFloat}}
    Gpa::FunctionWrapper{Array{RheoFloat,1},Tuple{Array{RheoFloat,1}}}
    Gpp::FunctionWrapper{RheoFloat,Tuple{RheoFloat}}
    Gppa::FunctionWrapper{Array{RheoFloat,1},Tuple{Array{RheoFloat,1}}}
    expressions::NamedTuple


    params::NamedTuple
    info::String
    log::Vector{String}

end



# function Base.show(io::IO, m::RheoModelClass)
#     print(io, "\nModel name: $(m.name)")
#     ps=join([string(s) for s in m.params], ", ", " and ")
#     print(io, "\n\nFree parameters: $ps\n")
#     print(io, m.info)
#     return
# end


# Cool replacement function inspired from
# https://stackoverflow.com/questions/29778698/julia-lang-metaprogramming-turn-expression-into-function-with-expression-depend
# This probably could go in its own module as this is very useful.
#

expr_replace!(ex, s, v) = ex
expr_replace!(ex::Symbol, s, v) = s == ex ? v : ex

function expr_replace!(ex::Expr, s, v)
    for i=1:length(ex.args)
        ex.args[i] = expr_replace!(ex.args[i], s, v)
    end
    return ex
end

expr_replace(ex, s, v) = expr_replace!(copy(ex), s, v)

function expr_replace(ex, nt)
    e = copy(ex)
    k=keys(nt)
    v=values(nt)
    for i in 1:length(nt)
        expr_replace!(e, k[i], v[i])
    end
    return e
end



function model_parameters(nt::NamedTuple, params::Vector{Symbol}, err_string::String)
    # check that every parameter in m exists in the named tuple nt
    @assert all( i-> i in keys(nt),params) "Missing parameter(s) in " * err_string
    # check that variableno extra parameters have been provided
    @assert length(params) == length(nt) "Mismatch number of model parameters and parameters provided in " * err_string

    p=map(i->RheoFloat(nt[i]),params)
    p = convert(Array{RheoFloat,1},p)
end




function RheoModel(m::RheoModelClass, nt0::NamedTuple; log_add::Array{String} = [" "])

    # check all parameters are provided and create a well ordered named tuple
    p = model_parameters(nt0, m.params,"model definition")
    nt=NamedTuple{Tuple(m.params)}(p)
    # string printed
    info = string("\nModel: $(m.name)\n\nParameter values: $nt \n",m.info)

    # This attaches expressions to the models with parameters substituted by values.
    # Future development: this could be disabled with a key word if processing time is an issue
    G = expr_replace(m.expressions.G, nt)
    J = expr_replace(m.expressions.J, nt)
    Gp = expr_replace(m.expressions.Gp, nt)
    Gpp = expr_replace(m.expressions.Gpp, nt)

    expressions=NamedTuple{(:G,:J,:Gp,:Gpp)}( ( G, J, Gp, Gpp ) )

    @eval return( RheoModel(
    (t -> begin $G; end) |> FunctionWrapper{RheoFloat,Tuple{RheoFloat}},
    (ta -> begin [$G for t in ta]; end) |> FunctionWrapper{Array{RheoFloat,1},Tuple{Array{RheoFloat,1}}},
    (t -> begin $J; end) |> FunctionWrapper{RheoFloat,Tuple{RheoFloat}},
    (ta -> begin [$J for t in ta]; end) |> FunctionWrapper{Array{RheoFloat,1},Tuple{Array{RheoFloat,1}}},
    (t -> begin $Gp; end) |> FunctionWrapper{RheoFloat,Tuple{RheoFloat}},
    (ta -> begin [$Gp for t in ta]; end) |> FunctionWrapper{Array{RheoFloat,1},Tuple{Array{RheoFloat,1}}},
    (t -> begin $Gpp; end) |> FunctionWrapper{RheoFloat,Tuple{RheoFloat}},
    (ta -> begin [$Gpp for t in ta]; end) |> FunctionWrapper{Array{RheoFloat,1},Tuple{Array{RheoFloat,1}}},
    $expressions, $nt, $info, $log_add) )
end


#
#   This form of the constructor does not work within function, due to world age restrictions
#   Could be reintegrated later as an option for people who know what they are doing.
#   https://discourse.julialang.org/t/how-to-bypass-the-world-age-problem/7012/10
#
# function RheoModel(m::RheoModelClass, nt0::NamedTuple; log_add::Array{String} = [" "])
#
#     # check all parameters are provided and create a well ordered named tuple
#     p = model_parameters(nt0, m.params,"model definition")
#     nt=NamedTuple{Tuple(m.params)}(p)
#     # string printed
#     info = string("\nModel: $(m.name)\n\nParameter values: $nt \n",m.info)
#
#
#     # This section creates moduli functions with material parameters
#     # replaced by specific values.
#       #
#       #  Need to add unique string to function name to make sure
#       #   different instances of the same model do not overwrite the functions.
#       #
#     gs=Symbol("G_"*m.name)
#     js=Symbol("J_"*m.name)
#     gps=Symbol("Gp_"*m.name)
#     gpps=Symbol("Gpp_"*m.name)
#     ineq=Symbol("ineq_"*m.name)
#
#     G = expr_replace(m.expressions.G, nt)
#     J = expr_replace(m.expressions.J, nt)
#     Gp = expr_replace(m.expressions.Gp, nt)
#     Gpp = expr_replace(m.expressions.Gpp, nt)
#     Ineq = expr_replace(m.expressions.ineq, nt)
#
#     expressions=NamedTuple{(:G,:J,:Gp,:Gpp)}( ( G, J, Gp, Gpp ) )
#
#     @eval $gs(t::RheoFloat) = begin $G; end
#     @eval $gs(ta::Array{RheoFloat,1}) = broadcast($gs, ta)
#     @eval $gs(t::Union{Array{T1,1},T1}) where T1<:Real = $gs(rheoconv(t))
#
#     @eval $js(t::RheoFloat) = begin $J; end
#     @eval $js(ta::Array{RheoFloat,1}) = broadcast($js, ta)
#     @eval $js(t::Union{Array{T1,1},T1}) where T1<:Real = $js(rheoconv(t))
#
#     @eval $gps(ω::RheoFloat) = begin $Gp; end
#     @eval $gps(ωa::Array{RheoFloat,1}) = broadcast($gps, ωa)
#     @eval $gps(ω::Union{Array{T1,1},T1}) where T1<:Real = $gps(rheoconv(ω))
#
#     @eval $gpps(ω::RheoFloat) = begin $Gpp; end
#     @eval $gpps(ωa::Array{RheoFloat,1}) = broadcast($gpps, ωa)
#     @eval $gpps(ω::Union{Array{T1,1},T1}) where T1<:Real = $gpps(rheoconv(ω))
#
#     @eval $ineq() = begin $Ineq; end
#
#     @eval tmpmod = RheoModel($gs,$js,$gps,$gpps, $expressions, $nt, $ineq, $info, $log_add)
#     return tmpmod
#     #return RheoModel(eval(gs),eval(js),eval(gps),eval(gpps), expressions, nt, eval(ineq), info, log_add)
# end


export freeze_params

"""
    freeze_params(m::RheoModelClass, name::String, nt0::NamedTuple)

Return a new RheoModelClass with some of the parameters frozen to specific values

# Fields

- m: original RheoModelClass
- name: name of the modified model
- nt0: named tuple with values for each parameter to freeze

# Example

julia> SLS2_mod = freeze_params( SLS2, (G₀=2,η₂=3.5))
[...]

julia> SLS2.G(1,[2,1,2,3,3.5])
3.8796492f0

julia> SLS2_mod.G(1,[1,2,3])
3.8796492f0

"""
function freeze_params(m::RheoModelClass, nt0::NamedTuple)

    # check that every parameter in m exists in the named tuple nt
    @assert all( i-> i in m.params,keys(nt0)) "A parameter to freeze is not involved in the model"
    # convert values format for consistency
    nt=NamedTuple{keys(nt0)}([RheoFloat(i) for i in nt0])

    # create array of remaining variables
    p = filter(s -> !(s in keys(nt)),m.params)

    name="$(m.name) with set parameters: $nt"
    info = m.info

    # Expression to unpack parameter array into suitably names variables in the moduli expressions
    unpack_expr = Meta.parse(string(join(string.(p), ","), ",=params"))

    # This section creates moduli expressions with material parameters
    # replaced by specific values.

    G = expr_replace(m.expressions.G, nt)
    J = expr_replace(m.expressions.J, nt)
    Gp = expr_replace(m.expressions.Gp, nt)
    Gpp = expr_replace(m.expressions.Gpp, nt)
    constraint = expr_replace(m.expressions.constraint, nt)

    expressions=NamedTuple{(:G,:J,:Gp,:Gpp)}( ( G, J, Gp, Gpp ) )


    @eval return( RheoModelClass($name, $p,
        ((t,params) -> begin $unpack_expr; $G; end) |> FunctionWrapper{RheoFloat,Tuple{RheoFloat,Array{RheoFloat,1}}},
        ((ta,params) -> begin $unpack_expr; [$G for t in ta]; end) |> FunctionWrapper{Array{RheoFloat,1},Tuple{Array{RheoFloat,1},Array{RheoFloat,1}}},
        ((t,params) -> begin $unpack_expr; $J; end) |> FunctionWrapper{RheoFloat,Tuple{RheoFloat,Array{RheoFloat,1}}},
        ((ta,params) -> begin $unpack_expr; [$J for t in ta]; end) |> FunctionWrapper{Array{RheoFloat,1},Tuple{Array{RheoFloat,1},Array{RheoFloat,1}}},
        ((t,params) -> begin $unpack_expr; $Gp; end) |> FunctionWrapper{RheoFloat,Tuple{RheoFloat,Array{RheoFloat,1}}},
        ((ta,params) -> begin $unpack_expr; [$Gp for t in ta]; end) |> FunctionWrapper{Array{RheoFloat,1},Tuple{Array{RheoFloat,1},Array{RheoFloat,1}}},
        ((t,params) -> begin $unpack_expr; $Gpp; end) |> FunctionWrapper{RheoFloat,Tuple{RheoFloat,Array{RheoFloat,1}}},
        ((ta,params) -> begin $unpack_expr; [$Gpp for t in ta]; end) |> FunctionWrapper{Array{RheoFloat,1},Tuple{Array{RheoFloat,1},Array{RheoFloat,1}}},
        (params -> begin $unpack_expr; $constraint; end) |> FunctionWrapper{Bool,Tuple{Array{RheoFloat,1}}},
        $info, $expressions)   )
end




function Base.show(io::IO, m::RheoModel)
    print(io,m.info)
end



export RheoModelClass, RheoModel, model_parameters


#
#  Do we still need this function?
#
function null_modulus(t::Vector{RheoFloat}, params::Vector{T}) where T<:Real
    return [-1.0]
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

    parameters::Vector{RheoFloat}

    log::Vector{String}

end


RheologyModel(;G::Function = null_modulus,
               J::Function = null_modulus,
               Gp::Function = null_modulus,
               Gpp::Function = null_modulus,
               params::Vector{T} where T<:Real = [-1.0],
               log::Vector{String} = ["model created by user with parameters $params"]) = RheologyModel(G, J, Gp, Gpp, convert(Vector{RheoFloat},params), log)

















struct RheologyData{T<:Real}

    σ::Vector{T}
    ϵ::Vector{T}
    t::Vector{T}

    sampling::String

    log::Vector{String}

end

eltype(::RheologyData{T}) where T = T




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

function -(self::RheologyData)

    ϵ = -self.ϵ
    σ = -self.σ

    # log
    log = vcat(self.log, ["multiplied by -1"])

    RheologyData(σ, ϵ, self.t, "constant", log)

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
struct RheologyDynamic{T<:Real}

    # original data
    Gp::Vector{T}
    Gpp::Vector{T}
    ω::Vector{T}

    # operations applied, stores history of which functions (including arguments)
    log::Vector{String}

end

eltype(::RheologyDynamic{T}) where T = T

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
