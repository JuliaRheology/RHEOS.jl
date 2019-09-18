#!/usr/bin/env julia

#=
-------------------------
Log related functionality
-------------------------
=#
struct RheoLogItem
   action      # Nothing, or NamedTuple with fields:
               #     type::Symbol, funct::Symbol, params::NamedTuple, keywords::NamedTuple
   info        # NamedTuple for process output, comments, etc.
end

const RheoLog = Vector{RheoLogItem}

# different types of actions
# source: create new data (import, model simulation, etc)
# process: takes data -> new data returned  (cut, resample, filter, etc)
# analysis: data provded, returns outcome of analysis, and log changes on original data (scale, modelfit, etc.)

# process=(type= :process, funct=:test, params=(x=1,y=1) , keywords=(k=1, l=2))



#
#  Constructors
#

# simply pass text as comment
function RheoLogItem(s::String)
   return(RheoLogItem(Nothing,(comment=s,)))
end

# store info data in log from arbitrary keyword arguments
function RheoLogItem(;kwargs...)
   return(RheoLogItem(Nothing, kwargs.data))
end

#
# function Base.show(io::IO, rli::RheoLogItem)
#     println()
#     print("action = "); println(rli.action)
#     print("info   = "); println(rli.info)
# end

# function Base.display(io::IO, rl::RheoLog)
#     println("Hello!")
#     for rli in rl
#         print(rli)
#     end
# end

#=
-------------------------------
Time data related functionality
-------------------------------
=#
"""
    RheoTimeData(;σ::Vector{T1}, ϵ::Vector{T2}, t::Vector{T3}) where {T1<:Real, T2<:Real, T3<:Real}

`RheoTimeData` struct contains stress, strain and time data.

If preferred, an instance can be generated manually by just providing the three data
vectors in the right order, sampling type will be checked automatically.

# Fields

- `σ`: stress
- `ϵ`: strain
- `t`: time
- `log`: a log of struct's events, e.g. preprocessing
"""
struct RheoTimeData

    σ::Vector{RheoFloat}
    ϵ::Vector{RheoFloat}
    t::Vector{RheoFloat}

    log::Union{RheoLog,Nothing}

end

function RheoTimeData(;ϵ::Vector{T1} = RheoFloat[], σ::Vector{T2} = RheoFloat[], t::Vector{T3} = RheoFloat[], comment="Created from generic constructor", savelog = true, log = savelog ? RheoLogItem(comment) : nothing)  where {T1<:Real, T2<:Real, T3<:Real}
    typecheck = check_time_data_consistency(t,ϵ,σ)
    RheoTimeData(convert(Vector{RheoFloat},σ), convert(Vector{RheoFloat},ϵ), convert(Vector{RheoFloat},t),
    log == nothing ? nothing : [ RheoLogItem(log.action,merge(log.info, (type=typecheck,)))]     )
end

@enum TimeDataType invalid_time_data=-1 time_only=0 strain_only=1 stress_only=2 strain_and_stress=3

function check_time_data_consistency(t,e,s)
    @assert (length(t)>0)  "Time data empty"

    sdef=(s != RheoFloat[])
    edef=(e != RheoFloat[])
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

function Base.show(io::IO, d::RheoTimeData)
    b = length(d.t) > 10
    n = b ? 10 : length(d.t)

    if d.t != RheoFloat[]
        print("t =\t")
        for i = 1:n
            print("$(d.t[i])\t")
        end
        println( b ? "..." : "")
    end
    if d.ϵ != RheoFloat[]
        print("ϵ =\t")
        for i = 1:n
            print("$(d.ϵ[i])\t")
        end
        println( b ? "..." : "")
    end
    if d.σ != RheoFloat[]
        print("σ =\t")
        for i = 1:n
            print("$(d.σ[i])\t")
        end
        println( b ? "..." : "")
    end
end

function +(d1::RheoTimeData, d2::RheoTimeData)

    type1 = RheoTimeDataType(d1)
    type2 = RheoTimeDataType(d2)
    @assert (type1!=invalid_time_data) "Addition error: first parameter invalid"
    @assert (type2!=invalid_time_data) "Addition error: second parameter invalid"
    @assert (type1==type2) "Addition error: parameters inconsistent"
    @assert (type1!=time_only) "Addition error: time only data cannot be added"
    @assert (d1.t == d2.t) "Addition error: timelines inconsistent"

    log =   if (d1.log == nothing) || (d2.log == nothing)
                nothing
            else
                [ RheoLogItem( (type = :source, funct = :+, params = (rl1 = d1.log, rl2 = d2.log), keywords = ()), ()  ) ]
            end

    if (type1==strain_only) && (type2==strain_only)
        return RheoTimeData(RheoFloat[], d1.ϵ+d2.ϵ, d1.t, log)
    end

    if (type1==stress_only) && (type2==stress_only)
        return RheoTimeData(d1.σ+d2.σ, RheoFloat[], d1.t, log)
    end

    if (type1==strain_and_stress) && (type2==strain_and_stress)
        return RheoTimeData(d1.σ+d2.σ, d1.ϵ+d2.ϵ, d1.t, log)
    end

end

function +(rl1::RheoLog, rl2::RheoLog)
    return(rheologrun(rl1) + rheologrun(rl2))
end

function -(d1::RheoTimeData, d2::RheoTimeData)

    type1 = RheoTimeDataType(d1)
    type2 = RheoTimeDataType(d2)
    @assert (type1!=invalid_time_data) "Subtraction error: first parameter invalid"
    @assert (type2!=invalid_time_data) "Subtraction error: second parameter invalid"
    @assert (type1==type2) "Subtraction error: parameters inconsistent"
    @assert (type1!=time_only) "Subtraction error: time only data cannot be added"
    @assert (d1.t == d2.t) "Subtraction error: timelines inconsistent"

    log =   if (d1.log == nothing) || (d2.log == nothing)
                nothing
            else
                [ RheoLogItem( (type = :source, funct = :-, params = (rl1 = d1.log, rl2 = d2.log), keywords = ()), ()  ) ]
            end

    if (type1==strain_only) && (type2==strain_only)
        return RheoTimeData(RheoFloat[], d1.ϵ-d2.ϵ, d1.t, log)
    end

    if (type1==stress_only) && (type2==stress_only)
        return RheoTimeData(d1.σ-d2.σ, RheoFloat[], d1.t, log)
    end

    if (type1==strain_and_stress) && (type2==strain_and_stress)
        return RheoTimeData(d1.σ-d2.σ, d1.ϵ-d2.ϵ, d1.t, log)
    end

end

function -(rl1::RheoLog, rl2::RheoLog)
    return(rheologrun(rl1) - rheologrun(rl2))
end

function -(d::RheoTimeData)

    type = RheoTimeDataType(d)
    @assert (type!=invalid_time_data) "unary - error: parameter invalid"
    @assert (type!=time_only) "unary - error: time only data cannot be manipulated this way"

    log = d.log == nothing ? nothing : [d.log; RheoLogItem( (type=:process, funct=:-, params=(), keywords=() ), () ) ]

    return RheoTimeData(-d.σ, -d.ϵ, d.t, log)

end

function *(operand::Real, d::RheoTimeData)

    type = RheoTimeDataType(d)
    @assert (type!=invalid_time_data) "* error: parameter invalid"
    @assert (type!=time_only) "* error: time only data cannot be manipulated this way"

    log = d.log == nothing ? nothing : [d.log; RheoLogItem( (type=:process, funct=:*, params=(operand = operand), keywords=() ), () ) ]

    return( RheoTimeData(operand .* d.σ, operand .* d.ϵ, d.t, log) )
end

function *(d::RheoTimeData, operand::Real)
    return(operand * d)
end

#  Time shift operator >>
#  shift time by a certain amount, trash the end and pad at the start with 0

#  Union operator |
#  Combine stress_only data and strain_only data into stress_strain

#=
------------------------------------
Frequency data related functionality
------------------------------------
=#
"""
    RheoFreqData(Gp::Vector{T1}, Gpp::Vector{T2}, ω::Vector{T3}, log::OrderedDict{Any,Any}) where {T1<:Real, T2<:Real, T3<:Real}

`RheoFreqData` contains storage modulus, loss modulus and frequency data.

If preferred, an instance can be generated manually by just providing the three data
vectors in the right order.

# Fields

- `Gp`: storage modulus
- `Gpp`: loss modulus
- `ω`: frequency
- `log`: a log of struct's events, e.g. preprocessing
"""
struct RheoFreqData

    # Complex modulus data
    Gp::Vector{RheoFloat}
    Gpp::Vector{RheoFloat}
    ω::Vector{RheoFloat}

    log::Union{RheoLog,Nothing}
end

function RheoFreqData(;Gp::Vector{T1} = RheoFloat[], Gpp::Vector{T2} = RheoFloat[], ω::Vector{T3} = RheoFloat[], comment="Created from generic constructor", savelog = true, log = savelog ? RheoLogItem(comment) : nothing)  where {T1<:Real, T2<:Real, T3<:Real}
    typecheck = check_freq_data_consistency(ω,Gp,Gpp)
    RheoFreqData(convert(Vector{RheoFloat},Gp), convert(Vector{RheoFloat},Gpp), convert(Vector{RheoFloat},ω),
    log == nothing ? nothing : [ RheoLogItem(log.action,merge(log.info, (type=typecheck,)))]     )
end

@enum FreqDataType invalid_freq_data=-1 freq_only=0 with_modulus=1

function check_freq_data_consistency(o,gp,gpp)
    @assert (length(o)>0)  "Freq data empty"

    gpdef=(gp != RheoFloat[])
    gppdef=(gpp != RheoFloat[])
    if (gpdef && gppdef)
        @assert (length(gp)==length(o)) && (length(gpp)==length(o)) "Data length inconsistent"
        return with_modulus
    end

    if ((!gpdef) && (!gppdef))
        return freq_only
    end

    return invalid_freq_data

end

function RheoFreqDataType(d::RheoFreqData)
    return check_freq_data_consistency(d.ω,d.Gp,d.Gpp)
end

function Base.show(io::IO, d::RheoFreqData)
    b = length(d.ω) > 10
    n = b ? 10 : length(d.ω)

    if d.ω != RheoFloat[]
        print("ω =\t")
        for i = 1:n
            print("$(d.ω[i])\t")
        end
        println( b ? "..." : "")
    end
    if d.Gp != RheoFloat[]
        print("Gp =\t")
        for i = 1:n
            print("$(d.Gp[i])\t")
        end
        println( b ? "..." : "")
    end
    if d.Gpp != RheoFloat[]
        print("Gpp =\t")
        for i = 1:n
            print("$(d.Gpp[i])\t")
        end
        println( b ? "..." : "")
    end
end

#=
---------------------------------------
Log-based data generation functionality
---------------------------------------
=#
function rheologrun(rli::RheoLogItem, d=nothing)

      if typeof(rli.action) <: NamedTuple && :type in keys(rli.action)
         type=rli.action.type
         if type==:source && d==nothing
            return(eval(rli.action.funct)(rli.action.params...;rli.action.keywords...))
         elseif type==:process && d!=nothing
            return(eval(rli.action.funct)(d,rli.action.params...;rli.action.keywords...))
         elseif type==:analysis && d!=nothing
            return(eval(rli.action.funct)(d,rli.action.params...;rli.action.keywords...))
         end
      end
   println(rli.info)
end

function rheologrun(arli::RheoLog, d::Union{RheoTimeData,RheoFreqData,Nothing} = nothing)

  # check first item is a source item
  if d == nothing
        @assert typeof(arli[1].action) <: NamedTuple && :type in keys(arli[1].action) && arli[1].action.type == :source
                "Source missing in RheoLog"
        end


  for rli in arli
      if typeof(rli.action) <: NamedTuple && :type in keys(rli.action)
          type=rli.action.type
          if type==:analysis && d!=nothing
              rheologrun(rli, d)
          else  d=rheologrun(rli, d)
          end
      else
          println(rli.info)
      end
  end
  return(d)
end

#=
---------------------------
Model related functionality
---------------------------
=#
"""
    RheoModelClass(name::String, params::Vector{Symbol}, G::FunctionWrapper{RheoFloat,Tuple{RheoFloat,Vector{RheoFloat}}}, Ga::FunctionWrapper{Vector{RheoFloat},Tuple{Vector{RheoFloat},Vector{RheoFloat}}}, J::FunctionWrapper{RheoFloat,Tuple{RheoFloat,Vector{RheoFloat}}}, Ja::FunctionWrapper{Vector{RheoFloat},Tuple{Vector{RheoFloat},Vector{RheoFloat}}}, Gp::FunctionWrapper{RheoFloat,Tuple{RheoFloat,Vector{RheoFloat}}}, Gpa::FunctionWrapper{Vector{RheoFloat},Tuple{Vector{RheoFloat},Vector{RheoFloat}}}, Gpp::FunctionWrapper{RheoFloat,Tuple{RheoFloat,Vector{RheoFloat}}}, Gppa::FunctionWrapper{Vector{RheoFloat},Tuple{Vector{RheoFloat},Vector{RheoFloat}}}, constraint::FunctionWrapper{Bool,Tuple{Vector{RheoFloat}}}, info::String, expressions::NamedTuple)

`RheoModelClass` contains a model name, it's symbolic parameters and all its moduli (both single-input and array-input versions).

It also contains information about any constraints that must be observed (e.g. the springpot coefficient being inbetween 0 and 1).

Lastly, it also contains additional info about the model which may include a text-art schematic.

Generally, users will want to use the `RheoModelClass` constructor function as shown in the 'Create Your Model' section of the documentation
rather than the default constructor.
"""
struct RheoModelClass

    name::String
    params::Vector{Symbol}

    G::FunctionWrapper{RheoFloat,Tuple{RheoFloat,Vector{RheoFloat}}}
    Ga::FunctionWrapper{Vector{RheoFloat},Tuple{Vector{RheoFloat},Vector{RheoFloat}}}
    J::FunctionWrapper{RheoFloat,Tuple{RheoFloat,Vector{RheoFloat}}}
    Ja::FunctionWrapper{Vector{RheoFloat},Tuple{Vector{RheoFloat},Vector{RheoFloat}}}
    Gp::FunctionWrapper{RheoFloat,Tuple{RheoFloat,Vector{RheoFloat}}}
    Gpa::FunctionWrapper{Vector{RheoFloat},Tuple{Vector{RheoFloat},Vector{RheoFloat}}}
    Gpp::FunctionWrapper{RheoFloat,Tuple{RheoFloat,Vector{RheoFloat}}}
    Gppa::FunctionWrapper{Vector{RheoFloat},Tuple{Vector{RheoFloat},Vector{RheoFloat}}}

    constraint::FunctionWrapper{Bool,Tuple{Vector{RheoFloat}}}

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

rheoconv(t::Real) = RheoFloat(t)
rheoconv(t::Array{T,1}) where T<:Real = convert(Vector{RheoFloat},t)

invLaplace(f::Function, t::Array{RheoFloat}) = InverseLaplace.talbotarr(f, t)
invLaplace(f::Function, t::RheoFloat) = InverseLaplace.talbot(f, t)

# Define Null Expression, and make it default parameter
# Model name and parameters should be required --> assert
# info should have generic default value, such as "Model with $n params named $s..."

function RheoModelClass(;name::String="Custom model",
                         p::Array{Symbol}=[],
                         G::Expr = quote NaN end,
                         J::Expr = quote NaN end,
                         Gp::Expr = quote NaN end,
                         Gpp::Expr = quote NaN end,
                         constraint::Expr = quote true end,
                         info=name)

    # Expression to unpack parameter array into suitably names variables in the moduli expressions
    unpack_expr = Meta.parse(string(join(string.(p), ","), ",=params"))
    expressions = (G=G,J=J,Gp=Gp,Gpp=Gpp,constraint=constraint)

    @eval return(RheoModelClass($name, $p,
        ((t,params) -> begin $unpack_expr; $G; end)                 |> FunctionWrapper{RheoFloat,Tuple{RheoFloat,Vector{RheoFloat}}},
        ((ta,params) -> begin $unpack_expr; [$G for t in ta]; end)  |> FunctionWrapper{Vector{RheoFloat},Tuple{Vector{RheoFloat},Vector{RheoFloat}}},
        ((t,params) -> begin $unpack_expr; $J; end)                 |> FunctionWrapper{RheoFloat,Tuple{RheoFloat,Vector{RheoFloat}}},
        ((ta,params) -> begin $unpack_expr; [$J for t in ta]; end)  |> FunctionWrapper{Vector{RheoFloat},Tuple{Vector{RheoFloat},Vector{RheoFloat}}},
        ((ω,params) -> begin $unpack_expr; $Gp; end)                |> FunctionWrapper{RheoFloat,Tuple{RheoFloat,Vector{RheoFloat}}},
        ((ωa,params) -> begin $unpack_expr; [$Gp for ω in ωa]; end) |> FunctionWrapper{Vector{RheoFloat},Tuple{Vector{RheoFloat},Vector{RheoFloat}}},
        ((ω,params) -> begin $unpack_expr; $Gpp; end)               |> FunctionWrapper{RheoFloat,Tuple{RheoFloat,Vector{RheoFloat}}},
        ((ωa,params) -> begin $unpack_expr; [$Gpp for ω in ωa]; end) |> FunctionWrapper{Vector{RheoFloat},Tuple{Vector{RheoFloat},Vector{RheoFloat}}},
        (params -> begin $unpack_expr; $constraint; end)            |> FunctionWrapper{Bool,Tuple{Vector{RheoFloat}}},
        $info, $expressions) )
end

# Cool replacement function inspired from
# https://stackoverflow.com/questions/29778698/julia-lang-metaprogramming-turn-expression-into-function-with-expression-depend
# This probably could go in its own module as this is very useful.
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
    @assert all(i->i in keys(nt),params) "Missing parameter(s) in " * err_string
    # check that no extra parameters have been provided
    @assert length(params) == length(nt) "Mismatch number of model parameters and parameters provided in " * err_string

    p = map(i->RheoFloat(nt[i]), params)
end

"""
    freeze_params(m::RheoModelClass, nt0::NamedTuple)

Return a new `RheoModelClass` with some of the parameters frozen to specific values

# Arguments

- `m`: original `RheoModelClass`
- `nt0`: named tuple with values for each parameter to freeze

# Example
```@example
julia> SLS2_mod = freeze_params( SLS2, (G₀=2,η₂=3.5))
[...]

julia> SLS2.G(1,[2,1,2,3,3.5])
3.8796492

julia> SLS2_mod.G(1,[1,2,3])
3.8796492
```

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
        ((t,params) -> begin $unpack_expr; $G; end) |> FunctionWrapper{RheoFloat,Tuple{RheoFloat,Vector{RheoFloat}}},
        ((ta,params) -> begin $unpack_expr; [$G for t in ta]; end) |> FunctionWrapper{Vector{RheoFloat},Tuple{Vector{RheoFloat},Vector{RheoFloat}}},
        ((t,params) -> begin $unpack_expr; $J; end) |> FunctionWrapper{RheoFloat,Tuple{RheoFloat,Vector{RheoFloat}}},
        ((ta,params) -> begin $unpack_expr; [$J for t in ta]; end) |> FunctionWrapper{Vector{RheoFloat},Tuple{Vector{RheoFloat},Vector{RheoFloat}}},
        ((ω,params) -> begin $unpack_expr; $Gp; end) |> FunctionWrapper{RheoFloat,Tuple{RheoFloat,Vector{RheoFloat}}},
        ((ωa,params) -> begin $unpack_expr; [$Gp for ω in ωa]; end) |> FunctionWrapper{Vector{RheoFloat},Tuple{Vector{RheoFloat},Vector{RheoFloat}}},
        ((ω,params) -> begin $unpack_expr; $Gpp; end) |> FunctionWrapper{RheoFloat,Tuple{RheoFloat,Vector{RheoFloat}}},
        ((ωa,params) -> begin $unpack_expr; [$Gpp for ω in ωa]; end) |> FunctionWrapper{Vector{RheoFloat},Tuple{Vector{RheoFloat},Vector{RheoFloat}}},
        (params -> begin $unpack_expr; $constraint; end) |> FunctionWrapper{Bool,Tuple{Vector{RheoFloat}}},
        $info, $expressions)   )
end


"""
    RheoModel(G::FunctionWrapper{RheoFloat,Tuple{RheoFloat}}, Ga::FunctionWrapper{Vector{RheoFloat},Tuple{Vector{RheoFloat}}}, J::FunctionWrapper{RheoFloat,Tuple{RheoFloat}}, Ja::FunctionWrapper{Vector{RheoFloat},Tuple{Vector{RheoFloat}}}, Gp::FunctionWrapper{RheoFloat,Tuple{RheoFloat}}, Gpa::FunctionWrapper{Vector{RheoFloat},Tuple{Vector{RheoFloat}}}, Gpp::FunctionWrapper{RheoFloat,Tuple{RheoFloat}}, Gppa::FunctionWrapper{Vector{RheoFloat},Tuple{Vector{RheoFloat}}}, expressions::NamedTuple)
    
`RheoModel` contains all known moduli of a particular model, as for a `RheoModelClass` model name. However, a `RheoModel`
has all it's parameters fixed to known values.

Generally, users will begin with a defined `RheoModelClass` (e.g. SLS) and then specialise the parameters,
rather than calling the default constructor explicitly.
"""
struct RheoModel

    G::FunctionWrapper{RheoFloat,Tuple{RheoFloat}}
    Ga::FunctionWrapper{Vector{RheoFloat},Tuple{Vector{RheoFloat}}}
    J::FunctionWrapper{RheoFloat,Tuple{RheoFloat}}
    Ja::FunctionWrapper{Vector{RheoFloat},Tuple{Vector{RheoFloat}}}
    Gp::FunctionWrapper{RheoFloat,Tuple{RheoFloat}}
    Gpa::FunctionWrapper{Vector{RheoFloat},Tuple{Vector{RheoFloat}}}
    Gpp::FunctionWrapper{RheoFloat,Tuple{RheoFloat}}
    Gppa::FunctionWrapper{Vector{RheoFloat},Tuple{Vector{RheoFloat}}}
    expressions::NamedTuple


    params::NamedTuple
    info::String    # do we need this?
    # log::OrderedDict{Any,Any}

end

function RheoModel(m::RheoModelClass; kwargs...)
    return(RheoModel(m,kwargs.data))
end

function RheoModel(m::RheoModelClass, nt0::NamedTuple)

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
    (ta -> begin [$G for t in ta]; end) |> FunctionWrapper{Vector{RheoFloat},Tuple{Vector{RheoFloat}}},
    (t -> begin $J; end) |> FunctionWrapper{RheoFloat,Tuple{RheoFloat}},
    (ta -> begin [$J for t in ta]; end) |> FunctionWrapper{Vector{RheoFloat},Tuple{Vector{RheoFloat}}},
    (ω -> begin $Gp; end) |> FunctionWrapper{RheoFloat,Tuple{RheoFloat}},
    (ωa -> begin [$Gp for ω in ωa]; end) |> FunctionWrapper{Vector{RheoFloat},Tuple{Vector{RheoFloat}}},
    (ω -> begin $Gpp; end) |> FunctionWrapper{RheoFloat,Tuple{RheoFloat}},
    (ωa -> begin [$Gpp for ω in ωa]; end) |> FunctionWrapper{Vector{RheoFloat},Tuple{Vector{RheoFloat}}},
    $expressions, $nt, $info) )
end

function Base.show(io::IO, m::RheoModel)
    print(io,m.info)
end
