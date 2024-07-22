#!/usr/bin/env julia

#===================================================
#
#  rheodata.jl
#
===================================================#



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

# simply pass text as comment
function RheoLogItem(s::String)
   return(RheoLogItem(Nothing,(comment=s,)))
end

# store info data in log from arbitrary keyword arguments
function RheoLogItem(;kwargs...)
   return(RheoLogItem(Nothing, kwargs.data))
end



function Base.show(io::IO, ::MIME"text/plain", rli::RheoLogItem)
    print("action = "); println(rli.action)
    print("info   = "); println(rli.info)
end


function Base.show(io::IO, ::MIME"text/plain", rl::RheoLog)
    println("RheoLog structure")
    i=0
    for rli in rl
       i+=1
       println("["*string(i)*"]")
       display(rli)
    end
end


#
# Internal functions to help manage logs
#
# not exported to users
#

function loginit(savelog, funct::Symbol; params=NamedTuple(), keywords=NamedTuple(), comment="Process added", info=(comment=comment,))
    if savelog
        [RheoLogItem( (type=:source, funct=funct, params=params, keywords=keywords), info )]
    else
        nothing
    end
end


function logadd_process(d, funct::Symbol; params=NamedTuple(), keywords=NamedTuple(), comment="Process added", info=(comment=comment,))
    d.log === nothing ? nothing : [d.log; RheoLogItem( (type=:process, funct=funct, params=params, keywords=keywords), info ) ]
end

function logadd_process!(d, funct::Symbol; params=NamedTuple(), keywords=NamedTuple(), comment="Process added", info=(comment=comment,))
    if d.log !== nothing
        push!(d.log, RheoLogItem( (type=:process, funct=funct, params=params, keywords=keywords), info) )
    end
end


function logadd_analysis!(d, funct::Symbol; params=NamedTuple(), keywords=NamedTuple(), comment="Process added", info=(comment=comment,))
    if d.log !== nothing
        push!(d.log, RheoLogItem( (type=:analysis, funct=funct, params=params, keywords=keywords), info) )
    end
end






# ----------
# Utilities
# ----------


"""
    rheoconvert(t)

Converts if required the scalar value or vector of real number `t` to the type `RheoFloat`.
"""
rheoconvert(t::Real) = RheoFloat(t)
rheoconvert(t::RheoFloat) = t
rheoconvert(t::Vector{T}) where T<:Real = convert(Vector{RheoFloat},t)
rheoconvert(t::Vector{RheoFloat}) = t



#
# Utilities for in-place functions
#
# Internal - Not be exposed to the user.
#

"""
   _setdata!(a::Vector{RheoFloat}, t::Vector{T}) where {T<:Real}
 
In place replacement of the values of a vector `a` with the values of another one.
"""
function _setdata!(a::Vector{RheoFloat}, t::Vector{T}) where {T<:Real}
    if length(a) == 0
        append!(a,t)
    elseif length(a) == length(t)
        a .= t
    else
        empty!(a)
        append!(a,t)
    end
end

"""
    _mapdata!(f::Function, a::Vector{RheoFloat}, t::Vector{T} where {T<:Real})

In place application of a function to the elements of a vector.
Useful for data generation through strainfunctions and related.
"""
function _mapdata!(f::Function, a::Vector{RheoFloat}, t::Vector{T} where {T<:Real})
    if length(a) != length(t)
        resize!(a,length(t))
    end
    map!(f,a,t)
end









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

function RheoTimeData(;strain = RheoFloat[], ϵ::Vector{T1} = strain, stress = RheoFloat[], σ::Vector{T2} = stress, t::Vector{T3} = RheoFloat[], 
                      comment="Created from generic constructor", savelog = true, log = nothing)  where {T1<:Real, T2<:Real, T3<:Real}
    typecheck = check_time_data_consistency(t,ϵ,σ)
         
    RheoTimeData(rheoconvert(σ), rheoconvert(ϵ), rheoconvert(t),
                 isnothing(log) ? loginit(savelog, :RheoTimeData, params = NamedTuple(), keywords = (ϵ=ϵ,σ=σ,t=t,comment=comment),
                                           info=(comment=comment, type=typecheck) )  
                                : log ) 
end




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




@enum TimeDataType invalid=-1 time_only=0 strain_only=1 stress_only=2 strain_and_stress=3



@enum RheosExceptionType timedata_invalid freqdata_invalid uniform_sampling_required model_parameters_invalid

struct RheosException <: Exception
    type::RheosExceptionType
    var::String
end

function Base.showerror(io::IO, err::RheosException)
    if err.type==timedata_invalid
        print(io, "Rheos Exception - Invalid Time data: $(err.var)")
    elseif err.type==freqdata_invalid 
        print(io, "Rheos Exception - Invalid Frequency data: $(err.var)")
    else
        print(io, "Rheos Exception: $(err.var)")
    end
       
end

# throw(RheosException("this code"))

function check_time_data_consistency(t,e,s)
    # @assert (length(t)>0)  "Time data empty"
    length(t)>0 || throw(RheosException(timedata_invalid,"Time data empty"))

    sdef=(s != RheoFloat[])
    edef=(e != RheoFloat[])
    if (sdef && edef)
        (length(s)==length(t)) && (length(e)==length(t))  || throw(RheosException("Time data length inconsistent"))
        return strain_and_stress
    end

    if (sdef && (!edef))
        (length(s)==length(t))  || throw(RheosException(timedata_invalid,"Time data length inconsistent"))
        return stress_only
    end

    if (edef && (!sdef))
        (length(e)==length(t))  || throw(RheosException(timedata_invalid,"Time data length inconsistent"))
        return strain_only
    end

    if ((!edef) && (!sdef))
        return time_only
    end

    return invalid

end

"""
    rheotimedatatype(d::RheoTimeData)

check the validity of the data and return information about the type of data (time only, with strain, with stress, or both.
"""
function rheotimedatatype(d::RheoTimeData)
    return check_time_data_consistency(d.t,d.ϵ,d.σ)
end

"""
    hastime(d::RheoTimeData)

returns 'true' if d contains a time array
"""
function hastime(d::RheoTimeData)
    return (d.t != RheoFloat[])
end

"""
    hasstress(d::RheoTimeData)

returns 'true' if d contains stress data
"""
function hasstress(d::RheoTimeData)
    return (d.σ != RheoFloat[])
end

"""
    hasstrain(d::RheoTimeData)

returns 'true' if d contains strain data
"""
function hasstrain(d::RheoTimeData)
    return (d.ϵ != RheoFloat[])
end


"""
    gettime(d::RheoTimeData)

returns the time vector.
"""
function gettime(d::RheoTimeData)
    return (d.t)
end



"""
    getstrain(d::RheoTimeData)

returns the strain vector if d contains strain data.
"""
function getstrain(d::RheoTimeData)
    return (d.ϵ)
end



"""
    getstress(d::RheoTimeData)

returns the stress vector if d contains stress data.
"""
function getstress(d::RheoTimeData)
    return (d.σ)
end





# Constants that are useful for certain processing function to determine which module to use.
@enum LoadingType strain_imposed=1 stress_imposed=2








function +(d1::RheoTimeData, d2::RheoTimeData)

    type1 = rheotimedatatype(d1)
    type2 = rheotimedatatype(d2)
    @assert (type1!=invalid) "Addition error: first parameter invalid"
    @assert (type2!=invalid) "Addition error: second parameter invalid"
    @assert (type1==type2) "Addition error: parameters inconsistent"
    @assert (type1!=time_only) "Addition error: time only data cannot be added"
    @assert (d1.t == d2.t) "Addition error: timelines inconsistent"

    log =   if (d1.log === nothing) || (d2.log === nothing)
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


function -(d1::RheoTimeData, d2::RheoTimeData)

    type1 = rheotimedatatype(d1)
    type2 = rheotimedatatype(d2)
    @assert (type1!=invalid) "Subtraction error: first parameter invalid"
    @assert (type2!=invalid) "Subtraction error: second parameter invalid"
    @assert (type1==type2) "Subtraction error: parameters inconsistent"
    @assert (type1!=time_only) "Subtraction error: time only data cannot be added"
    @assert (d1.t == d2.t) "Subtraction error: timelines inconsistent"

    log =   if (d1.log === nothing) || (d2.log === nothing)
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


function -(d::RheoTimeData)

    type = rheotimedatatype(d)
    @assert (type!=invalid) "unary - error: parameter invalid"
    @assert (type!=time_only) "unary - error: time only data cannot be manipulated this way"

    log = d.log === nothing ? nothing : [d.log; RheoLogItem( (type=:process, funct=:-, params=(), keywords=() ), () ) ]

    return RheoTimeData(-d.σ, -d.ϵ, d.t, log)

end


# Required to execute +/- operators for logs

function +(rl1::RheoLog, rl2::RheoLog)
    return(rheologrun(rl1) + rheologrun(rl2))
end

function -(rl1::RheoLog, rl2::RheoLog)
    return(rheologrun(rl1) - rheologrun(rl2))
end


function *(operand::Real, d::RheoTimeData)

    type = rheotimedatatype(d)
    @assert (type!=invalid) "* error: parameter invalid"
    @assert (type!=time_only) "* error: time only data cannot be manipulated this way"

    log = d.log === nothing ? nothing : [d.log; RheoLogItem( (type=:process, funct=:*, params=(operand = operand), keywords=() ), () ) ]

    return( RheoTimeData(operand .* d.σ, operand .* d.ϵ, d.t, log) )
end

function *(d::RheoTimeData, operand::Real)
    return(operand * d)
end

function |(d1::RheoTimeData, d2::RheoTimeData)
    # get types and type-check
    type1 = rheotimedatatype(d1)
    type2 = rheotimedatatype(d2)

    @assert ((type1==strain_only && type2==stress_only) || (type1==stress_only && type2==strain_only)) "One RheoTimeData must be stress_only and the other must be strain_only"

    # call recursively if not expected way round
    (type1==strain_only && type2==stress_only) && return d2|d1

    # expected case, type1==stress_only, type2==strain_only
    @assert d1.t==d2.t "Data to be combined must have same time-series"

    log = ((d1.log === nothing) || (d2.log === nothing)) ? nothing : [RheoLogItem((type = :source, funct = :|, params = (rl1 = d1.log, rl2 = d2.log), keywords = ()), ())]

    return RheoTimeData(d1.σ, d2.ϵ, d1.t, log)
end





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
    Gp::Vector{RheoFloat}
    Gpp::Vector{RheoFloat}
    ω::Vector{RheoFloat}

    log::Union{RheoLog,Nothing}
end

function RheoFreqData(;Gp::Vector{T1} = RheoFloat[], Gpp::Vector{T2} = RheoFloat[], omega = RheoFloat[], ω::Vector{T3} = omega, 
                       comment="Created from generic constructor", savelog = true, log = nothing)  where {T1<:Real, T2<:Real, T3<:Real}

    typecheck = check_freq_data_consistency(ω,Gp,Gpp)
    RheoFreqData(rheoconvert(Gp), rheoconvert(Gpp), rheoconvert(ω),
                 isnothing(log) ? loginit(savelog, :RheoFreqData, params = NamedTuple(), keywords = (ω=ω,Gp=Gp,Gpp=Gpp,comment=comment),
                                           info=(comment=comment, type=typecheck) )  
                                : log ) 
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

function rheofreqdatatype(d::RheoFreqData)
    return check_freq_data_consistency(d.ω,d.Gp,d.Gpp)
end




"""
    hasfreq(d::RheoFreqData)

returns `true` if d contains a frequency array
"""
function hasfreq(d::RheoFreqData)
    return (d.ω != RheoFloat[])
end

"""
    hasmodulus(d::RheoFreqData)

returns `true` if d contains a complex modulus function
"""
function hasmodulus(d::RheoFreqData)
    (d.Gp != RheoFloat[]) && (d.Gpp != RheoFloat[])
end



"""
    getfreq(d::RheoFreqData)

returns the frequency ω vector.
"""
function getfreq(d::RheoFreqData)
    return (d.ω)
end


"""
    getstorage(d::RheoFreqData)

returns the storage modulus vector.
"""
function getstorage(d::RheoFreqData)
    return (d.Gp)
end



"""
    getloss(d::RheoFreqData)

returns the loss modulus vector.
"""
function getloss(d::RheoFreqData)
    return (d.Gpp)
end












#=
---------------------------------------
Log-based data generation functionality
---------------------------------------
=#
function rheologrun(rli::RheoLogItem, d=nothing)

      if typeof(rli.action) <: NamedTuple && :type in keys(rli.action)
         type=rli.action.type
         if type==:source && d===nothing
            return(eval(rli.action.funct)(rli.action.params...;rli.action.keywords...))
         elseif type==:process && d!==nothing
            return(eval(rli.action.funct)(d,rli.action.params...;rli.action.keywords...))
         elseif type==:analysis && d!==nothing
            return(eval(rli.action.funct)(d,rli.action.params...;rli.action.keywords...))
         end
      end
   println(rli.info)
end

"""
    rheologrun(log::RheoLog, d::Union{RheoTimeData,RheoFreqData,Nothing} = nothing)

execute all actions from the log. It applies them to the data `d` provided, or use the log's first action to reload/recreate it otherwise.
"""
function rheologrun(arli::RheoLog, d::Union{RheoTimeData,RheoFreqData,Nothing} = nothing)

  # check first item is a source item
  if d === nothing
        @assert typeof(arli[1].action) <: NamedTuple && :type in keys(arli[1].action) && arli[1].action.type == :source
                "Source missing in RheoLog"
        end


  for rli in arli
      if typeof(rli.action) <: NamedTuple && :type in keys(rli.action)
          type=rli.action.type
          if type==:analysis && d!==nothing
              rheologrun(rli, d)
          else  d=rheologrun(rli, d)
          end
      else
          println(rli.info)
      end
  end
  return(d)
end



"""
    showlog(d::Union{RheoTimeData,RheoFreqData})

shows the record of operations (`RheoLog`) on rheological data.
"""
function showlog(d::Union{RheoTimeData,RheoFreqData})
    if d.log !== nothing
        display(d.log)
    else
        println("No log data available")
    end
end



"""
    extractfitdata(data)
Parse a log of actions from a `RheoTimeData` instance and extract model fitting entries into a dictionary.
# Arguments
- `data`: A `RheoTimeData` instance containing the log of actions.
# Returns
A dictionary where keys are model names and values are lists of named tuples,
each containing model parameters, error, info, and index from the log.
# Example
```julia
data.log = [...]  # Define your log data
models_data = extractfitdata(data)
println(models_data)
"""
function extractfitdata(data)
    models_dict = Dict{String, Vector{NamedTuple{(:params, :info, :index), Tuple{Any, Any, Int}}}}()

    for (idx, entry) in enumerate(log)
        if entry.action.funct == :modelfit
            model_name = entry.info.model_name
            model_params = entry.info.model_params
            error = entry.info.error
            info = entry.info
            
            model_entry = (params = model_params, info = info, index = idx)

            if haskey(models_dict, model_name)
                push!(models_dict[model_name], model_entry)
            else
                models_dict[model_name] = [model_entry]
            end
        end
    end

    return models_dict
end








