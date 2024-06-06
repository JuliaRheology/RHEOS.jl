#!/usr/bin/env julia

#===================================================
#
#  rheomodel.jl
#
===================================================#








#=
---------------------------
Model related functionality
---------------------------
=#

# Function types for the moduli with free params
const FWScaFree = FunctionWrapper{RheoFloat,Tuple{RheoFloat,Vector{RheoFloat}}}
const FWVecFree = FunctionWrapper{Vector{RheoFloat},Tuple{Vector{RheoFloat},Vector{RheoFloat}}}
# Function types for the moduli with fixed params
const FWScaFixed = FunctionWrapper{RheoFloat,Tuple{RheoFloat}}
const FWVecFixed = FunctionWrapper{Vector{RheoFloat},Tuple{Vector{RheoFloat}}}
# Function types for the constraint function
const FWConstraint = FunctionWrapper{Bool,Tuple{Vector{RheoFloat}}}



# constant NaN functions for default parameters and testing existence of moduli functions
const fwscafreeNaN = ((t,p)->NaN) |> FWScaFree
const fwvecfreeNaN = ((t,p)->fill(RheoFloat(NaN), length(t))) |> FWVecFree
const fwscafixedNaN = (t->NaN) |> FWScaFixed
const fwvecfixedNaN = (t->fill(RheoFloat(NaN), length(t))) |> FWVecFixed



struct _RheoModel{TSca,TVec}

	name::String
	freeparams::Tuple{Vararg{Symbol, N}} where N
	fixedparams::NamedTuple

	# Moduli functions specialised for scalar and array values
	_G::TSca
	_Ga::TVec
	_J::TSca
	_Ja::TVec
	_Gp::TSca
	_Gpa::TVec
	_Gpp::TSca
	_Gppa::TVec

	_constraint::FWConstraint
	_Gramp::Bool

	info::String
	expressions::NamedTuple

end


#
#  Specialised forms of _RheoModel
#  Dediated constructors are defined below
#


const RheoModelClass = _RheoModel{FWScaFree,FWVecFree}

const RheoModel = _RheoModel{FWScaFixed,FWVecFixed}


function Base.show(io::IO, m::_RheoModel)
    print(io, "\nModel name: $(m.name)")
    ps=join([string(s) for s in m.freeparams], ", ", " and ")
    print(io, "\n\nFree parameters: $ps\n")
    if length(m.fixedparams)>0
        print(io, "\n\nFixed parameters: $(m.fixedparams)\n")
    end
    print(io, m.info)
    return
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

function expr_replace(ex, s, v)
    expr_replace!(copy(ex), s, v)
end

function expr_replace(ex::Expr, nt::NamedTuple)  
    e = copy(ex)
    k=keys(nt)
    v=values(nt)
    for i in 1:length(nt)
        expr_replace!(e, k[i], v[i])
    end
    return e
end

function expr_replace(ex::Nothing, nt::NamedTuple)
	return nothing
end



#
#  Look for a specific symbol in an expression
#
expr_search(ex, s) = false
expr_search(ex::Symbol, s) = s == ex ? true : false

@inline function expr_search(ex::Expr, s)
        any( expr_search(e, s) for e in ex.args )
end


@inline function _replace_symbols_with_array(e::Expr, psymbs::Tuple)
    return( expr_replace(e,NamedTuple{psymbs}([Meta.parse(":(p_arr[$i])").args[1] for i in 1:length(psymbs)]) ) )
end

@inline function _expr_cleanup!(ex::Expr)
    filter!(e -> typeof(e)!=LineNumberNode,  ex.args)
    return()
end

@inline function _expr_cleanup(ex::Expr)
    e = copy(ex)
    _expr_cleanup!(e)
    return(e)
end

modulusexists(modulus::FWScaFree) =  (modulus !== fwscafreeNaN)
modulusexists(modulus::FWVecFree) =  (modulus !== fwvecfreeNaN)
modulusexists(modulus::FWScaFixed) =  (modulus !== fwscafixedNaN)
modulusexists(modulus::FWVecFixed) =  (modulus !== fwvecfixedNaN)

function modulusexists(modulus_single_input)
   # 5.0 is not meaningful, just an input to check that the function does not return NaN
    # at a non-zero value, which is the default function output for moduli that have
    # not been defined.
    return !isnan(modulus_single_input(5.0))
end



#
#  Create the scalar and vector function wrappers for the moduli.
#

function _buildmoduli_t(G::Expr, psymbs::Tuple)
    # Replace symbols by elements in array
	Ge = _replace_symbols_with_array(G, psymbs)
    _expr_cleanup!(Ge)
    isexprmultiline = (length(Ge.args) > 1)

    # Detect if expression depends on variable or if it is a constant
    # This is because the @. macro fails to vectorize if there is no vector to iterate on
    # It also fails on multiline expressions
    islaplace, isconstant = expr_search(Ge, :t) ? (false, false) : expr_search(Ge, :s) ? (true, false) : (false, true)

    if isconstant  # constant value returned
        @eval return( (     ((t,p_arr) -> $Ge)      |> FWScaFree,
                            ((t,p_arr) -> fill(rheoconvert($Ge), length(t)))  |> FWVecFree ) )
#    elseif islaplace
#        return()
    elseif isexprmultiline
	    @eval return( (     ((t,p_arr) -> $Ge)      |> FWScaFree,
		    			    ((ta,p_arr) -> [$Ge for t in ta] )  |> FWVecFree ) )
    else
	    @eval return( (     ((t,p_arr) -> $Ge)      |> FWScaFree,
		    			    ((t,p_arr) -> @. $Ge )  |> FWVecFree ) )
    end
end




function _buildmoduli_t(G::Expr)
    Ge = _expr_cleanup(G)
    isexprmultiline = (length(Ge.args) > 1)

    islaplace, isconstant = expr_search(Ge, :t) ? (false, false) : expr_search(Ge, :s) ? (true, false) : (false, true)

    if isconstant     # constant value returned
        @eval return( ( (t -> $Ge)                             |> FWScaFixed,
                        (t -> fill(rheoconvert($Ge), length(t)))  |> FWVecFixed ) ) 
#    elseif islaplace
#        return()
    elseif isexprmultiline
        @eval return( ( (t -> $Ge)          |> FWScaFixed,
                (ta -> [$Ge for t in ta] )  	|> FWVecFixed ) )
    else
        @eval return( ( (t -> $Ge)          |> FWScaFixed,
                        (t -> @. $Ge )  	|> FWVecFixed ) ) 
    end                          
end
                


function _buildmoduli_ω(G::Expr, psymbs::Tuple)
    # Replace symbols by elements in array
	Ge = _replace_symbols_with_array(G, psymbs)
    _expr_cleanup!(Ge)
    isexprmultiline = (length(Ge.args) > 1)

    # Detect if expression depends on variable or if it is a constant
    # This is because the @. macro fails to vectorize if there is no vector to iterate on
    # It also fails on multiline expressions
    isconstant = !expr_search(Ge, :ω) 

    if isconstant  # constant value returned
        @eval return( (     ((ω,p_arr) -> $Ge)      |> FWScaFree,
                            ((ω,p_arr) -> fill(rheoconvert($Ge), length(ω)))  |> FWVecFree ) )
    elseif isexprmultiline
	    @eval return( (     ((ω,p_arr) -> $Ge)      |> FWScaFree,
		    			    ((ωa,p_arr) -> [$Ge for ω in ωa] )  |> FWVecFree ) )
    else
	    @eval return( (     ((ω,p_arr) -> $Ge)      |> FWScaFree,
		    			    ((ω,p_arr) -> @. $Ge )  |> FWVecFree ) )
    end
end



function _buildmoduli_ω(G::Expr)
    Ge = _expr_cleanup(G)
    isexprmultiline = (length(Ge.args) > 1)

    isconstant = !expr_search(Ge, :ω)

    if isconstant     # constant value returned
        @eval return( ( (ω -> $Ge)                             |> FWScaFixed,
                        (ω -> fill(rheoconvert($Ge), length(ω)))  |> FWVecFixed ) ) 
#    elseif islaplace
#        return()
    elseif true # isexprmultiline
        @eval return( ( (ω -> $Ge)          |> FWScaFixed,
                (ωa -> [$Ge for ω in ωa] )  	|> FWVecFixed ) )
    else
        @eval return( ( (ω -> $Ge)          |> FWScaFixed,
                        (ω -> @. $Ge )  	|> FWVecFixed ) ) 
    end                          
end
         

function _buildmoduli_t(G::Nothing, psymbs::Tuple)
	return( (       fwscafreeNaN, fwvecfreeNaN ) )
end

function _buildmoduli_t(G::Nothing)
	return( (       fwscafixedNaN, fwvecfixedNaN ) )
end

function _buildmoduli_ω(G::Nothing, psymbs::Tuple)
	return( (       fwscafreeNaN, fwvecfreeNaN ) )
end

function _buildmoduli_ω(G::Nothing)
	return( (       fwscafixedNaN, fwvecfixedNaN ) )
end



function _buildconstraint(constraint::Expr, psymbs::Tuple)
    # Replace symbols by elements in array
	constraint_e = expr_replace(constraint,NamedTuple{psymbs}([Meta.parse(":(p_arr[$i])").args[1] for i in 1:length(psymbs)]) )
	@eval return( (     (p_arr -> $constraint_e)      |> FWConstraint ) )
end

function _buildconstraint(constraint::Expr)
	@eval return( (     (p_arr -> $constraint)      |> FWConstraint ) )
end


"""
    RheoModelClass(;name::String, p::Tuple, G,J,Gp,Gpp, constraint, info, G_ramp)
 
`RheoModelClass` is a complex data structure that contains all the relevant data to 
allow RHEOS to fit models and make predictions once parameters are set.
RHEOS possesses already a large number of such data structures to represent common rheological models,
such as Maxwell or KelvinVoigt.

Users are not expected to directly manipulate the content of a RheoModelClass object, but will 
pass it to relevant functions for fitting and numerically evaluating visco-elastic moduli.

The model name and its parameters can however be printed and display on the REPL.

# Example
```@example
julia> Maxwell

Model name: maxwell

Free parameters: η and k

                ___
            _____| |________╱╲  ╱╲  ╱╲  ___
                _|_|          ╲╱  ╲╱  ╲╱
                  η                  k
               
```
"""
function RheoModelClass(;name::String,
        p::Tuple,
        G = nothing,
        J = nothing,
        Gp = nothing,
        Gpp = nothing,
        constraint::Expr = quote true end,
        info="", 
        # flag to indicate use of integral forms of the relaxation modulus.
               # using the ramp response instead of the step response helps avoid singularities
        G_ramp::Bool = false,
        )

    # Building expressions tuple to store data provided to constructor
    expressions = (G=G,J=J,Gp=Gp,Gpp=Gpp,constraint=constraint)

    return(RheoModelClass(name, p, NamedTuple{}(),
        _buildmoduli_t(G,p)..., _buildmoduli_t(J,p)...,
        _buildmoduli_ω(Gp,p)..., _buildmoduli_ω(Gpp,p)...,
        _buildconstraint(constraint,p),
        G_ramp, info, expressions) )
end


#
#  Function used to check and convert any named tupple of parameters into an array ordered as expected by the moduli functions of the RheoModelClass
#

function check_and_reorder_parameters(nt::NamedTuple, params_ordered_list::Tuple; err_string::String="calling function")
	# check that every parameter in m exists in the named tuple nt
	@assert all(i->i in keys(nt),params_ordered_list) "Missing parameter(s) in " * err_string
	# check that no extra parameters have been provided
	@assert length(params_ordered_list) == length(nt) "Mismatch number of model parameters and parameters provided in " * err_string

	p = map(i->RheoFloat(nt[i]), params_ordered_list)
    rheoconvert([p...])
end




function _freeze_params(m::RheoModelClass, nt0::NamedTuple)

	# check that every parameter to fix actually exists in the model
	@assert all( i-> i in m.freeparams,keys(nt0)) "A parameter to freeze is not involved in the model"
	# convert values format for consistency
	nt=NamedTuple{keys(nt0)}([RheoFloat(i) for i in nt0])

	# update list of fixed parameters
	fixedparams = merge(m.fixedparams, nt)

	# create array of remaining variables
	freeparams = Tuple( filter(s -> !(s in keys(nt)),m.freeparams) )

	# This section creates moduli expressions with material parameters
	# replaced by specific values.

	G = expr_replace(m.expressions.G, nt)
	J = expr_replace(m.expressions.J, nt)
	Gp = expr_replace(m.expressions.Gp, nt)
	Gpp = expr_replace(m.expressions.Gpp, nt)
	constraint = expr_replace(m.expressions.constraint, nt)

	return( (freeparams=freeparams, fixedparams=fixedparams, G=G,J=J,Gp=Gp,Gpp=Gpp,constraint=constraint) )

end

"""
    freezeparams(m::RheoModelClass, nt0::NamedTuple)

Return a new `RheoModelClass` with some of the parameters frozen to specific values

# Arguments

- `m`: original `RheoModelClass`
- `nt0`: named tuple with values for each parameter to freeze

# Example
```@example
julia> SLS2_mod = freezeparams( SLS2, (G₀=2,η₂=3.5))
[...]

julia> SLS2.G(1,[2,1,2,3,3.5])
3.8796492

julia> SLS2_mod.G(1,[1,2,3])
3.8796492
```

"""
function freezeparams(m::RheoModelClass, nt0::NamedTuple)

	freeparams, fixedparams, G,J,Gp,Gpp,constraint = _freeze_params(m, nt0)

	# Check that some free params remains
	@assert length(freeparams) > 0  "All parameters are set. Build a RheoModel instead."

	# Building expressions tuple to store data provided to constructor
	expressions = (G=G,J=J,Gp=Gp,Gpp=Gpp,constraint=constraint)

	return(RheoModelClass(m.name, freeparams, fixedparams,
		_buildmoduli_t(G,freeparams)..., _buildmoduli_t(J,freeparams)...,
		_buildmoduli_ω(Gp,freeparams)..., _buildmoduli_ω(Gpp,freeparams)...,
		_buildconstraint(constraint,freeparams),
        m._Gramp, m.info, expressions) )
end


function freezeparams(m::RheoModelClass; kwargs...)
  return(freezeparams(m,symbol_to_unicode(values(kwargs))))
end


function freeze_params(m::RheoModelClass, nt0::NamedTuple)
    println("freeze_params is deprecated - use freezeparams instead.")
    return(freezeparams(m,nt0))
end

function freeze_params(m::RheoModelClass; kwargs...)
    println("freeze_params is deprecated - use freezeparams instead.")
    return(freezeparams(m,symbol_to_unicode(values(kwargs))))
end



"""
    RheoModel(m::RheoModelClass, nt0::NamedTuple)

`RheoModel` represents a rheological model with set parameters. They are obtained by fitting a model to data using `modelfit`, 
or by specialising the relevant RheoModelClass by prescribing its parameters. Parameters can be provided as a named tuple or keyword arguments.

`RheoModel` objects can then be used to simulate the response to an arbitrary input using `modelpredict`, and access the values of the moduli functions.

# Example
```@example
julia> model = RheoModel(Maxwell, (k=1, η=2.))
[...]

julia> model = RheoModel(Maxwell, k=1, η=2.)
[...]
```
"""
function RheoModel(m::RheoModelClass, nt0::NamedTuple)

	freeparams, fixedparams, G,J,Gp,Gpp,constraint = _freeze_params(m, nt0)
	
	# Check all free params are set
	@assert length(freeparams) == 0  "Some parameters need to be set: $freeparams"

	# Building expressions tuple to store data provided to constructor
	expressions = (G=G,J=J,Gp=Gp,Gpp=Gpp,constraint=constraint)

	return(RheoModel(m.name, freeparams, fixedparams,
		_buildmoduli_t(G)..., _buildmoduli_t(J)..., 
        _buildmoduli_ω(Gp)..., _buildmoduli_ω(Gpp)...,
		_buildconstraint(constraint),
		m._Gramp, m.info, expressions) )
end


function RheoModel(m::RheoModelClass; kwargs...)
    return(RheoModel(m,symbol_to_unicode(values(kwargs))))
end


"""
    getparams(m::RheoModel; unicode=true)

`getparams` return the list of model parameters with their values as a NamedTuple.
If `unicode` is set to `false`, the unicode symbols are converted their text equivalent.

# Example
```@example
julia> model = RheoModel(Maxwell, k=1, η=2.)
[...]

julia> getparams(m)
(k = 1.0, η = 2.0)

julia> getparams(m,unicode=false)
(k = 1.0, eta = 2.0)
``` 
"""
function getparams(m::RheoModel; unicode=true)
    if unicode
        return(m.fixedparams)
    else
        return(unicode_to_text(m.fixedparams))
    end
end

#=
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

    expressions=NamedTuple{(:G, :J, :Gp, :Gpp, :constraint, :Ga_safe, :Ja_safe)}( ( G, J, Gp, Gpp, constraint, m.expressions.Ga_safe, m.expressions.Ja_safe ) )

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







struct RheoModel

    _G::FunctionWrapper{RheoFloat,Tuple{RheoFloat}}
    _Ga::FunctionWrapper{Vector{RheoFloat},Tuple{Vector{RheoFloat}}}
    _J::FunctionWrapper{RheoFloat,Tuple{RheoFloat}}
    _Ja::FunctionWrapper{Vector{RheoFloat},Tuple{Vector{RheoFloat}}}
    _Gp::FunctionWrapper{RheoFloat,Tuple{RheoFloat}}
    _Gpa::FunctionWrapper{Vector{RheoFloat},Tuple{Vector{RheoFloat}}}
    _Gpp::FunctionWrapper{RheoFloat,Tuple{RheoFloat}}
    _Gppa::FunctionWrapper{Vector{RheoFloat},Tuple{Vector{RheoFloat}}}

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
    p = check_and_reorder_parameters(nt0, m.params, err_string = "model definition")
    nt=NamedTuple{Tuple(m.params)}(p)

    # string printed
    info = string("\nModel: $(m.name)\n\nParameter values: $nt \n",m.info)

    # This attaches expressions to the models with parameters substituted by values.
    # Future development: this could be disabled with a key word if processing time is an issue
    G = expr_replace(m.expressions.G, nt)
    J = expr_replace(m.expressions.J, nt)
    Gp = expr_replace(m.expressions.Gp, nt)
    Gpp = expr_replace(m.expressions.Gpp, nt)

    expressions=NamedTuple{(:G, :J, :Gp, :Gpp, :Ga_safe, :Ja_safe)}( ( G, J, Gp, Gpp, m.expressions.Ga_safe, m.expressions.Ja_safe ) )

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


=#



#
# These are exported functions enabling users to access the different moduli functions for each model
#



"""
    relaxmod(m[, t, params])

provide access to the relaxation modulus (G) of a given model m. It can be used with a broad range of inputs.

When a time value is provided (or an array of time values), the function returns the corresponding value(s) of the relaxation modulus.

    relaxmod(m::RheoModel, time (single value or array))

    relaxmod(m::RheoModelClass, time (single value or array), parameters (array, named tupple or keyword parameters))

Examples:

    relaxmod(Maxwell, 1, k=1., η=1)

    relaxmod(Maxwell, [1,2,3], [1,1])


When no time value is provided, relaxmod returns the relaxation function itself.

    relaxmod(m::RheoModel)

    relaxmod(m::RheoModelClass, parameters (array, named tupple or keyword parameters))

Examples:

    m = RheoModel(Maxwell, k=1., η=1)
    G = relaxmod(m)
    G([0,1,2])
"""
function relaxmod(m::RheoModel, t::Number)
    m._G(t)
end

function relaxmod(m::RheoModel, ta::Vector{T}) where T <: Number
    m._Ga(ta)
end

relaxmod(m::RheoModel) = x -> relaxmod(m,x)




function relaxmod(m::RheoModelClass, t::Number, params::Vector{T}) where T <: Number
    m._G(t, params)
end


function relaxmod(m::RheoModelClass, ta::Vector{T1}, params::Vector{T2}) where {T1 <: Number, T2 <: Number}
    m._Ga(ta, params)
end


function relaxmod(m::RheoModelClass, x, params::NamedTuple)
    # check all parameters are provided and create a well ordered named tuple
    p = check_and_reorder_parameters(params, m.freeparams, err_string = "relaxmod")
    relaxmod(m, x, p)
end


function relaxmod(m::RheoModelClass, x; kwargs...)
    relaxmod(m, x, symbol_to_unicode(values(kwargs)))
end

relaxmod(m::RheoModelClass, params::Vector{T}) where T <: Number =  x -> relaxmod(m, x, params)
relaxmod(m::RheoModelClass, params::NamedTuple) =  x -> relaxmod(m, x, params)
relaxmod(m::RheoModelClass; kwargs...) =  x -> relaxmod(m, x, symbol_to_unicode(values(kwargs)))















"""
    creepcomp(m[, t, params])

provide access to the creep compliance (J) of a given model m. It can be used with a broad range of inputs.

When a time value is provided (or an array of time values), the function returns the corresponding value(s) of the creep compliance.

    creepcomp(m::RheoModel, time (single value or array))

    creepcomp(m::RheoModelClass, time (single value or array), parameters (array, named tupple or keyword parameters))

Examples:

    creepcomp(Maxwell, 1, k=1., η=1)

    creepcomp(Maxwell, [1,2,3], [1,1])


When no time value is provided, creepcomp returns the creep compliance function itself.

    creepcomp(m::RheoModel)

    creepcomp(m::RheoModelClass, parameters (array, named tupple or keyword parameters))

Examples:

    m = RheoModel(Maxwell, k=1., η=1)
    J = creepcomp(m)
    J([0,1,2])
"""
function creepcomp(m::RheoModel, t::Number)
    m._J(t)
end

function creepcomp(m::RheoModel, ta::Vector{T}) where T <: Number
    m._Ja(ta)
end

creepcomp(m::RheoModel) = x -> creepcomp(m,x)




function creepcomp(m::RheoModelClass, t::Number, params::Vector{T}) where T <: Number
    m._J(t, params)
end


function creepcomp(m::RheoModelClass, ta::Vector{T1}, params::Vector{T2}) where {T1 <: Number, T2 <: Number}
    m._Ja(ta, params)
end


function creepcomp(m::RheoModelClass, x, params::NamedTuple)
    # check all parameters are provided and create a well ordered named tuple
    p = check_and_reorder_parameters(params, m.freeparams, err_string = "relaxmod")
    creepcomp(m, x, p)
end


function creepcomp(m::RheoModelClass, x; kwargs...)
  creepcomp(m, x, symbol_to_unicode(values(kwargs)))
end

creepcomp(m::RheoModelClass, params::Vector{T}) where T <: Number =  x -> creepcomp(m, x, params)
creepcomp(m::RheoModelClass, params::NamedTuple) =  x -> creepcomp(m, x, params)
creepcomp(m::RheoModelClass; kwargs...) =  x -> creepcomp(m, x, symbol_to_unicode(values(kwargs)))















"""
    storagemod(m[, ω, params])

provide access to the storage modulus (G') of a given model m.

When a frequency value is provided (or an array of frequency values), the function returns the corresponding value(s) of the storage modulus.

    storagemod(m::RheoModel, frequency (single value or array))

    storagemod(m::RheoModelClass, frequency (single value or array), parameters (array, named tupple or keyword parameters))

Examples:

    storagemod(Maxwell, 1, k=1., η=1)

    storagemod(Maxwell, [1,2,3], [1,1])


When no time value is provided, storagemod returns the storage modulus function itself.

    storagemod(m::RheoModel)

    storagemod(m::RheoModelClass, parameters (array, named tupple or keyword parameters))

Examples:

    m = RheoModel(Maxwell, k=1., η=1)
    Gp = storagemod(m)
    Gp([0,1,2])
"""
function storagemod(m::RheoModel, ω::Number)
    m._Gp(ω)
end

function storagemod(m::RheoModel, ωa::Vector{T}) where T <: Number
    m._Gpa(ωa)
end

storagemod(m::RheoModel) = x -> storagemod(m,x)




function storagemod(m::RheoModelClass, ω::Number, params::Vector{T}) where T <: Number
    m._Gp(ω, params)
end


function storagemod(m::RheoModelClass, ωa::Vector{T1}, params::Vector{T2}) where {T1 <: Number, T2 <: Number}
      m._Gpa(ωa, params)
end


function storagemod(m::RheoModelClass, x, params::NamedTuple)
    # check all parameters are provided and create a well ordered named tuple
    p = check_and_reorder_parameters(params, m.freeparams, err_string = "relaxmod")
    storagemod(m, x, p)
end


function storagemod(m::RheoModelClass, x; kwargs...)
    storagemod(m, x, symbol_to_unicode(values(kwargs)))
end

storagemod(m::RheoModelClass, params::Vector{T}) where T <: Number =  x -> storagemod(m, x, params)
storagemod(m::RheoModelClass, params::NamedTuple) =  x -> storagemod(m, x, params)
storagemod(m::RheoModelClass; kwargs...) =  x -> storagemod(m, x, symbol_to_unicode(values(kwargs)))





"""
    lossmod(m[, ω, params])

provide access to the loss modulus (G'') of a given model m.

When a frequency value is provided (or an array of frequency values), the function returns the corresponding value(s) of the loss modulus.

    lossmod(m::RheoModel, frequency (single value or array))

    lossmod(m::RheoModelClass, frequency (single value or array), parameters (array, named tupple or keyword parameters))

Examples:

    lossmod(Maxwell, 1, k=1., η=1)

    lossmod(Maxwell, [1,2,3], [1,1])


When no time value is provided, lossmod returns the loss modulus function itself.

    lossmod(m::RheoModel)

    lossmod(m::RheoModelClass, parameters (array, named tupple or keyword parameters))

Examples:

    m = RheoModel(Maxwell, k=1., η=1)
    Gpp = lossmod(m)
    Gpp([0,1,2])
"""
function lossmod(m::RheoModel, ω::Number)
    m._Gpp(ω)
end

function lossmod(m::RheoModel, ωa::Vector{T}) where T <: Number
    m._Gppa(ωa)
end

lossmod(m::RheoModel) = x -> lossmod(m,x)




function lossmod(m::RheoModelClass, ω::Number, params::Vector{T}) where T <: Number
    m._Gpp(ω, params)
end


function lossmod(m::RheoModelClass, ωa::Vector{T1}, params::Vector{T2}) where {T1 <: Number, T2 <: Number}
    m._Gppa(ωa, params)
end


function lossmod(m::RheoModelClass, x, params::NamedTuple)
    # check all parameters are provided and create a well ordered named tuple
    p = check_and_reorder_parameters(params, m.freeparams, err_string = "relaxmod")
    lossmod(m, x, p)
end


function lossmod(m::RheoModelClass, x; kwargs...)
    lossmod(m, x, symbol_to_unicode(values(kwargs)))
end

lossmod(m::RheoModelClass, params::Vector{T}) where T <: Number =  x -> lossmod(m, x, params)
lossmod(m::RheoModelClass, params::NamedTuple) =  x -> lossmod(m, x, params)
lossmod(m::RheoModelClass; kwargs...) =  x -> lossmod(m, x, symbol_to_unicode(values(kwargs)))






"""
    dynamicmod(m[, ω, params])

provide access to the complex dynamic modulus (G' + i G'') of a given model m.

* When a frequency value is provided (or an array of frequency values), the function returns the corresponding value(s) of the complex dynamic modulus.

    dynamicmod(m::RheoModel, frequency (single value or array))

    dynamicmod(m::RheoModelClass, frequency (single value or array), parameters (array, named tupple or keyword parameters))

Examples:

    dynamicmod(Maxwell, 1, k=1., η=1)

    dynamicmod(Maxwell, [1,2,3], [1,1])


* When no time value is provided, lossmod returns the complex dynamic modulus function itself.

    dynamicmod(m::RheoModel)

    dynamicmod(m::RheoModelClass, parameters (array, named tupple or keyword parameters))

Examples:

    m = RheoModel(Maxwell, k=1., η=1)
    Gd = dynamicmod(m)
    Gd([0,1,2])

Note: use abs() and angle() to get the magnitude and phase of the complex modulus.
"""
function dynamicmod(m::RheoModel, ω::Number)
    m._Gp(ω) + m._Gpp(ω) * im
end

function dynamicmod(m::RheoModel, ωa::Vector{T}) where T <: Number
    m._Gpa(ωa) + m._Gppa(ωa) * im
end

dynamicmod(m::RheoModel) = x -> dynamicmod(m,x)




function dynamicmod(m::RheoModelClass, ω::Number, params::Vector{T}) where T <: Number
    m._Gp(ω, params) + m._Gpp(ω, params) * im
end


function dynamicmod(m::RheoModelClass, ωa::Vector{T1}, params::Vector{T2}) where {T1 <: Number, T2 <: Number}
    m._Gpa(ωa, params) + m._Gppa(ωa, params) * im
end


function dynamicmod(m::RheoModelClass, x, params::NamedTuple)
    # check all parameters are provided and create a well ordered named tuple
    p = check_and_reorder_parameters(params, m.freeparams, err_string = "relaxmod")
    dynamicmod(m, x, p)
end


function dynamicmod(m::RheoModelClass, x; kwargs...)
    dynamicmod(m, x, symbol_to_unicode(values(kwargs)))
end

dynamicmod(m::RheoModelClass, params::Vector{T}) where T <: Number =  x -> dynamicmod(m, x, params)
dynamicmod(m::RheoModelClass, params::NamedTuple) =  x -> dynamicmod(m, x, params)
dynamicmod(m::RheoModelClass; kwargs...) =  x -> dynamicmod(m, x, symbol_to_unicode(values(kwargs)))

