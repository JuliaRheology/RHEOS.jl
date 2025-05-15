
"""
    Interface(uϵ::Symbol,uσ::Symbol,from_ϵσ::Function,to_ϵσ::Function)

The `Interface` struct contain information to map stress and strain to other related variables such as force and displacement in experimental systems such as atomic force microscopes or tweezers.

"""
struct Interface
		  uϵ::Symbol
		  uσ::Symbol
		  from_ϵσ::Function
		  to_ϵσ::Function
end




"""
    AFM(R::Real)

Create and return an `Interface` struct using the Hertz contact model with indenter radius `R`.
The AFM interface assumes an incompressible material (Poisson's ratio = 0.5)
"""
function AFM(R::Real)
  Interface(:d, :f, (ϵ,σ)->(d = (R.*(3. /8.).^(2. /3)) .* ϵ.^(2 ./3), f = R^2 .* σ ), (d,f)->(ϵ = 8 ./(3 .*R.^1.5) .* d.^1.5, σ = f./R^2) )
end


# function AFM_correction(R::Real, h::Real)
#   Interface(:d, :f, (ϵ,σ,d)->(d = (R.*(3. /8.).^(2. /3)) .* ϵ.^(2 ./3), f = R^2 .* σ ), (d,f)->(ϵ = 8 ./(3 .*R.^1.5) .* d.^1.5, σ = f./R^2) )
# end

"""
    Tweezers(R::Real, slip::Real)

Create and return an `Interface` struct for a spherical probe of radius `R` embeded in a material, with probe/material viscosity ratio `slip`.
"""
function Tweezers(R::Real, slip::Real=Inf)
  # slip can go between 0 and Inf 
  # no slip is at Inf, frictionless slip is at 0
  R >= 0 || error("sphere must have a positive radius")
  slip >= 0 || error("slip less than 0 is not physical")

  if slip == Inf
    hadamard = 3. / 2.
  else
    hadamard = (((3 * slip) + 2) / (2 * (slip + 1)))
  end

  stress_conversion = hadamard * 4 * π * R^2

  Interface(:d, :f, (ϵ,σ)->(d = R .* ϵ, f = stress_conversion .* σ ), (d,f)->(ϵ = d ./ R, σ = f ./ stress_conversion) )
end


function Base.getindex(d::RheoTimeData, i::Interface)
	return( NamedTuple{(i.uϵ, i.uσ),Tuple{Vector{RheoFloat},Vector{RheoFloat}}}( i.from_ϵσ(d.ϵ,d.σ) )   )
end



function RheoTimeData(interface::Interface ; t::Vector{T} = RheoFloat[], comment="Created from generic constructor", savelog = true, log = savelog ? RheoLogItem(comment) : nothing, kwargs...)  where {T<:Real}

   if interface.uϵ in keys(kwargs)
	  uϵ = convert(Vector{RheoFloat}, kwargs[interface.uϵ])
   else
	  uϵ = RheoFloat[]
   end

   if interface.uσ in keys(kwargs)
	  uσ = convert(Vector{RheoFloat}, kwargs[interface.uσ])
   else
	  uσ = RheoFloat[]
   end

  typecheck = check_time_data_consistency(t,uϵ,uσ)
  ϵ,σ = interface.to_ϵσ(uϵ, uσ)

  RheoTimeData(convert(Vector{RheoFloat},σ), convert(Vector{RheoFloat},ϵ), convert(Vector{RheoFloat},t),
  log === nothing ? nothing : [ RheoLogItem(log.action,merge(log.info, (type=typecheck, interface = interface)))]     )

end




"""
    importcsv(filepath::String, interface::Interface; delimiter = ',', header = false, comment = "Imported from csv file", savelog = true, kwargs...)

Import function for raw data to be transformed to stress and strain using an `Interface`.
Column numbers in the csv containging the force/displacement data need to be indicated as keyword arguments using the corresponding symbols in the `interface`..

# Examples
```julia-repl
julia> importcsv("myfile.csv", AFM(2e-6), t=1, d=2, f=3)
```
"""
function importcsv(filepath::String, interface::Interface; delimiter = ',', header = false, comment = "Imported from csv file", savelog = true, kwargs...)
    cols=symbol_to_unicode(values(kwargs))

    # legacy keyword arguments in place until v 1.0
   uϵ_col_sym = Symbol(string(interface.uϵ) * "_col")
   uσ_col_sym = Symbol(string(interface.uσ) * "_col")

   if uϵ_col_sym in keys(cols)
	  uϵ_col = cols[uϵ_col_sym]
   elseif interface.uϵ in keys(cols)
	  uϵ_col = cols[interface.uϵ]
   else
	  uϵ_col = nothing  
   end

   if uσ_col_sym in keys(cols)
	  uσ_col = cols[uσ_col_sym]
    elseif interface.uσ in keys(cols)
      uσ_col = cols[interface.uσ]
    else
	  uσ_col = nothing
   end


   data = importcsv(filepath, t = cols.t, σ = uσ_col, ϵ = uϵ_col, delimiter = delimiter, header = header, comment = comment, savelog = false)

   # need to make sure that missing columns are properly transfered.
   ϵ,σ = interface.to_ϵσ(data.ϵ, data.σ)

   log = if savelog
			 info = (comment=comment, folder=pwd(), stats=(t_min=data.t[1],t_max=data.t[end], n_sample=size(data.t[:])))
			 kwds = NamedTuple{(:t, interface.uϵ, interface.uσ, :header),Tuple{Any,Any,Any,Bool}}( (cols.t, uϵ_col, uσ_col, header) )
			 RheoLogItem( (type=:source, funct=:importcsv, params=(filepath=filepath, interface=interface), keywords=kwds), info )
		 else
			 RheoLogItem(type = nothing, info = "nothing")
		 end

   return RheoTimeData(σ, ϵ, data.t, [log])

end
