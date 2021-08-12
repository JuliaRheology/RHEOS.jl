
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
"""
function AFM(R::Real)
  Interface(:d, :f, (ϵ,σ)->(d = (R.*(3. /4.).^(2. /3)) .* ϵ.^(2. ./3), f = R^2 .* σ ), (d,f)->(ϵ = 4 ./(3 .*R.^1.5) .* d.^1.5, σ = f./R^2) )
end

"""
    Tweezers(R::Real)

Create and return an `Interface` struct for a spherical probe of radius `R` embeded in a material.
"""
function Tweezers(R::Real)
  Interface(:d, :f, (ϵ,σ)->(d = R .* ϵ, f = R^2 .* σ ), (d,f)->(ϵ = d ./ R, σ = f ./ R^2) )
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
   cols=symbol_to_unicode(kwargs.data)

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
			 kwds = NamedTuple{(:t, interface.uϵ, interface.uσ, header),Tuple{IntOrNone,IntOrNone,IntOrNone,Bool}}( (cols.t, uϵ_col, uσ_col, header) )
			 RheoLogItem( (type=:source, funct=:importcsv, params=(filepath=filepath, interface=interface), keywords=kwds), info )
		 else
			 RheoLogItem(type = nothing, info = "nothing")
		 end

   return RheoTimeData(σ, ϵ, data.t, [log])

end
