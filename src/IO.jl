#!/usr/bin/env julia

"""
    nanremove(arr::Array{T,2}) where T<:Real

For a real 2D array with columns of e.g. time, stress and strain, or frequency, G' and G'',
find any row which contains NaN values and remove that row completely.
"""
function nanremove(arr::Array{T,2}) where T<:Real
    # get 1D array which is NaN for any corresponding row with NaNs in
    rowsums = sum(arr, dims=2)

    # get list of NaNs / not NaNs (true for not NaN, false for NaN)
    notnanlist = vec(.!isnan.(rowsums))

    if false in notnanlist
        # raise warning
        @warn "Please note that NaN data rows are not included in resultant data struct."
    end

    # extract only those rows whose sum did not evaluate to NaN
    arr[notnanlist, :]
end



function look_for_csv_col(var, cols, header_cells, default_options)
	col = nothing 

	if var in keys(cols)
		if typeof(cols[var])<:String
			col = findfirst(e -> e == cols[var], header_cells)
		elseif typeof(cols[var])<:Integer
			col = cols[var]
		end

	else
		col = findfirst(e -> lowercase(e) in default_options, header_cells)
	end
	return(col)
end



"""
    importcsv(filepath::String; delimiter=',', header=false, comment = "Imported from csv file", savelog = true, kwargs...)

Load data from a CSV file (two/three columns, comma separated by default but
delimiter can be specified in the `delimiter` keyword argument). 

The function can be used to construct either a RheoTimeData instance or a RheoFreqData
instance, depending on the parameters provided or column names;the function detects whether time or frequency has been included
and proceeds accordingly. For oscillatory data, all three columns (Gp, Gpp, Frequency)
must be provided. For regular viscoelastic data only time, or time-stress, or time-strain or
time-stress-strain data can be provided.

Column can be specified in the csv file by providing their values or header name (case insensitive) as keyword parameters. 

`importcsv("filename.csv", time=1, strain=2, stress=3)` loads a `RheoTimeData` from a csv file with time in the first column, strain in the second, and stress in the third.

`importcsv("filename.csv", omega="Frequency", Gp="Storage", Gpp="Loss")` would detect the column numbers by reading the headers on the csv file.

If no information is provided, it would try to detect columns from header names, expecting standard names. Some of the recognised keywords are `time`, `stress`, `strain`, `frequency`, `storage modulus`, `loss modulus`.

If the csv file contains headers, but it is prefered to indicate columns by their numbers rather than header strings, the keyword `header` must be set to `true`.
"""
function importcsv(filepath::String; delimiter = ',', header = false, comment = "Imported from csv file", savelog = true, kwargs...)
# Convert parameter names to standard values.
cols=symbol_to_unicode(values(kwargs))

	use_header = header
	if length(cols)==0
		use_header = true
	elseif any([typeof(e)<:String for e in cols])
		use_header = true
	end

	# read data from file
    if use_header
		dataraw,h = readdlm(filepath, delimiter, header=use_header)
        header_cells = [strip(string(e)) for e in [h...]]
	else
		dataraw = readdlm(filepath, delimiter)
        header_cells=[]
	end


    # remove and rows with NaN values in any of the columns
    data = nanremove(dataraw)

	# is it time data?
	t_col = look_for_csv_col(:t, cols, header_cells, ("t", "time"))

	if !isnothing(t_col)
		σ_col = look_for_csv_col(:σ, cols, header_cells, ("stress", "sigma", "σ"))
        ϵ_col = look_for_csv_col(:ϵ, cols, header_cells, ("strain", "epsilon", "ϵ"))
              
		
    #  @assert (!isnothing(t_col) && isnothing(ω_col)) || (isnothing(t_col) && !isnothing(ω_col)) "Data must contain either \"time\" or \"frequency\" "
    #  @assert endswith(lowercase(filepath), ".csv") "filepath must point to a .csv file."

        log = loginit(savelog, :importcsv, params=(filepath=filepath,), keywords=(t=t_col, σ=σ_col, ϵ=ϵ_col, header=header),
                      info = (comment=comment, folder=pwd(), stats=(t_min=data[1,t_col],t_max=data[end,t_col], n_sample=size(data[:,t_col]))))

        # generate RheologyData struct and output
        if !isnothing(ϵ_col) && !isnothing(σ_col)
            return RheoTimeData(t = data[:, t_col], ϵ = data[:, ϵ_col], σ = data[:, σ_col], log = log)
        elseif !isnothing(ϵ_col) && isnothing(σ_col)
            return RheoTimeData(t = data[:, t_col], ϵ = data[:, ϵ_col], log = log)
        elseif !isnothing(σ_col) && isnothing(ϵ_col)
            return RheoTimeData(t = data[:, t_col], σ = data[:, σ_col], log = log)
        end
    end

	ω_col = look_for_csv_col(:ω, cols, header_cells, ("ω", "omega", "frequency", "freq"))

    if !isnothing(ω_col)

        # check colnames length is correct
        #@assert isnothing(σ_col) || isnothing(ϵ_col) "Stress and strain not allowed for frequency data"
        #@assert !isnothing(Gp_col) && !isnothing(Gpp_col) "\"Gp\" and \"Gpp\" are required."
        Gp_col = look_for_csv_col(:Gp, cols, header_cells, ("Gp", "G'", "storage modulus"))
        Gpp_col = look_for_csv_col(:Gpp, cols, header_cells, ("Gpp", "G''", "loss modulus"))

        # Add error message for missing loss moduli / error loading
 
        log = loginit(savelog, :importcsv, params=(filepath=filepath,), keywords=(ω=ω_col, Gp=Gp_col, Gpp=Gpp_col, header=header),
                      info = (comment=comment, folder=pwd(), stats=(ω_min=data[1,ω_col],ω_max=data[end,ω_col], n_sample=size(data[:,ω_col]))))

      
        return RheoFreqData(ω = data[:,ω_col], Gp = data[:,Gp_col], Gpp = data[:,Gpp_col], log = log)

    end
end







"""
    exportcsv(self::Union{RheoTimeData, RheoFreqData}, filedir::String; delimiter=',', colorder=nothing)

Export `RheoTimeData` or `RheoFreqData` type to csv format. May be useful for plotting/analysis in other software.
By default, full time data will be exported with columns ordered as (t, σ, ϵ). Partial time data will be ordered
as either (t, σ) or (t, ϵ). Full frequency data will be ordered as (ω, Gp, Gpp). The order of columns can be customised
by passing a NamedTuple to the `colorder` arguments. For example (σ = 1, t = 3, ϵ = 2) would export the columns in the
order (σ, ϵ, t). As with `importcsv`, the delimiter can be set by keyword argument.
"""
function exportcsv(self::Union{RheoTimeData, RheoFreqData}, filedir::String; delimiter=',', colorder=nothing)

    #@assert endswith(lowercase(filedir), ".csv") "filedir must be a .csv file."

    # if ordering of columns not provided, get default ordering depending on data contained in struct
    if isnothing(colorder)
        if typeof(self)==RheoTimeData
            datacontained = rheotimedatatype(self)
            if datacontained == stress_only
                colorder = (t=1, σ=2)
            elseif datacontained == strain_only
                colorder = (t=1, ϵ=2)
            elseif datacontained == strain_and_stress
                colorder = (t=1, σ=2, ϵ=3)
            end

        elseif typeof(self)==RheoFreqData
            # invalid_freq_data=-1 freq_only=0 with_modulus=1
            datacontained = rheofreqdatatype(self)
            colorder = (ω = 1, Gp = 2, Gpp = 3)
        end
    end

    cols = collect(keys(colorder))
    sortedcols = sort(cols, by=i->getfield(colorder, i))
    dataraw = getfield.((self,), sortedcols)
    dataout = hcat(dataraw...)
    open(filedir, "w") do io
        writedlm(io, dataout, delimiter)
    end
end

# """
#     savedata(self::Union{RheoTimeData, RheoFreqData}, filedir::String)

# Save RheoTimeData or RheoFreqData object using JLD2 format for re-use in
# a later Julia session.
# """
# function savedata(self::Union{RheoTimeData, RheoFreqData}, filedir::String)

#     @assert endswith(lowercase(filedir), ".jld2") "filedir must be a .jld2 file."

#     @save filedir self

# end

# """
#     loaddata(filedir::String)

# Loads RheoTimeData or RheoFreqData object from a jld2 file.
# """
# function loaddata(filedir::String)

#     @assert endswith(lowercase(filedir), ".jld2") "filedir must be a .jld2 file."

#     @load filedir self

#     return self

# end

# """
#     savemodel(self::RheologyModel, filedir::String; ext = ".jld2")

# Save RheologyModel object using JLD2 format reuse in a later Julia session.
# """
# function savemodel(self::RheologyModel, filedir::String; ext = ".jld2")

#     fulldir = string(filedir, ext)

#     @save fulldir self

# end

# """
#     loaddata(filedir::String)

# Loads RheologyModel from a JLD2 file.
# """
# function loadmodel(filedir::String)

#     @load filedir self

#     return self

# end
