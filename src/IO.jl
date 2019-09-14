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

    if any(i -> !i, notnanlist)
        # raise warning
        @warn "Please note that NaN data rows are not included in resultant data struct."
    end

    # extract only those rows whose sum did not evaluate to NaN
    arr[notnanlist, :]
end

"""
    importcsv(filepath::String; t_col::IntOrNone = nothing, σ_col::IntOrNone = nothing, ϵ_col::IntOrNone = nothing, ω_col::IntOrNone = nothing, Gp_col::IntOrNone = nothing, Gpp_col::IntOrNone = nothing, delimiter=',')

Load data from a CSV file (two/three columns, comma seperated by default but
delimiter can be specified in the `delimiter` keyword argument). Arguments must be
identified by providing the number of the column in which they are contained.

Can be used to construct either a RheoTimeData instance or a RheoFreqData
instance. Function detects whether time or frequency has been included
and proceeds accordingly. For oscillatory data, all three columns (Gp, Gpp, Frequency)
must be provided. For regular viscoelastic data only time, or time-stress, or time-strain or
time-stress-strain data can be provided.
"""
function importcsv(filepath::String; t_col::IntOrNone = nothing, σ_col::IntOrNone = nothing, ϵ_col::IntOrNone = nothing, ω_col::IntOrNone = nothing, Gp_col::IntOrNone = nothing, Gpp_col::IntOrNone = nothing, delimiter = ',', comment = "Imported from csv file", savelog = true)

    @assert (!isnothing(t_col) && isnothing(ω_col)) || (isnothing(t_col) && !isnothing(ω_col)) "Data must contain either \"time\" or \"frequency\" "
    @assert endswith(lowercase(filepath), ".csv") "filepath must point to a .csv file."

    if !isnothing(t_col)

        @assert isnothing(Gp_col) && isnothing(Gpp_col) "Loss and storage modulus not allowed for time data"
        # read data from file
        dataraw = readdlm(filepath, delimiter)

        # remove and rows with NaN values in any of the columns
        data = nanremove(dataraw)

        #info=(comment=comment, folder=pwd(), stats=(t_min=data[1,t_col],t_max=data[end,t_col], n_sample=size(data[:,t_col])))
        #log = RheoLogItem( (type=:source, funct=:importcsv, params=(filepath=filepath,), keywords=(t_col=t_col, σ_col=σ_col, ϵ_col=ϵ_col)), info )

        log = if savelog
                info = (comment=comment, folder=pwd(), stats=(t_min=data[1,t_col],t_max=data[end,t_col], n_sample=size(data[:,t_col])))
                RheoLogItem( (type=:source, funct=:importcsv, params=(filepath=filepath,), keywords=(t_col=t_col, σ_col=σ_col, ϵ_col=ϵ_col)), info )
              else
                nothing
              end

        # generate RheologyData struct and output
        if !isnothing(ϵ_col) && !isnothing(σ_col)
            return RheoTimeData(t = data[:, t_col], ϵ = data[:, ϵ_col], σ = data[:, σ_col], log = log)
        elseif !isnothing(ϵ_col) && isnothing(σ_col)
            return RheoTimeData(t = data[:, t_col], ϵ = data[:, ϵ_col], log = log)
        elseif !isnothing(σ_col) && isnothing(ϵ_col)
            return RheoTimeData(t = data[:, t_col], σ = data[:, σ_col], log = log)
        end

    elseif !isnothing(ω_col)

        # check colnames length is correct
        @assert isnothing(σ_col) || isnothing(ϵ_col) "Stress and strain not allowed for frequency data"
        @assert !isnothing(Gp_col) && !isnothing(Gpp_col) "\"Gp\" and \"Gpp\" are required."

        # read data from file
        dataraw = readdlm(filepath, delimiter)
        # remove and rows with NaN values in any of the columns
        data = nanremove(dataraw)

        log = if savelog
                info = (comment=comment, folder=pwd(), stats=(ω_min=data[1,ω_col],ω_max=data[end,ω_col], n_sample=size(data[:,ω_col])))
                RheoLogItem( (type=:source, funct=:importcsv, params=(filepath=filepath,), keywords=(ω_col=ω_col, Gp_col=Gp_col, Gpp_col=Gpp_col)), info )
              else
                nothing
              end

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

    @assert endswith(lowercase(filedir), ".csv") "filedir must be a .csv file."

    # if ordering of columns not provided, get default ordering depending on data contained in struct
    if isnothing(colorder)
        if typeof(self)==RheoTimeData
            datacontained = RheoTimeDataType(self)
            if datacontained == stress_only
                colorder = (t=1, σ=2)
            elseif datacontained == strain_only
                colorder = (t=1, ϵ=2)
            elseif datacontained == strain_and_stress
                colorder = (t=1, σ=2, ϵ=3)
            end

        elseif typeof(self)==RheoFreqData
            # invalid_freq_data=-1 freq_only=0 with_modulus=1
            datacontained = RheoFreqDataType(self)
            colorder = (ω = 1, Gp = 2, Gpp = 3)
        end
    end

    cols = collect(keys(colorder))
    sortedcols = sort(cols, by=i->getfield(colorder, i))
    dataraw = getfield.((self,), sortedcols)
    dataout = hcat(dataraw...)
    open(filedir, "w") do io
        writedlm(filedir, dataout, delimiter)
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
