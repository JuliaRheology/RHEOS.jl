#!/usr/bin/env julia

"""
    nanremove(arr::Array{T,2}) where T<:Real

For a real 2D array with columns of e.g. time, stress and strain, or frequency, G' and G'',
find any row which contains NaN values and remove that row completely.
"""
function nanremove(arr::Array{T,2}) where T<:Real
    # get 1D array which is NaN for any corresponding row with NaNs in
    rowsums = sum(arr, dims=2)

    # extract only those rows whose sum did not evaluate to NaN
    arr[vec(.!isnan.(rowsums)), :]
end

"""
    importcsv(filedir::String; t_col::IntOrNone= nothing, σ_col::IntOrNone= nothing, ϵ_col::IntOrNone=nothing, ω_col::IntOrNone= nothing, Gp_col::IntOrNone = nothing, Gpp_col::IntOrNone = nothing, delimiter=',')

Load data from a CSV file (two/three columns, comma seperated by default but
delimiter can be specified in the `delimiter` keyword argument). Arguments must be
identified by providing the number of the column in which they are contained.

Can be used to construct either a RheoTimeData instance or a RheoFreqData
instance. Function detects whether "time" or "frequency" has been included
and proceeds accordingly. For oscillatory data, all three columns (Gp, Gpp, Frequency)
must be provided. For regular viscoelastic data only time, or time-stress, or time-strain or
time-stress-strain data can be provided.
"""
function importcsv(filedir::String; t_col::IntOrNone= nothing, σ_col::IntOrNone= nothing, ϵ_col::IntOrNone=nothing, ω_col::IntOrNone= nothing, Gp_col::IntOrNone = nothing, Gpp_col::IntOrNone = nothing, delimiter=',')

    @assert (!isnothing(t_col) && isnothing(ω_col)) || (isnothing(t_col) && !isnothing(ω_col)) "Data must contain either \"time\" or \"frequency\" "
    @assert endswith(filedir, ".csv") "filedir must be a .csv file."

    if !isnothing(t_col)

        @assert isnothing(Gp_col) && isnothing(Gpp_col) "Loss and storage modulus not allowed for time data"
        # read data from file
        dataraw = readdlm(filedir, delimiter)

        # remove and rows with NaN values in any of the columns
        data = nanremove(dataraw)

        source = filedir
        # generate RheologyData struct and output
        if !isnothing(ϵ_col) && !isnothing(σ_col)
            return RheoTimeData(t = data[:, t_col], ϵ = data[:, ϵ_col], σ = data[:, σ_col], source = source)
        elseif !isnothing(ϵ_col) && isnothing(σ_col)
            return RheoTimeData(t = data[:, t_col], ϵ = data[:, ϵ_col], source = source)
        elseif !isnothing(σ_col) && isnothing(ϵ_col)
            return RheoTimeData(t = data[:, t_col], σ = data[:, σ_col], source = source)
        end

    elseif !isnothing(ω_col)
        # check colnames length is correct
        @assert isnothing(σ_col) || isnothing(ϵ_col) "Stress and strain not allowed for frequency data"
        @assert !isnothing(Gp_col) && !isnothing(Gpp_col) "\"Gp\" and \"Gpp\" are required."

        # read data from file
        dataraw = readdlm(filedir, delimiter)
        # remove and rows with NaN values in any of the columns
        data = nanremove(dataraw)

        # generate RheologyDynamic struct and output
        return RheoFreqData(ω = data[:,ω_col], Gp = data[:,Gp_col], Gpp = data[:,Gpp_col], source = filedir)

    end
end

"""
    exportcsv(self::Union{RheoTimeData, RheoFreqData}, filedir::String; ext = ".csv", delimiter=',', colorder=nothing)

Export RheoTimeData or RheoFreqData type to csv format. May be useful for plotting/analysis in other software.
By default, full time data will be exported with columns ordered as (σ, ϵ, t). Partial time data will be ordered
as either (σ, t) or (ϵ, t). Full frequency data will be ordered as (Gp, Gpp, ω). The order of columns can be customised
by passing a NamedTuple to the `colorder` arguments. For example (σ = 2, t = 1, ϵ = 3) would export the columns in the
order (t, σ, ϵ). As with `importcsv`, the delimiter can be set by keyword argument.
"""
function exportcsv(self::Union{RheoTimeData, RheoFreqData}, filedir::String; delimiter=',', colorder=nothing)

    @assert endswith(filedir, ".csv") "filedir must be a .csv file."

    # if ordering of columns not provided, get default ordering depending on data contained in struct
    if isnothing(colorder)
        if typeof(self)==RheoTimeData
            datacontained = RheoTimeDataType(self)
            if datacontained == stress_only
                colorder = (σ=1, t=2)
            elseif datacontained == strain_only
                colorder = (ϵ=1, t=2)           
            elseif datacontained == strain_and_stress
                colorder = (σ=1, ϵ=2, t=3)
            end

        elseif typeof(self)==RheoFreqData
            # invalid_freq_data=-1 frec_only=0 with_modulus=1
            datacontained = RheoFreqDataType(self)
            colorder = (Gp = 1, Gpp = 2, ω = 3)
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
#     savedata(self::Union{RheologyData, RheologyDynamic}, filedir::String; ext = ".jld2")

# Save RheologyData or RheologyDynamic object using JLD2 format reuse in
# a later Julia session.
# """
# function savedata(self::Union{RheologyData, RheologyDynamic}, filedir::String; ext = ".jld2")

#     fulldir = string(filedir, ext)

#     @save fulldir self

# end

# """
#     loaddata(filedir::String)

# Loads RheologyData or RheologyDynamic object from a jld2 file.
# """
# function loaddata(filedir::String)

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
