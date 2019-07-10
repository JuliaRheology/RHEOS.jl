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
    importdata(filedir::String; t_col::Integer= -1, σ_col::Integer= -1, ϵ_col::Integer=-1, ω_col::Integer= -1, Gp_col::Integer = -1, Gpp_col::Integer = -1, delimiter=',')

Load data from a CSV file (two/three columns, comma seperated by default but
delimiter can be specified in the `delimiter` keyword argument). Arguments must be
identified by providing the number of the column in which they are contained.

Can be used to construct either a RheoTimeData instance or a RheoFreqData
instance. Function detects whether "time" or "frequency" has been included
and proceeds accordingly. For oscillatory data, all three columns (Gp, Gpp, Frequency)
must be provided. For regular viscoelastic data only time, or time-stress, or time-strain or
time-stress-strain data can be provided.
"""
function importdata(filedir::String; t_col::Union{Integer, Nothing}= nothing, σ_col::Union{Integer, Nothing}= nothing, ϵ_col::Union{Integer, Nothing}=nothing, ω_col::Union{Integer, Nothing}= nothing, Gp_col::Union{Integer, Nothing} = nothing, Gpp_col::Union{Integer, Nothing} = nothing, delimiter=',')

    @assert (!isnothing(t_col) && isnothing(ω_col)) | (isnothing(t_col==-1) && !isnothing(ω_col)) "Data must contain either \"time\" or \"frequency\" "

    if !isnothing(t_col)

        @assert isnothing(Gp_col) & isnothing(Gpp_col) "Loss and storage modulus not allowed for time data"
        # read data from file
        dataraw = readdlm(filedir, delimiter)

        # remove and rows with NaN values in any of the columns
        data = nanremove(dataraw)

        source = filedir
        # generate RheologyData struct and output
        if (ϵ_col!=-1) & (σ_col!=-1)
            return RheoTimeData(t = data[:, t_col], ϵ = data[:, ϵ_col], σ = data[:, σ_col], source = source)
        elseif (ϵ_col!=-1) & (σ_col==-1)
            return RheoTimeData(t = data[:, t_col], ϵ = data[:, ϵ_col], source = source)
        elseif (σ_col!=-1) & (ϵ_col==-1)
            return RheoTimeData(t = data[:, t_col], σ = data[:, σ_col], source = source)
        elseif (t_col!=-1) & (σ_col==-1) &  (ϵ_col==-1)
            return RheoTimeData(t = data[:, t_col], source = source)
        end

    elseif ω_col!=-1
        # check colnames length is correct
        @assert (σ_col==-1) & (ϵ_col ==-1) "Stress and strain not allowed for frequency data"
        @assert Gp_col!=-1 & Gpp_col!=-1 "\"Gp\" and \"Gpp\" are required."

        # read data from file
        dataraw = readdlm(filedir, delimiter)
        # remove and rows with NaN values in any of the columns
        data = nanremove(dataraw)

        # generate RheologyDynamic struct and output
        return RheoFreqData(ω = data[:,ω_col], Gp = data[:,Gp_col], Gpp = data[:,Gpp_col], source = filedir)

    end
end

"""
    exportdata(self::Union{RheologyData, RheologyDynamic}, filedir::String; ext = ".csv", delimiter=',')

Export RheologyData or RheologyDynamic type to csv format. Exports three columns in order:
stress, strain, time for standard viscoelastic data. Or: Gp, Gpp, frequency for oscillatory data.
File extension can be modified using the `ext` keyword argument. As with `importdata`, the delimiter
can also be set by keyword argument.

Useful for plotting/analysis in other software.
"""
function exportdata(self::Union{RheologyData, RheologyDynamic}, filedir::String; ext = ".csv", delimiter=',')

    fulldir = string(filedir, ext)

    subtype = eltype(self)

    if typeof(self)==RheologyData{subtype}
        fulldata_array = hcat(self.σ, self.ϵ, self.t)
        open(fulldir, "w") do io
            writedlm(fulldir, fulldata_array, delimiter)
        end

    elseif typeof(self)==RheologyDynamic{subtype}
        fulldata_array = hcat(self.Gp, self.Gpp, self.ω)
        open(fulldir, "w") do io
            writedlm(fulldir, fulldata_array, delimiter)
        end
    end
end

"""
    savedata(self::Union{RheologyData, RheologyDynamic}, filedir::String; ext = ".jld2")

Save RheologyData or RheologyDynamic object using JLD2 format reuse in
a later Julia session.
"""
function savedata(self::Union{RheologyData, RheologyDynamic}, filedir::String; ext = ".jld2")

    fulldir = string(filedir, ext)

    @save fulldir self

end

"""
    loaddata(filedir::String)

Loads RheologyData or RheologyDynamic object from a jld2 file.
"""
function loaddata(filedir::String)

    @load filedir self

    return self

end

"""
    savemodel(self::RheologyModel, filedir::String; ext = ".jld2")

Save RheologyModel object using JLD2 format reuse in a later Julia session.
"""
function savemodel(self::RheologyModel, filedir::String; ext = ".jld2")

    fulldir = string(filedir, ext)

    @save fulldir self

end

"""
    loaddata(filedir::String)

Loads RheologyModel from a JLD2 file.
"""
function loadmodel(filedir::String)

    @load filedir self

    return self

end
