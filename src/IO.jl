#!/usr/bin/env julia

"""
    importdata(colnames::Array{String,1}, filedir::String; delimiter=',')

Load data from a CSV file (two/three columns, comma seperated by default but
delimiter can be specified in the `delimiter` keyword argument). Columns must
be identified by providing an array of strings which tell the function
which data is contained in each column. 

Can be used to construct either a RheologyData instance or a RheologyDynamic
instance. Function detects whether "time" or "frequency" has been included
and proceeds accordingly. For oscillatory data, all three columns (Gp, Gpp, Frequency)
must be provided. For regular viscoelastic data two or three columns can be provided
but time must always be provided. I.e. either stress/strain/time, stress/time or 
strain/time.
"""
function importdata(colnames::Array{String,1}, filedir::String; delimiter=',')

    # get length 
    N = length(colnames)

    @assert ("frequency" in colnames)||("time" in colnames) "Data must contain \"time\" or \"frequency\" "

    if "time" in colnames
        # check colnames length is correct
        @assert N==3 || N==2 "Two or three column names required, one of each: \"stress\", \"strain\" and \"time\"."

        # read data from file
        data = readdlm(filedir, delimiter)

        # generate RheologyData struct and output
        if N==3
            return RheologyData(colnames, data[:, 1], data[:, 2], data[:, 3]; filedir = filedir)
        elseif N==2
            return RheologyData(colnames, data[:, 1], data[:, 2]; filedir = filedir)
        end
    
    elseif "frequency" in colnames
        # check colnames length is correct
        @assert N==3 "Three column names required, one of each: \"Gp\", \"Gpp\" and \"frequency\"."

        # read data from file
        data = readdlm(filedir, delimiter)

        # generate RheologyDynamic struct and output
        return RheologyDynamic(colnames, data[:, 1], data[:, 2], data[:, 3]; filedir=filedir)

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
