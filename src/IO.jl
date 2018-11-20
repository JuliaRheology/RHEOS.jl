#!/usr/bin/env julia

"""
    fileload(colnames::Array{String,1}, filedir::String)

Load data from a CSV file (three columns, comma seperated). Columns must
be identified by providing an array of strings which tell the function
which data (stress, strain or time) is contained in each column. This is
used to construct a RheologyData struct, which provides the basis for
subsequent high level operations within RHEOS.

# Example

```jldoctest
# directory path to the file
fileDir = "../data/rheologyData1.csv"

# load the data into RheologyData struct
dataforprocessing = fileload(["time","stress","strain"], fileDir)
```
"""
function fileload(colnames::Array{String,1}, filedir::String; delimiter=',')

    # get length 
    N = length(colnames)

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
    savedata(self::RheologyData; filedir::String = "", ext = "_RheologyData.jld")

Save RheologyData object using JLD format. Save file directory
must be specified. If data was loaded from disk using fileload
or the full RheologyData constructor then filedir argument can 
be set to empty string "" which will try to use the original file 
dir concatenated with the optional ext string argument  - e.g. 
"/originalpathto/file.csv_RheologyData.jld".
"""
function savedata(self::Union{RheologyData, RheologyDynamic}, filedir::String; ext = ".jld2")

    fulldir = string(filedir, ext)

    @save fulldir self 

end

"""
    loaddata(filedir::String)

Convenience function loads RheologyData.
"""
function loaddata(filedir::String)

    @load filedir self

    return self

end

"""
    savemodel(self::RheologyModel; filedir::String = "", ext = "")
"""
function savemodel(self::RheologyModel, filedir::String; ext = ".jld2")

    fulldir = string(filedir, ext)

    @save fulldir self 

end

function loadmodel(filedir::String)

    @load filedir self

    return self

end

"""
    exportdata(self::RheologyData; filedir::String = "", ext = ".csv")

Export RheologyData to csv format. Exports three columns in order: stress, strain, time. Or: Gp, Gpp, frequency for oscillatory data.
Useful for plotting/analysis in other software.
"""
function exportdata(self::Union{RheologyData, RheologyDynamic}, filedir::String; ext = ".csv", delimiter=',')

    fulldir = string(filedir, ext)

    if typeof(self)==RheologyData
        fulldata_array = hcat(self.σ, self.ϵ, self.t)
        open(fulldir, "w") do io
            writedlm(fulldir, fulldata_array, delimiter)
        end
            

    elseif typeof(self)==RheologyDynamic
        fulldata_array = hcat(self.Gp, self.Gpp, self.ω)
        open(fulldir, "w") do io
            writedlm(fulldir, fulldata_array, delimiter)
        end
    end
end
