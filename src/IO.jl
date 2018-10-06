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
function fileload(colnames::Array{String,1}, filedir::String)

    # get length 
    N = length(colnames)

    # check colnames length is correct
    @assert N==3 || N==2 "Two or three column names required, one of each: 'stress', 'strain' and 'time'."

    # init types helper
    types = [Float64 for i = 1:N]

    # read data from file
    (datacsv, head_out) = uCSV.read(filedir; delim=',', types=types)

    # generate RheologyData struct and output
    if N==3
        return RheologyData(colnames, datacsv[1], datacsv[2], datacsv[3]; filedir = filedir)
    elseif N==2
        return RheologyData(colnames, datacsv[1], datacsv[2]; filedir = filedir)
    end
end

# """
#     savedata(self::RheologyData; filedir::String = "", ext = "_RheologyData.jld")

# Save RheologyData object using JLD format. Save file directory
# must be specified. If data was loaded from disk using fileload
# or the full RheologyData constructor then filedir argument can 
# be set to empty string "" which will try to use the original file 
# dir concatenated with the optional ext string argument  - e.g. 
# "/originalpathto/file.csv_RheologyData.jld".
# """
# function savedata(self::RheologyData; filedir::String = "", ext = "_RheologyData.jld")

#     if filedir == ""
#         if self.log[2] == "none"
#             error("No file directory provided in savedata method or embedded in RheologyData object.")

#         else
#             _filedir = self.log[2]

#         end

#     else
#         _filedir = filedir

#     end

#     fulldir = string(_filedir, ext)

#     jldopen(fulldir, "w") do file
#         # register RHEOS module (with RheologyData type) to JLD
#         addrequire(file, RHEOS)
#         # write self to file
#         write(file, "self", self)
#     end

# end

# """
#     loaddata(filedir::String)

# Convenience function loads RheologyData.
# """
# function loaddata(filedir::String)
#     # load in result
#     loaded_data = load(filedir)
#     # return
#     loaded_data["self"]

# end

# """
#     savemodel(self::RheologyModel; filedir::String = "", ext = "")
# """
# function savemodel(self::RheologyModel; filedir::String = "", ext = "")

#     #########################################################################
#     #=
#     TEMPORARY STRUCT AS A WORKAROUND FOR THIS JLD ISSUE, FUNCTIONS or STRUCTS CONTAINING
#     FUNCTIONS CANNOT BE SAVED. SEE https://github.com/JuliaIO/JLD.jl/issues/57 FOR MORE
#     INFORMATION. 
#     =#
#     self = RheologyModelTemp(string(self.modulus), self.parameters, self.log)
#     #########################################################################

#     if filedir == ""
#         if self.log[2] == "none"
#             error("No file directory provided in savedata method or embedded in RheologyModel object.")

#         else
#             _filedir = self.log[2]

#         end

#     else
#         _filedir = filedir

#     end

#     if ext == ""
#         _ext = string("_", string(self.modulus), ".jld")

#     else
#         _ext = ext

#     end

#     fulldir = string(_filedir, _ext)

#     jldopen(fulldir, "w") do file
#         # register RHEOS module (with RheologyModel type) to JLD
#         addrequire(file, RHEOS)
#         # write self to file
#         write(file, "self", self)
#     end

# end

# """
#     loadmodel(filedir::String)

# Convenience function loads RheologyModel from disk.
# """
# function loadmodel(filedir::String)

#     # load in result
#     loaded_data = load(filedir)

#     #########################################################################
#     #=
#     TEMPORARY STRUCT AS A WORKAROUND FOR THIS JLD ISSUE, FUNCTIONS or STRUCTS CONTAINING
#     FUNCTIONS CANNOT BE SAVED. SEE https://github.com/JuliaIO/JLD.jl/issues/57 FOR MORE
#     INFORMATION. 
#     =#
#     self_text = loaded_data["self"]

#     self_proper = RheologyModel(getfield(Main, Symbol(self_text.modulus[7:end])), self_text.parameters, self_text.log)

#     # can substitute for below line when issue is resolved
#     # loaded_data["self"]
#     ######################################################################### 

# end

# """
#     exportdata(self::RheologyData; filedir::String = "", ext = "_mod.csv")

# Export RheologyData to csv format. Exports three columns in order: stress, strain, time.
# Useful for plotting/analysis in other software.
# """
# function exportdata(self::RheologyData; filedir::String = "", ext = "_mod.csv")

#     if filedir == ""
#         if self.log[2] == "none"
#             error("No file directory provided in savedata method or embedded in RheologyData object.")

#         else
#             _filedir = self.log[2]

#         end

#     else
#         _filedir = filedir

#     end

#     fulldir = string(_filedir, ext)

#     fulldata_array = hcat(self.σ, self.ϵ, self.t)

#     fulldata_frame = convert(DataFrame, fulldata_array)

#     names!(fulldata_frame, [:stress, :strain, :time])

#     # array to write
#     # fulldata = [self.σ, self.ϵ, self.t]

#     uCSV.write(fulldir, fulldata_frame)

# end
