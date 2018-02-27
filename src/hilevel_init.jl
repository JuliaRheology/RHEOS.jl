#!/usr/bin/env julia

"""
    RheologyData(stress::Array{Float64,1}, strain::Array{Float64,1}, time::Array{Float64,1})

RheologyData mutable struct used for high level interaction with RHEOS
preprocessing and fitting functions. Initialise an instance directly or
indirectly. If data is in three column, comma seperated CSV file then
fileload function can be used, which calls the RheologyData constructor.
If not, load data in according to format and call RheologyData constructor.
"""
mutable struct RheologyData

    # insight parameter, tells functions whether to plot result or not
    # during preprocessing. False by default.
    insight::Bool

    # sampling type. "constant" by default. Use of `var_resample`, `downsample`
    # with more than one section or `fixed_resample` with more than one section
    # overrides to "variable".
    sampling::String

    # test type is either "strlx" (stress relaxation, strain controlled)
    # or "creep" (creep, stress controlled).
    test_type::String

    # original data
    σ::Array{Float64,1}
    ϵ::Array{Float64,1}
    t::Array{Float64,1}

    # stress and strain numerically differentiated WRT to time
    dσ::Array{Float64,1}
    dϵ::Array{Float64,1}

    # models
    fittedmodels::Dict

    # inner constructor forms duplicates of data. All processing is done
    # to duplicates so comparisons can be made between preprocessed and original.
    # If not processing then original data remains preserved ready for use
    # without reptitive coding requirements.
    RheologyData(σ, ϵ, t, test_type) = derivconstruct!(new(false, "constant", test_type, σ, ϵ, t), σ, ϵ, t)

    # inner constructor for resampled data
    RheologyData(_insight, _sampling, _test_type, _σ, _ϵ, _t, _dσ, _dϵ) =
                new(_insight, _sampling, _test_type, _σ, _ϵ, _t, _dσ, _dϵ, Dict())
end

"""
Inner constructor for RheologyData struct, providing gradients of stress
and strain. If stress/strain arrays have NaN values at the beginning (some datasets
have 1 or 2 samples of NaN at beginning) then deletes these and starts at the first
non-NaN sample, also readjusts time start to t = 0 to account for NaNs and and
negative time values at beginning of data recording.
"""
function derivconstruct!(self::RheologyData, σ::Array{Float64,1}, ϵ::Array{Float64,1}, t::Array{Float64,1})

    # define as local so it can be accessed in subsequent scopes
    local newstartingval::T where T<:Integer

    # test for NaNs
    for i in 1:length(self.σ)
        if !isnan(self.σ[i]) && !isnan(self.ϵ[i])
            newstartingval = i
            break
        end
    end

    # adjust starting point accordingly to remove NaNs in σ, ϵ
    initfields = [:σ, :ϵ, :t]
    for n in initfields
        setfield!(self, n, getfield(self, n)[newstartingval:end])
    end

    # derivative of strain WRT time
    self.dσ = deriv(self.σ, self.t)

    # derivative of stress WRT time
    self.dϵ = deriv(self.ϵ, self.t)

    # readjust time to account for NaN movement and/or negative time values
    self.t = self.t - minimum(self.t)

    # initialise empty dictionary for model fit results
    self.fittedmodels = Dict()

    # return class with all fields initialised
    self
end

"""
    fileload(filedir::String, colnames::Array{String,1}, test_type::String)

Load data from a CSV file (three columns, comma seperated). Columns must
be identified by providing an array of strings which tell the function
which data (stress, strain or time) is contained in each column. This is
used to construct a RheologyData struct, which provides the basis for
subsequent high level operations within RHEOS. test_type is either "strlx"
for a stress-relaxation test (strain controlled) or "creep" for a creep
test (stress controlled).

# Example

```jldoctest
# directory path to the file
fileDir = "../data/rheologyData1.csv"

# load the data into RheologyData struct
dataforprocessing = fileload(fileDir, ["time","stress","strain"], "strlx")
```
"""
function fileload(filedir::String, colnames::Array{String,1}, test_type::String)::RheologyData

    # check colnames length is correct
    @assert length(colnames) == 3 "Only three column names required, 'stress', 'strain' and 'time'."

    # read data from file
    (data, head_out) = uCSV.read(filedir; delim=',',
                                 types=[Float64,Float64,Float64])

    # get column numbers in order of RheologyData struct
    cols = zeros(Int8, 3)
    for (i, v) in enumerate(colnames)
        if v == "stress"
            cols[1] = i
        elseif v == "strain"
            cols[2] = i
        elseif v == "time"
            cols[3] = i
        else
            @assert false "Incorrect Column Names"
        end
    end

    # generate RheologyData struct and output
    data = RheologyData(data[cols[1]], data[cols[2]], data[cols[3]], test_type)
end
