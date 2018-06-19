#!/usr/bin/env julia

"""
    RheologyData(σ::Array{Float64,1}, ϵ::Array{Float64,1}, t::Array{Float64,1}, sampling::String, log::Array{String,1})

RheologyData struct used for high level interaction with RHEOS
preprocessing and fitting functions. Initialise an instance directly or
indirectly. If data is in three column, comma separated CSV file then
fileload function can be used, which calls the constructRheologyData function. 
If not, load data in according to format and call constructRheologyData function.
"""
struct RheologyData

    # original data
    σ::Array{Float64,1}
    ϵ::Array{Float64,1}
    t::Array{Float64,1}

    # sampling type. "constant" by default. Use of `var_resample`, `downsample`
    # with more than one section or `fixed_resample` with more than one section
    # overrides to "variable".
    sampling::String

    # operations applied, stores history of which functions (including arguments)
    log::Array{String,1}

end

"""
    constructRheologyData(colnames::Array{String,1}, data1::Array{Float64,1}, data2::Array{Float64,1}[, data3::Array{Float64,1}, filedir::String="none"])::RheologyData

Constructor function for RheologyData struct, if stress/strain arrays have NaN values at the beginning (some datasets
have 1 or 2 samples of NaN at beginning) then deletes these and starts at the first non-NaN sample, also readjusts time start to 
t = 0 to account for NaNs and and negative time values at beginning of data recording.
"""
function constructRheologyData(colnames::Array{String,1}, data1::Array{Float64,1}, data2::Array{Float64,1}, data3::Array{Float64,1}=zeros(length(data2)), filedir::String="none")::RheologyData

    # get data in correct variables
    data = [data1, data2, data3]
    local σ::Array{Float64,1}
    local ϵ::Array{Float64,1}
    local t::Array{Float64,1}

    for (i, v) in enumerate(colnames)
        if v == "stress"
            σ = data[i]
        elseif v == "strain"
            ϵ = data[i]
        elseif v == "time"
            t = data[i]
        else
            @assert false "Incorrect Column Names"
        end
    end

    # define as local so it can be accessed in subsequent scopes
    local newstartingval::T where T<:Integer

    # test for NaNs
    for i in 1:length(σ)
        if !isnan(σ[i]) && !isnan(ϵ[i])
            newstartingval = i
            break
        end
    end

    # adjust starting point accordingly to remove NaNs in σ, ϵ
    σ = σ[newstartingval:end]
    ϵ = ϵ[newstartingval:end]
    t = t[newstartingval:end]

    # Check if time vector is equally spaced
    diff = round.(t[2:end]-t[1:end-1], 4)
    check = any(x->x!=diff[1], diff)
    if check == true
       sampling = "variable"
    else
       sampling = "constant"
    end

    # readjust time to account for NaN movement and/or negative time values
    t = t - minimum(t)

    # initialise empty array to record operations applied during preprocessing
    log = ["loaded file from:", filedir]

    # return class with all fields initialised
    RheologyData(σ, ϵ, t, sampling, log)
end

"""
    fileload(filedir::String, colnames::Array{String,1})

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
dataforprocessing = fileload(fileDir, ["time","stress","strain"])
```
"""
function fileload(filedir::String, colnames::Array{String,1})::RheologyData

    # check colnames length is correct
    @assert length(colnames) == 3 "Three column names required, 'stress', 'strain' and 'time'."

    # read data from file
    (data, head_out) = uCSV.read(filedir; delim=',', types=[Float64,Float64,Float64])

    # generate RheologyData struct and output
    data = constructRheologyData(colnames, data[1], data[2], data[3], filedir)

end

"""
    function stepdata_generate(t_total::Float64, t_on::Float64, t_off::Float64, t_transition::Float64, amplitude::Float64, test_type::String; step_size::Union{String, Float64} = "auto")

Generate partially filled RheologyData struct with a step function approximation
(2 logisitc functions).
"""
function stepdata_generate(t_total::Float64, t_on::Float64, t_off::Float64, t_transition::Float64, amplitude::Float64, test_type::String; step_size::Union{String, Float64} = "auto" )

    if step_size == "auto"
        step_size = t_total/5000.0
    end

    t = collect(0:step_size:t_total)

	k = 10.0/t_transition

	data_controlled = amplitude./(1 + exp.(-k*(t-t_on))) - amplitude./(1 + exp.(-k*(t-t_off)))

    RheologyData(data_controlled, t, test_type, "INTERNAL STEP")

end

"""
    AFMData(stress::Array{Float64,1}, strain::Array{Float64,1}, time::Array{Float64,1}, test_type::String)

AFMData mutable struct used for high level interaction with RHEOS.
Initialise an instance directly or indirectly using the AFMfileload
function.
"""
mutable struct AFMData

    # filedir is location of file on disk, for reference and saving
    filedir::String

    # insight parameter, tells functions whether to plot result or not
    # during preprocessing. False by default.
    insight::Bool

    # sampling type. "constant" by default. Use of `var_resample`, `downsample`
    # with more than one section or `fixed_resample` with more than one section
    # overrides to "variable".
    sampling::String

    # test type is either "strlx" (stress relaxation, displacement controlled)
    # or "creep" (creep, force controlled).
    test_type::String

    # radius of sphere
    R::Float64

    # original data
    f::Array{Float64,1}
    δ::Array{Float64,1}
    t::Array{Float64,1}

    # measured variable to appropriate power
    measured::Array{Float64,1}

    # controlled variable to appropriate power and multiplied by prefactor,
    # numerically differentiated WRT to time
    dcontrolled ::Array{Float64,1}

    # models
    fittedmodels::Dict

    # operations applied, stores history of which functions (including arguments)
    appliedops::Array{String,1}

    # to resample if "all"
    numericdata::Array{Symbol,1}

    # for initial loading of full dataset (measured + prescribed)
    AFMData(f::Array{Float64,1}, δ::Array{Float64,1}, t::Array{Float64,1}, R::Float64, test_type::String, filedir::String) = AFMconstruct!(new(filedir, false, "constant", test_type, R, f, δ, t))

    # inner constructor for incomplete data; data_part should generally
    # be controlled variable (stress for creep, strain for strlx).
    # How to deal with contact point detection?

end

"""
Inner constructor for AFMData struct, fills in fields as appropriate.
"""
function AFMconstruct!(self::AFMData)

    # assume spherical Hertz model
    prefactor = 8.0*sqrt(self.R)/3.0

    # define cm_measured and cm_dcontrolled as appropriate, REF Oyen paper
    if self.test_type == "strlx"

        self.measured = self.f

        self.dcontrolled = prefactor*deriv((self.δ).^(3/2), self.t)

    elseif self.test_type == "creep"

        self.measured = (self.δ).^(3/2)

        self.dcontrolled = deriv(self.f, self.t)/prefactor

    end

    # initialise empty dictionary for model fit results
    self.fittedmodels = Dict()

    # initialise empty array to record operations applied during preprocessing
    self.appliedops = String[]

    # initialise array of sybols for use if resampling/smoothing all numerical data
    self.numericdata = [:f, :δ, :t, :measured, :dcontrolled]

    # return class with appropriate fields initialised
    self
end

"""
Used in AFMdataget pipeline - clean up version number line.
Identify beginning and end of AFM data sections, fancy names
line and version number of JPK data processing software.
"""
function AFMgetlines(filedir::String)

    # initialise array for storing line numbers of section beginnings and endings
    sectionbreaks = Int32[]

    # initialise other local variables
    local fancynamesraw::String
    local vernumraw::String
    local lastlinenum::Int32

    # iterate through lines
    open(filedir) do filJ
        for (i, line) in enumerate(eachline(filJ))

            if startswith(line, "# units")
                append!(sectionbreaks, i+1)

            elseif line == ""
                append!(sectionbreaks, i)

            elseif startswith(line, "# fancyNames:")
                # Get names of columns
                fancynamesraw = line

            elseif startswith(line, "# data-description.modification-software:")
                # get software version number
                vernumraw = line

            elseif eof(filJ)
                lastlinenum = i+1

            end
        end
        # append last line
        append!(sectionbreaks, lastlinenum)
    end

    # get section starts and stops, i.e. shift section break line numbers to account comment lines not read by readtable
    numsections = Int(length(sectionbreaks)/2)

    local blocklen = [0 for i in 1:numsections]
    local blockmod = [0 for i in 1:numsections]
    local blockcount = 0

    for num in (1:numsections)
        # get lengths of each section
        blocklen[num] = sectionbreaks[2:2:end][num] - sectionbreaks[1:2:end][num]
        # accumulate lengths
        blockcount += blocklen[num]
        # store and account for extra values from boundaries
        blockmod[num] = blockcount - num + 1
    end

    # add first line
    sectionbreaksproc = prepend!(blockmod, [1])

    # get sections as unit ranges for ease of use
    local sectionArray  = [0:1 for i in 1:(length(sectionbreaksproc)-1)]

    for i in 1:(length(sectionbreaksproc)-1)
        sectionArray[i] = sectionbreaksproc[i]:(sectionbreaksproc[i+1]-1)
    end

    return sectionArray, fancynamesraw, vernumraw
end

"""
Used in AFMdataget pipeline - clean up version number line.
"""
function versionnumstrip(vernumraw::String)

    return String(split(vernumraw)[3])
end

"""
Used in AFMdataget pipeline - clean up fancy names string so column
titles are in array format.
"""
function fancynamestrip(fancynamesraw::String, vernum::String)

    fancynames_split = split(fancynamesraw,"\"")

    fancynames1 = deleteat!(fancynames_split, 1)

    fancynames = fancynames1[1:2:end]

    local forcecol::Int64
    local timecol::Int64
    local dispcol::Int64

    # use fancynames and vernum to get column numbers of relevant data
    forcecol = find(fancynames .== "Vertical Deflection")[1]

    # if no time column present then return timecol = 0
    if length(find(fancynames .== "Series Time")) == 1
        timecol = find(fancynames .== "Series Time")[1]
    else
        timecol = 0
    end

    _vernum = parse(Int16, vernum[5])

    if _vernum < 6
        dispcol = find(fancynames .== "Tip-Sample Separation")[1]
    elseif _vernum >= 6
        dispcol = find(fancynames .== "Vertical Tip Position")[1]
    end
    return forcecol, dispcol, timecol
end

"""
convenience function for loading in JPK AFM data with relevant metadata
and split in to sections.
"""
function AFMdataget(filedir::String)

    # get section break line numbers, column names and version info
    (sectionarray, fancynamesraw, vernumraw) = AFMgetlines(filedir)

    vernum = versionnumstrip(vernumraw)

    (forcecol, dispcol, timecol) = fancynamestrip(fancynamesraw, vernum)

    # get sectioned data with correct column names
    (data, header) = uCSV.read(filedir; delim=' ', comment='#')

    data = hcat(data...)

    return data, sectionarray, forcecol, dispcol, timecol
end

"""
    AFMfileload(filedir::String, test_type::String; visco::Bool = true, cpfind::String = "hertz", param::Float64 = NaN64)

Load data from a JPK AFM (plaintext) file format into an AFMData struct.
AFMData struct provides the basis for subsequent high level operations
within RHEOS. test_type is either "strlx" for a stress-relaxation test
(strain controlled) or "creep" for a creep test (stress controlled).
If `visco = true`, only approach and hold sections are included. If
`visco = false` then only approach seciton included, which may be useful
for elastic analysis only. Retraction portion of data is always discarded.

`cpfind` determines contact point detection method used, if any. Options are:

- `hertz`: Used by default, fits the approach section of the data to an elastic
           Hertz model and then infers the contact point.

- `threshold`: Takes contact point as element after force exceeds `param`.

- `none`: No contact point detection algortihm is used.

# Example

```jldoctest
# directory path to the file
filedir = "../data/afmData1.txt"

# load the data into AFMData struct
dataforprocessing = AFMfileload(filedir, "strlx")
```
"""
function AFMfileload(filedir::String, test_type::String; visco::Bool = true, cpfind::String = "hertz", R::Float64 = NaN64, param::Float64 = NaN64)

    # dictionary for referencing contact model functions
    contactmodels = Dict( "hertz" => contact_hertz,
                          "threshold" => contact_threshold,
                          "none" => contact_none )

    # check that params provided where necessary
    if cpfind == "threshold"
        @assert !isnan(param) "If force threshold method is used, a force threshold must be provided using param keyword argument."
    elseif cpfind == "hertz"
        @assert !isnan(R) "If Hertz contact point method detection is used, radius of spherical indenter must be provided."
        param = R
    end

    # get data from JPK formatted file using functions built above,
    # add ability to parse other formats in future
    (data, sec, fcol, δcol, tcol) = AFMdataget(filedir)

    # get contact contact point info using only approach section
    cp_index = contactmodels[cpfind](data[sec[1], fcol], data[sec[1], δcol]; _param = param)

    # get data
    if visco
        f = vcat(data[sec[1], fcol][cp_index:end], data[sec[2], fcol])
        δ = vcat(data[sec[1], δcol][cp_index:end], data[sec[2], δcol])
        t = vcat(data[sec[1], tcol][cp_index:end], data[sec[2], tcol])
    else
        f = data[sec[1], fcol][cp_index:end]
        δ = data[sec[1], δcol][cp_index:end]
    end

    # regularise data
    f = f - minimum(f)
    δ = -δ - minimum(-δ)
    if visco
        t = t - minimum(t)
    end

    # send to constructor and return
    AFMData(f, δ, t, R, test_type, filedir)

end

"""
    RheologyType = Union{RheologyData, AFMData}

A superset of the two possible data types accepted by preprocessing, processing and postprocessing functions.
"""
RheologyType = Union{RheologyData, AFMData}
