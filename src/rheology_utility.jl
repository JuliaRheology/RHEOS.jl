#!/usr/bin/env julia

"""
    RheologyData(stress::Array{Float64,1}, strain::Array{Float64,1}, time::Array{Float64,1}, test_type::String)

RheologyData mutable struct used for high level interaction with RHEOS
preprocessing and fitting functions. Initialise an instance directly or
indirectly. If data is in three column, comma seperated CSV file then
fileload function can be used, which calls the RheologyData constructor.
If not, load data in according to format and call RheologyData constructor.
"""
mutable struct RheologyData

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

    # for initial loading of full dataset (measured + prescribed)
    RheologyData(σ::Array{Float64,1}, ϵ::Array{Float64,1}, t::Array{Float64,1}, ) = completeconstruct!(new(σ, ϵ, t))

    RheologyData(σ::Array{Float64,1}, ϵ::Array{Float64,1}, t::Array{Float64,1}, filedir::String) = completeconstruct!(new(σ, ϵ, t), filedir)

    # inner constructor for incomplete data; data_controlled should generally
    # be controlled variable (stress for creep, strain for strlx).
    RheologyData(data_controlled::Array{Float64,1}, t::Array{Float64,1}, test_type::String) =
                partialconstruct!(new(filedir, false, "constant", test_type), data_controlled, t, test_type)

end

"""
Inner constructor for RheologyData struct, providing gradients of stress
and strain. If stress/strain arrays have NaN values at the beginning (some datasets
have 1 or 2 samples of NaN at beginning) then deletes these and starts at the first
non-NaN sample, also readjusts time start to t = 0 to account for NaNs and and
negative time values at beginning of data recording.
"""
function completeconstruct!(self::RheologyData)

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

    # Check if time vector is equally spaced
    diff = round.(self.t[2:end]-self.t[1:end-1],4)
    check = any(x->x!=diff[1], diff)
    if check == true
       self.sampling = "variable"
    else
       self.sampling = "constant"
    end

    # readjust time to account for NaN movement and/or negative time values
    self.t = self.t - minimum(self.t)

    # initialise empty array to record operations applied during preprocessing
    self.log = []

    # initialise array of sybols for use if resampling/smoothing all numerical data
    self.numericdata = [:σ, :ϵ, :t, :measured, :dcontrolled]

    # return class with all fields initialised
    self
end

function completeconstruct!(self::RheologyData, filedir::String)

    self = completeconstruct!(self)

    self.log.append("file loaded from $filedir")

end

"""
Inner constructor completion function for case when only partial data is provided.
"""
function partialconstruct!(self::RheologyData, data_controlled::Array{Float64,1}, t::Array{Float64,1}, test_type::String)

    # define as local so it can be accessed in subsequent scopes
    local newstartingval::T where T<:Integer

    # test for NaNs
    for i in 1:length(data_controlled)
        if !isnan(data_controlled[i])
            newstartingval = i
            break
        end
    end

    # adjust starting point accordingly to remove NaNs
    data_controlled = data_controlled[newstartingval:end]
    t = t[newstartingval:end]

    # readjust time to account for NaN movement and/or negative time values
    t = t - minimum(t)

    # set time
    self.t = t

    # test dependent
    if test_type=="strlx"
        # set controlled variable as strain
        self.ϵ = data_controlled
        # derivative of stress WRT time
        self.dcontrolled = deriv(self.ϵ, self.t)

    elseif test_type=="creep"
        # set controlled variable as stress
        self.σ = data_controlled
        # derivative of strain WRT time
        self.dcontrolled = deriv(self.σ, self.t)
    end

    # initialise empty dictionary for model fit results
    self.fittedmodels = Dict()

    # initialise empty array for operations applied
    self.appliedops = []

    # initialise array of sybols for use if resampling/smoothing all numerical data
    self.numericdata = [:σ, :ϵ, :t, :measured, :dcontrolled]

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
    data = RheologyData(data[cols[1]], data[cols[2]], data[cols[3]], filedir)

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
