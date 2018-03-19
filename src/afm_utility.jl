#!/usr/bin/env julia

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

    # original data
    f::Array{Float64,1}
    δ::Array{Float64,1}
    t::Array{Float64,1}

    # force and displacement to appropriate power and multiplied by prefactor,
    # numerically differentiated WRT to time
    dfᵩ::Array{Float64,1}
    dδᵩ::Array{Float64,1}

    # models
    fittedmodels::Dict

    # operations applied, stores history of which functions (including arguments)
    appliedops::Array{String,1}

    # for initial loading of full dataset (measured + prescribed)
    AFMData(f::Array{Float64,1}, δ::Array{Float64,1}, t::Array{Float64,1}, test_type::String, filedir::String; visco::Bool = true) 
            = AFMconstruct!(new(filedir, false, "constant", test_type, f, δ, t), f, δ, t, visco)

    # # inner constructor for resampled data
    # AFMData(_filedir::String, _insight::Bool, _sampling::String, _test_type::String,
    #              _f::Array{Float64,1}, _δ::Array{Float64,1}, _t::Array{Float64,1},
    #              _df::Array{Float64,1}, _dδ::Array{Float64,1}, _appliedops::Array{String,1}) =
    #              new(_filedir, _insight, _sampling, _test_type,
    #              _f, _δ, _t, _df, _dδ, Dict(), _appliedops)

    # # inner constructor for incomplete data; data_part should generally
    # # be controlled variable (stress for creep, strain for strlx).
    # AFMData(data_part::Array{Float64,1}, t::Array{Float64,1}, test_type::String, filedir::String) =
    #             partialconstruct!(new(filedir, false, "constant", test_type), data_part, t, test_type)

end

"""
Inner constructor for AFMData struct, providing gradients of stress and strain. `visco` argument
determines whether approach and hold regions are used (`visco=true`) or just approach (`visco=false`).
The second case is useful if only elastic analysis is required.
"""
function AFMconstruct!(self::AFMData, f::Array{Float64,1}, δ::Array{Float64,1}, t::Array{Float64,1}, visco)

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

    # initialise empty array to record operations applied during preprocessing
    self.appliedops = []

    # return class with all fields initialised
    self
end

"""
Inner constructor completion function for case when only partial data is provided.
"""
function partialconstruct!(self::AFMData, data_part::Array{Float64,1}, t::Array{Float64,1}, test_type::String)

    # define as local so it can be accessed in subsequent scopes
    local newstartingval::T where T<:Integer

    # test for NaNs
    for i in 1:length(data_part)
        if !isnan(data_part[i])
            newstartingval = i
            break
        end
    end

    # adjust starting point accordingly to remove NaNs
    data_part = data_part[newstartingval:end]
    t = t[newstartingval:end]

    # readjust time to account for NaN movement and/or negative time values
    t = t - minimum(t)

    # set time
    self.t = t

    # test dependent
    if test_type=="strlx"
        # set controlled variable as strain
        self.ϵ = data_part
        # derivative of stress WRT time
        self.dϵ = deriv(self.ϵ, self.t)

    elseif test_type=="creep"
        # set controlled variable as stress
        self.σ = data_part
        # derivative of strain WRT time
        self.dσ = deriv(self.σ, self.t)
    end

    # initialise empty dictionary for model fit results
    self.fittedmodels = Dict()

    # initialise empty array for operations applied
    self.appliedops = []

    # return class with all fields initialised
    self
end


"""
    AFMfileinspect(filedir::String, test_type::String)


"""
function AFMgetlines(filedir::String, test_type::String)::AFMData
    #function identifies beginning and end of AFM data sections,
    #returns fancy names line,
    #and version number of JPK data processing software.

    #initialise array for storing line numbers of section beginnings and endings
    sectionBreaks = Int32[]

    #initialise other local variables
    local fancynamesraw::String
    local vernumraw::String
    local lastLineNum::Int32

    #iterate through lines
    open(filedir) do filJ
        for (i, line) in enumerate(eachline(filJ))

            if startswith(line, "# units")
                append!(sectionBreaks, i+1)

            elseif line == ""
                append!(sectionBreaks, i)

            elseif startswith(line, "# fancyNames:")
                #Get names of columns
                fancynamesraw = line

            elseif startswith(line, "# data-description.modification-software:")
                #get software version number
                vernumraw = line

            elseif eof(filJ)
                lastLineNum = i+1

            end
        end
        #append last line
        append!(sectionBreaks, lastLineNum)
    end

    #get section starts and stops, i.e. shift section break line numbers to account comment lines not read by readtable
    numSections = Int(length(sectionBreaks)/2)

    local blockLen = [0 for i in 1:numSections]
    local blockMod = [0 for i in 1:numSections]
    local blockCount = 0

    for num in (1:numSections)
        #get lengths of each section
        blockLen[num] = sectionBreaks[2:2:end][num] - sectionBreaks[1:2:end][num]
        #accumulate lengths
        blockCount += blockLen[num]
        #store and account for extra values from boundaries
        blockMod[num] = blockCount - num + 1
    end

    #add first line
    sectionBreaksProc = prepend!(blockMod, [1])

    # println(sectionBreaksProc)

    #get sections as unit ranges for ease of use
    local sectionArray  = [0:1 for i in 1:(length(sectionBreaksProc)-1)]

    for i in 1:(length(sectionBreaksProc)-1)
        sectionArray[i] = sectionBreaksProc[i]:(sectionBreaksProc[i+1]-1)
    end

    return sectionArray, fancynamesraw, vernumraw
end

function versionNumStrip(vernumraw::String)
    #clean up version number line
    return split(vernumraw)[3]
end

function fancynamestrip(fancynamesraw::String, vernum::SubString{String})
    #clean up fancy names string so column titles are in array format
    fancyNamesSplit = split(fancynamesraw,"\"")

    fancyNames1 = deleteat!(fancyNamesSplit, 1)

    fancyNames = fancyNames1[1:2:end]

    local forceCol::Int64
    local timeCol::Int64
    local dispCol::Int64

    #use fancyNames and vernum to get column numbers of relevant data
    forceCol = find(fancyNames .== "Vertical Deflection")[1]

    timeCol = find(fancyNames .== "Series Time")[1]

    if Int(vernum[5]) < 6
        dispCol = find(fancyNames .== "Tip-Sample Separation")[1]
    elseif Int(vernum[5]) >= 6
        dispCol = find(fancyNames .== "Vertical Tip Position")[1]
    end
    return forceCol, dispCol, timeCol
end

function AFMdataprocess(filedir::String)
    #function calls other functions to get column + version info, and data.

    #get section break line numbers, column names and version info
    (sectionArray, fancynamesraw, vernumraw) = AFMgetlines(filedir)

    vernum = versionNumStrip(vernumraw)

    (forceCol, dispCol, timeCol) = fancynamestrip(fancynamesraw, vernum)

    #get sectioned data with correct column names
    (data, header) = uCSV.read(filedir; delim=' ', comment='#', types=[Float64,Float64,Float64])

    data = hcat(data...)

    return data, sectionArray, forceCol, dispCol, timeCol
end

"""
    AFMfileload(filedir::String, test_type::String; visco::Bool = true, cpfind::String = "hertz", cpthresh::Float64 = NaN64)

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

- `threshold`: Takes contact point as element after force exceeds `cpthresh`.

- `none`: No contact point detection algortihm is used.

# Example

```jldoctest
# directory path to the file
filedir = "../data/afmData1.txt"

# load the data into AFMData struct
dataforprocessing = AFMfileload(filedir, "strlx")
```

"""
function AFMfileload(filedir::String, test_type::String; visco::Bool = true, cpfind::String = "hertz", cpthresh::Float64 = NaN64)

    # check that cpthresh provided
    if cpfind=="threshold"
        @assert !isnan(cpthresh) "If force threshold method is used, a cpthresh must be provided."
    end

    # call to get data from file using functions built above

    # add calls to cp methods
    # using only approach section

    # concat sections if visco, don't if not and send data to constructor

end