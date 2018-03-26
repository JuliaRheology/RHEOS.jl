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

    # radius of sphere
    R::Float64

    # original data
    f::Array{Float64,1}
    δ::Array{Float64,1}
    t::Array{Float64,1}

    # measured variable to appropriate power
    cm_measured::Array{Float64,1}

    # controlled variable to appropriate power and multiplied by prefactor,
    # numerically differentiated WRT to time
    cm_dcontrolled ::Array{Float64,1}

    # models
    fittedmodels::Dict

    # operations applied, stores history of which functions (including arguments)
    appliedops::Array{String,1}

    # for initial loading of full dataset (measured + prescribed)
    AFMData(f::Array{Float64,1}, δ::Array{Float64,1}, t::Array{Float64,1}, R::Float64, test_type::String, filedir::String) = AFMconstruct!(new(filedir, false, "constant", test_type, R, f, δ, t))

    # inner constructor for resampled data
    AFMData(filedir::String, insight::Bool, sampling::String, test_type::String, R::Float64, f::Array{Float64,1}, δ::Array{Float64,1}, t::Array{Float64,1}, cm_measured::Array{Float64,1},
            cm_dcontrolled::Array{Float64,1}, fittedmodels::Dict, appliedops::Array{String,1}) = 
            new(filedir, insight, sampling, test_type, R, f, δ, t, cm_measured, cm_dcontrolled, fittedmodels, appliedops)

    # inner constructor for incomplete data; data_part should generally
    # be controlled variable (stress for creep, strain for strlx) TODO
    
end

"""
Inner constructor for AFMData struct, fills in fields as appropriate.
"""
function AFMconstruct!(self::AFMData)

    # assume spherical Hertz model
    prefactor = 8.0*sqrt(self.R)/3.0

    # define cm_measured and cm_dcontrolled as appropriate, REF Oyen paper
    if self.test_type == "strlx"

        self.cm_measured = self.f

        self.cm_dcontrolled = prefactor*deriv((self.δ).^(3/2), self.t)

    elseif self.test_type == "creep"

        self.cm_measured = (self.δ).^(3/2)

        self.cm_dcontrolled = deriv(self.f, self.t)/prefactor

    end

    # initialise empty dictionary for model fit results
    self.fittedmodels = Dict()

    # initialise empty array to record operations applied during preprocessing
    self.appliedops = String[]

    # return class with appropriate fields initialised
    self
end

# """
# Inner constructor completion function for case when only partial data is provided.
# """
# function partialconstruct!(self::AFMData, data_part::Array{Float64,1}, t::Array{Float64,1}, test_type::String)

#     # define as local so it can be accessed in subsequent scopes
#     local newstartingval::T where T<:Integer

#     # test for NaNs
#     for i in 1:length(data_part)
#         if !isnan(data_part[i])
#             newstartingval = i
#             break
#         end
#     end

#     # adjust starting point accordingly to remove NaNs
#     data_part = data_part[newstartingval:end]
#     t = t[newstartingval:end]

#     # readjust time to account for NaN movement and/or negative time values
#     t = t - minimum(t)

#     # set time
#     self.t = t

#     # test dependent
#     if test_type=="strlx"
#         # set controlled variable as strain
#         self.ϵ = data_part
#         # derivative of stress WRT time
#         self.dϵ = deriv(self.ϵ, self.t)

#     elseif test_type=="creep"
#         # set controlled variable as stress
#         self.σ = data_part
#         # derivative of strain WRT time
#         self.dσ = deriv(self.σ, self.t)
#     end

#     # initialise empty dictionary for model fit results
#     self.fittedmodels = Dict()

#     # initialise empty array for operations applied
#     self.appliedops = []

#     # return class with all fields initialised
#     self
# end

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