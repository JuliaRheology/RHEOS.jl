
"""
    JPKgetlines(filedir)
Identify beginning and end of AFM data sections, fancy names
line and version number of JPK data processing software.
Used in JPKconvertRHEOS.
"""
function JPKgetlines(filedir::String)

    # initialise array for storing line numbers of section beginnings and endings
    sectionbreaks = Int32[]

    # initialise other local variables
    fancynamesraw = ""
    vernumraw = ""
    lastlinenum = -1

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

    blocklen = [0 for i in 1:numsections]
    blockmod = [0 for i in 1:numsections]
    blockcount = 0

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
    sectionArray  = [0:1 for i in 1:(length(sectionbreaksproc)-1)]

    for i in 1:(length(sectionbreaksproc)-1)
        sectionArray[i] = sectionbreaksproc[i]:(sectionbreaksproc[i+1]-1)
    end

    return sectionArray, fancynamesraw, vernumraw
end

"""
    JPKversionnumstrip(vernumraw::String)
Used in JPKconvertRHEOS pipeline - clean up version number line.
"""
function JPKversionnumstrip(vernumraw::String)

    return String(split(vernumraw)[3])
end

"""
    JPKfancynames(fancynamesraw::String, vernum::String)
Used in JPKconvertRHEOS pipeline - clean up fancy names string so column
titles are in array format.
"""
function JPKfancynamestrip(fancynamesraw::String, vernum::String)

    fancynames_split = split(fancynamesraw,"\"")

    fancynames1 = deleteat!(fancynames_split, 1)

    fancynames = fancynames1[1:2:end]

    forcecol = -1;
    timecol = -1;
    dispcol = -1;

    # use fancynames and vernum to get column numbers of relevant data
    forcecol = findall(fancynames .== "Vertical Deflection")[1]

    # if no time column present then return timecol = 0
    if length(findall(fancynames .== "Series Time")) == 1
        timecol = findall(fancynames .== "Series Time")[1]
    else
        timecol = 0
    end

    _vernum = parse(Int16, vernum[1])

    if _vernum < 6
        dispcol = findall(fancynames .== "Tip-Sample Separation")[1]
    elseif _vernum >= 6
        dispcol = findall(fancynames .== "Vertical Tip Position")[1]
    end
    return forcecol, dispcol, timecol
end

"""
    JPKconverRHEOS(filedir::String)
convenience function for loading in JPK AFM data into RHEOS with relevant metadata
and split in to sections.
"""
function JPKconvertRHEOS(filedir::String)

    # get section break line numbers, column names and version info
    (sectionarray, fancynamesraw, vernumraw) = JPKgetlines(filedir)

    vernum = JPKversionnumstrip(vernumraw)

    (forcecol, dispcol, timecol) = JPKfancynamestrip(fancynamesraw, vernum)

    # get sectioned data with correct column names
    data = readdlm(filedir, comments=true, comment_char='#')

    time_force_disp = hcat(data[:,timecol],data[:,forcecol],data[:,dispcol])

    return time_force_disp, sectionarray
end
