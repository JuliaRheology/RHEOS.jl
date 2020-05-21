
"""
    JPKgetlines(filedir)
Identify beginning and end of AFM data sections, fancy names
line and version number of JPK data processing software and x and y position of the curve in a force map.
Used in JPKconvertRHEOS.
"""
function JPKgetlines(filedir::String)

    # initialise array for storing line numbers of section beginnings and endings
    sectionbreaks = Int32[]

    # initialise other local variables
    fancynamesraw = ""
    vernumraw = ""
    xposition = ""
    yposition = ""
    segments = Array{String}(undef, 0)
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
            elseif startswith(line, "# xPosition:")
                xposition = line
            elseif startswith(line, "# yPosition:")
                yposition = line  
            elseif startswith(line, "# segment:") 
                push!(segments, line)
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

    return sectionArray, fancynamesraw, vernumraw, xposition, yposition, segments
end

"""
    JPKversionnumstrip(vernumraw::String)
Used in JPKconvertRHEOS pipeline - clean up version number line.
"""
function JPKversionnumstrip(vernumraw::String)

    return String(split(vernumraw)[3])
end

"""
    JPKposition(xposition::String, yposition::String)
Used in JPKconvertRHEOS pipeline - obtain the x and y position of the pixel from a force map.
"""
function JPKposition(xposition::String, yposition::String)
    x = split(xposition)[3];
    y = split(yposition)[3];
    return parse(RheoFloat, x), parse(RheoFloat, y)
end

"""
    JPKsegments(segments::Array{String})
Used in JPKconvertRHEOS pipeline - obtain the description of each segment.
"""
function JPKsegments(segments::Array{String})
    list_segments = Array{String}(undef, 0);

    for i = 1:1:size(segments,1)
        push!(list_segments,segments[i][12:end])
    end
    return list_segments
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


function IndexToExtract(list_segments::Array{String}, sections::Array{String})

    j = []

    for i = 1:1:size(sections,1)
        index = findall(x -> x == sections[i], list_segments)
        
        if !isempty(index)
            push!(j,index[1])
        end
    end
    j = sort(j)

    return j

end

"""
    JPKconverRHEOS(filedir::String)
convenience function for loading in JPK AFM data into RHEOS with relevant metadata
and split in to sections.
"""
function importJPK(filepath::String, interface::Interface; sections::Array{String} = nothing, comment = "Imported from JPK txt file", savelog = true)

    # get section break line numbers, column names and version info, and x and y positions of the pixel, and the names of the segments
    (sectionarray, fancynamesraw, vernumraw, xposition, yposition, segments) = JPKgetlines(filepath);

    vernum = JPKversionnumstrip(vernumraw);

    x_coord, y_coord = JPKposition(xposition, yposition);

    list_segments = JPKsegments(segments);

    (f_col, d_col, t_col) = JPKfancynamestrip(fancynamesraw, vernum)

    # get sectioned data with correct column names
    data_raw = readdlm(filepath, comments=true, comment_char='#')

                                # TO CHECK: displacement correction!

    # extract only the desired sections (["extend", "hold", "retract"])
    
    if sections == nothing

        disp = -data_raw[:,d_col] 
        min_disp = findmin(disp)[1]

        data = RheoTimeData(t = data_raw[:, t_col], ϵ = disp .- min_disp, σ = data_raw[:, f_col])

    else

        index = IndexToExtract(list_segments, sections)

        i_start = sectionarray[index[1]][1]
        i_end = sectionarray[index[end]][end]

        disp = -data_raw[i_start:i_end,d_col]
 
        min_disp = findmin(disp)[1]

        data = RheoTimeData(t = data_raw[i_start:i_end, t_col], ϵ = disp .- min_disp[1], σ = data_raw[i_start:i_end, f_col])

    end

    
    # convert force and displacement into stress and strain
    ϵ,σ = interface.to_ϵσ(data.ϵ, data.σ)

    # create log 
    log = if savelog
			 info = (comment=comment, folder=pwd(), stats=(t_min=data.t[1],t_max=data.t[end], n_sample=size(data.t[:])))
			 RheoLogItem( (type=:source, funct=:importJPK, params=(filepath=filepath, interface=interface), sections = sections, comment = comment), info )
		   else
			 RheoLogItem(type = nothing, info = "nothing")
		   end

   
    # save extra information into named tuple
    extra_info = (x = x_coord, y = y_coord, sections = sectionarray, sections_names = list_segments)


    return RheoTimeData(σ, ϵ, data.t, [log]), extra_info
end
