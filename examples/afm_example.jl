#!/usr/bin/env julia
include("../src/RHEOS.jl")
using RHEOS
using PyPlot

fdir = "../data/afmData1.txt"

test_fdirs = [  "C:\\Projects\\11_Louis\\GelsProject\\20170205_mixture41\\b3p2\\rs -2250--2017.02.05-13.54.58.932.txt"
                "C:\\Projects\\11_Louis\\GelsProject\\20161207_20161209_20161216_data\\DM41\\07b1p2_down\\RS -2000--2016.12.07-15.07.49.907.txt"
                "C:\\Projects\\11_Louis\\GelsProject\\20170219 plastic DM41 testing\\1-\\p01_4000nN_15sCreep_00s.txt"
                "C:\\Projects\\11_Louis\\GelsProject\\20161207_20161209_20161216_data\\DM41\\07b2p2_up\\RS -0500--2016.12.07-15.38.36.219.txt"
                "C:\\Projects\\11_Louis\\GelsProject\\20161207_20161209_20161216_data\\DM41\\16b2p2_down\\rs -0500--2016.12.16-12.15.16.169.txt"
                "C:\\Projects\\11_Louis\\GelsProject\\20161207_20161209_20161216_data\\DM41\\16b3p2_up\\rs -0500--2016.12.16-13.03.12.129.txt"
                "C:\\Projects\\11_Louis\\ALGAE CHAPTER\\2017 03 09 Data Processed\\Brown Fucus\\S1_1ASW_A1\\SP-150 RSP-0100-2017.02.01-14.02.24.523.txt"
                "C:\\Projects\\11_Louis\\ALGAE CHAPTER\\2017 03 09 Data Processed\\Brown Fucus\\S1_1ASW_A2\\SP-150 RSP-0100-2017.02.01-14.19.37.647.txt"
                "C:\\Projects\\11_Louis\\GelsProject\\20161207_20161209_20161216_data\\DM33\\09b2p2\\RS -0500--2016.12.09-12.33.50.934.txt"
                "C:\\Projects\\11_Louis\\GelsProject\\20161207_20161209_20161216_data\\DM33\\09b1p4\\RS -5750--2016.12.09-11.39.26.802.txt"
                "C:\\Projects\\11_Louis\\GelsProject\\20161207_20161209_20161216_data\\DM33\\09b2p3\\RS -5750--2016.12.09-12.44.58.763.txt"
                "C:\\Projects\\11_Louis\\GelsProject\\20161207_20161209_20161216_data\\DM33\\09b2p4\\RS -4500--2016.12.09-13.01.01.034.txt"
                ]

# data = AFMfileload(fdir, "strlx"; cpfind = "threshold", param = 1e-8)

data = AFMfileload(fdir, "strlx"; cpfind = "hertz", param = 25e-6)

# for fdir in test_fdirs
#     data = AFMfileload(fdir, "strlx"; cpfind = "hertz", param = 25e-6)
# end

# plot(data.t, data.f)
# plot(data.t, data.δ, "--")
# show()