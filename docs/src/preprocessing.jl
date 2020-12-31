# # Preprocessing of Data

# RHEOS offers several functions for sampling and filtering data; this page is intended to be a brief tutorial of their use. For detailed descriptions of functions and their optional arguments, see the [API](@ref) section.

using RHEOS
using PyPlot 

# ## Upsampling and Downsampling

# We generate a simple sinusoid with constant sampling rate (more details about data generation are discussed in the [Generating Data](@ref) section).

## Generate a sinusoidal data set
foo = timeline(t_start = 0, t_end = 10, step = 0.4)
foo = strainfunction(foo, t -> sin(t))
fig, ax = subplots(1, 1, figsize = (5, 5))
ax.plot(foo.t, foo.ϵ, "--", marker = "o", markersize = 8, color = "blue")
fig # for showing figure in docs

# To downsample the full data set by taking every nth sample, the function [`resample`](@ref) is defined with negative argument. Similarly, to increase the sample rate, the function [`resample`](@ref) is defined with positive argument.

#md # !!! note "Note"
#md #     When 1 or -1 is defined, the function returns the original [`RheoTimeData`](@ref), whilst 0 is not accepted as a valid argument.

#nb # Note: when 1 or -1 is defined, the function returns the original RheoTimeData, whilst 0 is not accepted as a valid argument.

## Downsample
foo_dsamp = resample(foo, -2)
## Upsample
foo_usamp = resample(foo, 2)

## Plotting
fig, ax = subplots(1,2, figsize=(10,5))
ax[1].set_title("Downsampling")
ax[1].plot(foo.t, foo.ϵ, "--", marker = "o", markersize = 8, color = "blue")
ax[1].plot(foo_dsamp.t, foo_dsamp.ϵ, "--", marker = "x", markersize = 10, markeredgewidth=2, color = "orange")
ax[2].set_title("Upsampling")
ax[2].plot(foo.t, foo.ϵ, "--", marker = "o", markersize = 8, color = "blue")
ax[2].plot(foo_usamp.t, foo_usamp.ϵ, "--", marker = "x", markersize = 10, markeredgewidth=2, color = "orange")
fig

# RHEOS also allows us to define time regions with different sampling rate. This requires the definition of the time boundaries to delimit where one sample rate finishes and another begins, it also includes the beginning and end points so the length of the time boundaries array argument should be one longer than the number of sample rates provided.

## Variable sampling: Downsampling region [0,5] and upsampling [5,10]
foo_samp = resample(foo, [-2,2], time_boundaries = [0.0, 5.0, 10.0])

## Plotting
fig, ax = subplots(1,1, figsize=(5,5))
ax.plot(foo.t, foo.ϵ, "--", marker = "o", markersize = 8, color = "blue")
ax.plot(foo_samp.t, foo_samp.ϵ, "--", marker = "x", markersize = 10, markeredgewidth = 2, color = "orange")
fig

# ## Cutting

# RHEOS provides a dedicated function to remove the data outside a specified time interval.
