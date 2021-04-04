# # Preprocessing of Data
#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/preprocessing.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/preprocessing.ipynb)

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
#!nb fig #hide

# To downsample the full data set , the function [`resample`](@ref) is used with a `scale` argument less than 1. Similarly, to increase the sample rate, the `scale` argument should be greater than 1.

#md # !!! note "Note"
#md #     When scaling, the argument can be an integer, floating point, or rational type. The function contains greater functionality and versatility than shown here. See the [`resample`](@ref) docstring for more information.

#nb # Note: when scaling, the argument can be an integer, floating point, or rational type. The function contains greater functionality and versatility than shown here. See the `resample` docstring for more information.

## Downsample
foo_dsamp = resample(foo, scale = 1//2)
## Upsample
foo_usamp = resample(foo, scale = 2)

## Plotting
fig, ax = subplots(1, 2, figsize = (10, 5))
ax[1].set_title("Downsampling")
ax[1].plot(foo.t, foo.ϵ, "--", marker = "o", markersize = 8, color = "blue")
ax[1].plot(foo_dsamp.t, foo_dsamp.ϵ, "--", marker = "x", markersize = 10, markeredgewidth = 2, color = "orange")
ax[2].set_title("Upsampling")
ax[2].plot(foo.t, foo.ϵ, "--", marker = "o", markersize = 8, color = "blue")
ax[2].plot(foo_usamp.t, foo_usamp.ϵ, "--", marker = "x", markersize = 10, markeredgewidth = 2, color = "orange")
#!nb fig #hide

# ## Cutting

# RHEOS provides a dedicated function, [`cutting`](@ref), to remove the data outside a specified time interval.

foo_cut = cutting(foo, 2.0, 8.0)

## Plotting
fig, ax = subplots(1, 1, figsize = (5, 5))
ax.plot(foo.t, foo.ϵ, "--", marker = "o", markersize = 8, color = "blue")
ax.plot(foo_cut.t, foo_cut.ϵ, "--", marker = "x", markersize = 10, markeredgewidth = 2, color = "orange")
#!nb fig #hide

# ## Smoothing

# Lastly, RHEOS provides a smoothing function, [`smooth`](@ref). The first argument is the data to smooth and the second argument is the (very) approximate time scale of smoothing. (It uses Gaussian smoothing and can be thought of as a low pass filter for information occuring on time scales shorter than the 2nd argument). The padding can be changed using a keyword argument if desired, see [API](@ref) and [ImageFiltering.jl Documentation](https://juliaimages.github.io/ImageFiltering.jl/stable/) for more details.

foo_s = timeline(t_start = 0, t_end = 10, step = 0.02)
foo_s = strainfunction(foo_s, t -> 3*sin(t))
noise = strainfunction(foo_s, t -> rand())

foo_noisy = foo_s + noise
foo_smooth = smooth(foo_noisy, 1)

fig, ax = subplots(1, 1, figsize = (5, 5))
ax.plot(foo_noisy.t, foo_noisy.ϵ, color = "blue")
ax.plot(foo_smooth.t, foo_smooth.ϵ, color = "orange")
#!nb fig #hide
