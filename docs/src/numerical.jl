# # Numerical Aspects
#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/numerical.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/numerical.ipynb)

#md # !!! note "Note"
#md #     This section is still under development.

using RHEOS
#-

# ## Constant vs. Variable Sample Rate
# #### Weight indices rather than using a variable sample rate
# There may be situations when it is appropriate to focus the model fitting process on particular regions of the data. This may occur, for example, when the small time scale loading region is of particular interest. To achieve this, we could send data to fit that is sampled variably, with a greater sample density near the region of interest. However, the variable sampling rate means that the method RHEOS uses to compute the model predictions is slower and less accurate.

# A more optimal solution is to send constant sample rate data to the fitting function, but ask RHEOS to modify the final cost function to focus on specific points. In this way, RHEOS will be able to benefit from the much faster and more accurate constant sample rate model simulation methods.

# In concrete terms, consider the following array of time samples:

t = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5];

# We can send another array to the fitting process, `weightings`, that tells RHEOS to weight the first three time samples as doubly important, and the subsequent samples as half as important. The array we send would look like this:

weightings = [1, 1, 2, 2, 3, 3, 4, 6];

# In words, the elements in the time region we are particularly concerned about (0 ≤ t ≤ 0.2) are considered twice in the final fitting cost function, whereas we only consider every other element for times after 0.2 seconds.

# For practical situations, when working with data embedded in [`RheoTimeData`](@ref) structs, RHEOS has a convenience function [`indexweight`](@ref) which makes it easy to generate a weightings array that can then be passed to either [`modelfit`](@ref) or [`modelstepfit`](@ref). The full description of [`indexweight`](@ref) can be found in the [API](@ref) section. Below, we include a typical example of its use, based on the same data used in [Fitting and Predicting - Time Data](@ref). 

## As in original example, load in the data
data = importcsv("assets/data_time.csv", t_col = 1, ϵ_col = 2, σ_col = 3)
rheotimedatatype(data)
#-

# Now instead of going straight to fitting, we will construct a weightings array using [`indexweight`](@ref).

## construct weighted array
weightings = indexweight(data; elperiods = [2, 3, -2], time_boundaries = [0.0, 0.5, 1.0, data.t[end]])

## fit model with weighted array
modelfit(data, Maxwell, strain_imposed, weights = weightings)

# We note that the results are slightly different from the original example in [Fitting and Predicting - Time Data](@ref) due to the element weighting provided to the fitting function.
