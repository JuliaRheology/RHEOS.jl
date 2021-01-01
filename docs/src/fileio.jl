# # File I/O
#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/fileio.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/fileio.ipynb)

# ## Import Data

using RHEOS

# RHEOS has a convenience function for importing data from CSV files: [`importcsv`](@ref). The default column delimiter is ',' but an alternative can be specified as a keyword argument. The row delimiter is a newline character ('\n'). For standard time-domain viscoelastic testing data RHEOS expects either stress, strain and time data, just stress and time, or just strain and time. Arguments must be identified by providing the number of the column in which they are contained. The function returns a [`RheoTimeData`](@ref) object.

## Import file
data_1 = importcsv("assets/data_time.csv", t_col=1, ϵ_col=2, σ_col=3)
## Check data type
RheoTimeDataType(data_1)

## Import incomplete data
data_2 = importcsv("assets/data_time.csv", t_col=1, ϵ_col=2)
## Check data type
RheoTimeDataType(data_2)

# The function [`importcsv`](@ref) can also be used to import frequency, storage modulus, and loss modulus data (as a complete set). In this case, the function returns a [`RheoFreqData`](@ref) object.

## Import file
data_f = importcsv("assets/data_freq.csv", ω_col=1, Gp_col=2, Gpp_col=3)
## Check data type
RheoFreqDataType(data_f)

# ## Export Data

# If you want to analyse or plot your data in software other than Julia you will likely want to export it to a CSV file. To export [`RheoTimeData`](@ref) and [`RheoFreqData`](@ref) objects to CSV files we can use the [`exportcsv`](@ref) function. For the two complete data-sets we imported above, we can export them into new files in the following way. As with [`importcsv`](@ref), the order of the columns can be specified by the user. 

## Export file
exportcsv(data_1,"assets/export_timedata.csv")
exportcsv(data_f,"assets/export_frequdata.csv")
