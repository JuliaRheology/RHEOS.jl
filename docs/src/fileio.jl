# # File I/O
#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/simpleexamples.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/simpleexamples.ipynb)

# ## Import Data

using RHEOS

# RHEOS has a convenience function for importing data from CSV files: importcsv. The default column delimiter is ',' but an alternative can be specified as a keyword argument. The row delimiter is a newline character ('\n'). For standard time-domain viscoelastic testing data RHEOS expects either stress, strain and time data, just stress and time, or just strain and time. Arguments must be identified by providing the number of the column in which they are contained. The function returns a RheoTimeData object.