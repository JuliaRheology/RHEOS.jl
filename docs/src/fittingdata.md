# Fitting Data

This page is a tutorial on how to fit viscoelastic models to data using RHEOS. If you want to try out the code below, it can all be run from the Julia REPL but **note that the importing data functions will only work if you are using the 'RHEOS/examples' folder as your working directory as that's where the example data files are stored**.

## Tensile/Compression Data

This section is for standard viscoelastic tensile or compression tests where time, stress and strain data are available.

First, we need to load in RHEOS
```
using RHEOS
```
RHEOS has a convenience function for importing data from CSV files. The default column delimiter is ',' but an alternative can be specified as a keyword argument. The row delimiter is a newline character ('\n'). For standard viscoelastic testing data RHEOS expects either stress, strain and time data, just stress and time, or just strain and time. The order of the columns is specified as the first argument in the function `importdata`. The second argument is the directory of the file, as shown below.
```
data = importdata(["stress","strain", "time"], "DataComplete.csv")
```
In this tutorial, our data file `DataComplete.csv` is in the same directory as our Julia script so we can just use its relative directory. 

Let's fit a Standard Linear Solid viscoelastic model via its relaxation modulus, G, as our data is from a stress relaxation test. The first argument is our data, the second argument tells RHEOS which model to fit and the final argument tells RHEOS whether to fit the model using a relaxation modulus (:G) or creep modulus (:J).
```
fitted_SLS_model = modelfit(data, SLS(), :G)
```
Our first fitted model is not contained in the `fitted_SLS_model` variable.

Next, we'll fit a fractional Standard Linear Solid model (the only difference from the above model is that the dash-pot is replaced by a spring-pot). This time we'll also add upper and lower bounds on the model parameters. This is highly recommended for fractional models in particular as values less than 0 or greater than 1 for the spring-pot parameter are unphysical and can cause errors in the Mittag-Leffler function used.
```
lb = [0.1, 0.01, 0.1, 0.1]
ub = [Inf, 0.99, Inf, Inf]
fitted_fractSLS_model = modelfit(data, FractionalSLS(), :G; lo=lb, hi=ub)
```
Note the two keyword arguments used -- `lo` and `hi` for the lower and upper parameter boundaries respectively. The special argument `Inf` for the three of the parameters' upper bounds represent a type of infinity such that the parameters can be as large as required by the optimisation algorithm.

For a full list of keyword arguments and features of the `modelfit` function, see the relevant part of the API section. Models included in RHEOS are also listed in the API section, and discussed in more detail in the Models section.

## Complex Modulus and Frequency Data

RHEOS can also fit models to dynamic mechanical analysis data from oscillatory tests. The `importdata` function can again be used here but with the column names `Gp` for the storage modulus, `Gpp` for the loss modulus and `frequency` for the frequency column. RHEOS will detect the `frequency` string and know to load the data into a `RheologyDynamic` data type. Assuming RHEOS has already been imported, let's load in our data file:
```
data = importdata(["Gp","Gpp", "Frequency"], "FrequencyData.csv")
```
As before, we'll try and fit the data to a Standard Linear Solid model -- but this time we need to use the `dynamicmodelfit` function.
```
fitted_SLS_model = dynamicmodelfit(data, SLS())
```
Note that we do not have to specify whether we are fitting via the creep or relaxation modulus as RHEOS always fits dynamic data to the complex (dynamic) modulus. Let's also try to fit the data to a fractional Standard Linear Solid model:
```
lb = [0.1, 0.01, 0.1, 0.1]
ub = [Inf, 0.99, Inf, Inf]
fitted_fractSLS_model = dynamicmodelfit(data, FractionalSLS(); lo=lb, hi=ub)
```
For more information on the `dynamicmodelfit` function, see the API section.
