# File I/O

## CSV Import/Export

In the [Fitting Data](@ref) section we saw how to import data from csv files. As a brief overview of that functionality, the 3 lines below demonstrate how to import stress/strain/time data, partial strain/time data, and G'/G''/frequency data. The columns can be in any order but need to specified by the first argument array of strings shown below. Note that if you want to try and run this code exactly as shown below, the data files used in these examples are stored in the 'RHEOS/examples' directory.
```
data = importdata(["stress","strain", "time"], "DataRelaxation.csv")

data_incomplete = importdata(["strain", "time"], "DataIncomplete.csv")

data_dynamic = importdata(["Gp","Gpp", "Frequency"], "FrequencyData.csv")
```
Where the first two function calls return [`RheologyData`](@ref) objects and the last function call returns a [`RheologyDynamic`](@ref) object.

If you want to analyse or plot your data in software other than Julia you will likely want to export it to a CSV file. To export [`RheologyData`](@ref) and [`RheologyDynamic`](@ref) objects to CSV files we can use the [`exportdata`](@ref) function. For the two complete data-sets we imported above, we can export them into new files in the following way.
```
exportdata(data, "exported_data")

exportdata(data_dynamic, "exported_dynamic_data")
```
The '.csv' extension will automatically be added (and can be modified if needs be through use of a keyword argument). The delimiter can also be changed if necessary. See [API](@ref) section for information on this.

## Native Julia data files (JLD2)

When using RHEOS and doing all subsequent analysis and plotting in Julia, it is convenient to be able to store native Julia objects to disk so that they can be loaded in to subsequent Julia sessions. RHEOS provides two pairs of convenience functions to facilitate this: [`savedata`](@ref) and [`loaddata`](@ref), [`savemodel`](@ref) and [`loadmodel`](@ref). Using our data imported above we can demonstrate saving and loading data (it is exactly the same process for both RHEOS data types).
```
savedata(data, "imported_data")

loaded_data = loaddata("imported_data.jld2")
```
The [`savedata`](@ref) function will automatically append a '.jld2' to the filename when saving but this can be modified if necessary through a keyword argument. [`savedata`](@ref) and [`loaddata`](@ref) work on both [`RheologyData`](@ref) and  [`RheologyDynamic`](@ref) objects.

Lastly, let's create a Standard Linear Solid model object with default parameters, save that to disk and load it again.
```
model = SLS()

savemodel(model, "saved_model")

loaded_model = loadmodel("saved_model.jld2")
```