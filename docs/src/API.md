# API
```@meta
EditURL = "https://github.com/JuliaRheology/RHEOS.jl/edit/master/docs/src/API.md"
```

## Main RHEOS data structures
```@docs
RheoTimeData
RheoFreqData
gettime
getstress
getstrain
getfreq
getstorage
getloss
```

## Moduli/compliance functions
```@docs
RheoModelClass
RheoModel
getparams
freezeparams
relaxmod
creepcomp
storagemod
lossmod
dynamicmod
```

## Sampling and Filtering Functions
```@docs
resample
indexweight
cutting
smooth
onlytime
onlystrain
onlystress
onlyfreq
extract
```

## Fitting and Predicting Functions
```@docs
modelfit
modelpredict
modelstepfit
modelsteppredict
dynamicmodelfit
dynamicmodelpredict
```

## Data Generation Functions
```@docs
timeline
strainfunction
stressfunction
hstep
ramp
stairs
square
sawtooth
triangle
frequencyspec
```

## Data IO
```@docs
importcsv
exportcsv
```
