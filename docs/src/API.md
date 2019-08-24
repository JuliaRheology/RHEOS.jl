# API

## Main RHEOS Structs
```@docs
RheoTimeData
RheoFreqData
RheoModelClass
RheoModel
```

## Sampling and Filtering Functions
```@docs
resample
cutting
smooth
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
freeze_params
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