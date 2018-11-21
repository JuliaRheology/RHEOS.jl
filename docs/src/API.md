# API

### Data Structure
```@docs
RheologyData(σ::Array{Float64,1}, ϵ::Array{Float64,1}, t::Array{Float64,1}, sampling::String, log::Array{String,1})
```

```@docs
RheologyData(colnames::Array{String,1}, data1::Array{Float64,1}, data2::Array{Float64,1}[, data3::Array{Float64,1}; filedir::String="none", log::Array{String,1}=Array{String}(0)])
```

### Model Structure
```@docs
RheologyModel(name::Function, parameters::Array{Float64,1}[, log::Array{String,1}])
```

### Generated Data Structure
```@docs
RheologyArtificial(data::Array{Float64,1}, t::Array{Float64,1}, stepsize::Float64, log::Array{String,1})
```

```@docs
modelfit(data::RheologyData, modulus::Function[, p0::Array{Float64,1}, lo::Array{Float64,1}, hi::Array{Float64,1}; verbose::Bool = false])
```

```@docs
modelpredict(data::RheologyData, model::RheologyModel)
```

```@docs
modelfit(data::RheologyData, modulus::Function[, p0::Array{Float64,1}, lo::Array{Float64,1}, hi::Array{Float64,1}; verbose::Bool = false])
```

```@docs
modelpredict(data::RheologyData, model::RheologyModel)
```

```@docs
var_resample(self::RheologyData, refvar::Symbol, pcntdownsample::Float64; mapback::Bool = false)
```

```@docs
downsample(self::RheologyData, time_boundaries::Array{Float64,1}, elperiods::Array{Int64,1})

```

```@docs
fixed_resample(self::RheologyData, time_boundaries::Array{Float64,1}, elperiods::Array{Int64,1}, direction::Array{String,1})
```

```@docs
smooth(self::RheologyData, τ::Float64)
```

```@docs
function mapbackdata(self_new::RheologyData, self_original::RheologyData)
```

```@docs
savedata(self::RheologyData; filedir::String = "", ext = "_RheologyData.jld")
```

```@docs
loaddata(filedir::String)
```

```@docs
savemodel(self::RheologyModel; filedir::String = "", ext = "")
```

```@docs
loadmodel(filedir::String)
```

```@docs
exportdata(self::RheologyData; filedir::String = "", ext = "_mod.csv")
```