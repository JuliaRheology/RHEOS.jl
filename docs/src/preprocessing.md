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
