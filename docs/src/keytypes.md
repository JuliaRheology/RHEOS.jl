# Rheology Data

### Data Structure
```@docs
RheologyData(σ::Array{Float64,1}, ϵ::Array{Float64,1}, t::Array{Float64,1}, sampling::String, log::Array{String,1})
```

```@docs
RheologyData(colnames::Array{String,1}, data1::Array{Float64,1}, data2::Array{Float64,1}[, data3::Array{Float64,1}; filedir::String="none", log::Array{String,1}=Array{String}(0)])::RheologyData
```

### Model Structure
```@docs
RheologyModel(name::Function, parameters::Array{Float64,1}[, log::Array{String,1}])
```

### Generated Data Structure
```@docs
RheologyArtificial(data::Array{Float64,1}, t::Array{Float64,1}, stepsize::Float64, log::Array{String,1})
```