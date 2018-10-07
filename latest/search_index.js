var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#RHEOS.jl-Documentation-1",
    "page": "Home",
    "title": "RHEOS.jl Documentation",
    "category": "section",
    "text": "RHEOS, an abbreviation of Rheology Open Source, is a Julia package which provides tools for analysis of rheology data and generation of new viscoelastic models. It features a large array of standard and fractional linear viscoelastic models - new ones can be easily added by the user. It can also be used as an educational tool to demonstrate the qualitative behaviour of various viscoelastic models. "
},

{
    "location": "index.html#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "Install Julia, version 1.0.1\nFrom Julia REPL, enter pkg mode by pressing ]\nRun the command add \"https://github.com/JuliaRheology/RHEOS.jl\""
},

{
    "location": "index.html#Program-Structure-1",
    "page": "Home",
    "title": "Program Structure",
    "category": "section",
    "text": ""
},

{
    "location": "keytypes.html#",
    "page": "Rheology Types",
    "title": "Rheology Types",
    "category": "page",
    "text": ""
},

{
    "location": "keytypes.html#Rheology-Data-1",
    "page": "Rheology Types",
    "title": "Rheology Data",
    "category": "section",
    "text": ""
},

{
    "location": "keytypes.html#Data-Structure-1",
    "page": "Rheology Types",
    "title": "Data Structure",
    "category": "section",
    "text": "RheologyData(σ::Array{Float64,1}, ϵ::Array{Float64,1}, t::Array{Float64,1}, sampling::String, log::Array{String,1})RheologyData(colnames::Array{String,1}, data1::Array{Float64,1}, data2::Array{Float64,1}[, data3::Array{Float64,1}; filedir::String=\"none\", log::Array{String,1}=Array{String}(0)])"
},

{
    "location": "keytypes.html#Model-Structure-1",
    "page": "Rheology Types",
    "title": "Model Structure",
    "category": "section",
    "text": "RheologyModel(name::Function, parameters::Array{Float64,1}[, log::Array{String,1}])"
},

{
    "location": "keytypes.html#Generated-Data-Structure-1",
    "page": "Rheology Types",
    "title": "Generated Data Structure",
    "category": "section",
    "text": "RheologyArtificial(data::Array{Float64,1}, t::Array{Float64,1}, stepsize::Float64, log::Array{String,1})"
},

{
    "location": "preprocessing.html#RHEOS.var_resample-Tuple{RheologyData,Symbol,Float64}",
    "page": "Preprocessing",
    "title": "RHEOS.var_resample",
    "category": "method",
    "text": "var_resample(self::RheologyData, refvar::Symbol, pcntdownsample::Float64; mapback::Bool = false)\n\nConvert a fixed sample rate array to a variable sample rate, with sampling points added according to a relative change in chosen variable refvar, 1st derivative of refvar and 2nd derivative of refvar (WRT time). Usually chosen as the measured variable, so :σ for a stress relaxation test and :ϵ for a creep test.\n\nCurrently only variable downsampling supported. pcntdown sample is approximate, works well in some cases and very poorly in others. If required, compare resampled length vs original length after processing has finished. If data is noisy, may benefit from sending smoothed signal to this algorithm and either using mapback function or interpolating onto unsmoothed data.\n\nSee help docstring for var_resample for more details on algorithm implementation.\n\n\n\n\n\n"
},

{
    "location": "preprocessing.html#RHEOS.downsample-Tuple{RheologyData,Array{Float64,1},Array{Int64,1}}",
    "page": "Preprocessing",
    "title": "RHEOS.downsample",
    "category": "method",
    "text": "downsample(self::RheologyData, time_boundaries::Array{Float64,1}, elperiods::Array{Int64,1})\n\nHigh-level RheologyData interface to downsample in base.jl. Boundaries are floating point times which are then converted to the closest elements.\n\n\n\n\n\n"
},

{
    "location": "preprocessing.html#RHEOS.fixed_resample-Tuple{RheologyData,Array{Float64,1},Array{Int64,1},Array{String,1}}",
    "page": "Preprocessing",
    "title": "RHEOS.fixed_resample",
    "category": "method",
    "text": "fixed_resample(self::RheologyData, time_boundaries::Array{Float64,1}, elperiods::Array{Int64,1}, direction::Array{String,1})\n\nHigh-level RheologyData interface to fixed_resample in base.jl\n\n\n\n\n\n"
},

{
    "location": "preprocessing.html#RHEOS.smooth-Tuple{RheologyData,Float64}",
    "page": "Preprocessing",
    "title": "RHEOS.smooth",
    "category": "method",
    "text": "smooth(self::RheologyData, τ::Float64)\n\nSmooth data using a Gaussian Kernel to time scale τ (approximately half power).\n\nSmooths both σ and ϵ.\n\n\n\n\n\n"
},

{
    "location": "preprocessing.html#",
    "page": "Preprocessing",
    "title": "Preprocessing",
    "category": "page",
    "text": "var_resample(self::RheologyData, refvar::Symbol, pcntdownsample::Float64; mapback::Bool = false)downsample(self::RheologyData, time_boundaries::Array{Float64,1}, elperiods::Array{Int64,1})\nfixed_resample(self::RheologyData, time_boundaries::Array{Float64,1}, elperiods::Array{Int64,1}, direction::Array{String,1})smooth(self::RheologyData, τ::Float64)function mapbackdata(self_new::RheologyData, self_original::RheologyData)"
},

{
    "location": "fitpredict.html#",
    "page": "Fit and Predict",
    "title": "Fit and Predict",
    "category": "page",
    "text": "modelfit(data::RheologyData, modulus::Function[, p0::Array{Float64,1}, lo::Array{Float64,1}, hi::Array{Float64,1}; verbose::Bool = false])modelpredict(data::RheologyData, model::RheologyModel)"
},

{
    "location": "loadsave.html#",
    "page": "Save and Load",
    "title": "Save and Load",
    "category": "page",
    "text": "savedata(self::RheologyData; filedir::String = \"\", ext = \"_RheologyData.jld\")loaddata(filedir::String)savemodel(self::RheologyModel; filedir::String = \"\", ext = \"\")loadmodel(filedir::String)exportdata(self::RheologyData; filedir::String = \"\", ext = \"_mod.csv\")"
},

]}
