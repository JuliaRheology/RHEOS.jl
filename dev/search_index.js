var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#RHEOS.jl-1",
    "page": "Home",
    "title": "RHEOS.jl",
    "category": "section",
    "text": ""
},

{
    "location": "#Overview-1",
    "page": "Home",
    "title": "Overview",
    "category": "section",
    "text": "RHEOS, an abbreviation of Rheology Open Source, is a software package written in the Julia programming language that provides tools for analyzing rheological data. Features include:Stress/Strain/Time data can be easily be fitted to a viscoelastic model\nG\'/G\'\'/Frequency data can easily be fitted to a viscoelastic model\nMany standard and fractional viscoelastic models have already been implemented within RHEOS new ones can easily be added by users\nA fitted model can be used to predict the behaviour of the material under other loading conditions, enabling the fit/predict paradigm of model selection\nArtificial loading conditions can be generated within RHEOS to better understand a model\'s response"
},

{
    "location": "#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "Install Julia, version 1.0.2\nFrom Julia REPL, enter pkg mode by pressing ]\n(Optional) Enable desired Project.toml environment\nRun the command add \"https://github.com/JuliaRheology/RHEOS.jl\""
},

{
    "location": "#Included-Dependencies-1",
    "page": "Home",
    "title": "Included Dependencies",
    "category": "section",
    "text": ""
},

{
    "location": "#[FastConv.jl](https://github.com/aamini/FastConv.jl)-1",
    "page": "Home",
    "title": "FastConv.jl",
    "category": "section",
    "text": ""
},

{
    "location": "#[MittagLeffler.jl](https://github.com/jlapeyre/MittagLeffler.jl)-1",
    "page": "Home",
    "title": "MittagLeffler.jl",
    "category": "section",
    "text": ""
},

{
    "location": "#Citation-1",
    "page": "Home",
    "title": "Citation",
    "category": "section",
    "text": "If you use RHEOS in your work, please consider citing the following paper TBA"
},

{
    "location": "#References-1",
    "page": "Home",
    "title": "References",
    "category": "section",
    "text": "W. N. Findley, J. S. Lai, K. Onaran — Creep and Relaxation of Nonlinear Viscoelastic Materials (with an Introduction to Linear Viscoelasticity), Dover Publications, New York. (1976)\nS. G. Johnson — The NLopt nonlinear-optimization package, http://ab-initio.mit.edu/nlopt\nJ. Bezanson, A. Edelman, S. Karpinski, V. B. Shah — Julia: A Fresh Approach to Numerical Computing, SIAM Review, doi: 10.1137/141000671. (2017)"
},

{
    "location": "keytypes/#",
    "page": "Rheology Types",
    "title": "Rheology Types",
    "category": "page",
    "text": ""
},

{
    "location": "keytypes/#Rheology-Data-1",
    "page": "Rheology Types",
    "title": "Rheology Data",
    "category": "section",
    "text": ""
},

{
    "location": "keytypes/#Data-Structure-1",
    "page": "Rheology Types",
    "title": "Data Structure",
    "category": "section",
    "text": "RheologyData(σ::Array{Float64,1}, ϵ::Array{Float64,1}, t::Array{Float64,1}, sampling::String, log::Array{String,1})RheologyData(colnames::Array{String,1}, data1::Array{Float64,1}, data2::Array{Float64,1}[, data3::Array{Float64,1}; filedir::String=\"none\", log::Array{String,1}=Array{String}(0)])"
},

{
    "location": "keytypes/#Model-Structure-1",
    "page": "Rheology Types",
    "title": "Model Structure",
    "category": "section",
    "text": "RheologyModel(name::Function, parameters::Array{Float64,1}[, log::Array{String,1}])"
},

{
    "location": "keytypes/#Generated-Data-Structure-1",
    "page": "Rheology Types",
    "title": "Generated Data Structure",
    "category": "section",
    "text": "RheologyArtificial(data::Array{Float64,1}, t::Array{Float64,1}, stepsize::Float64, log::Array{String,1})"
},

{
    "location": "preprocessing/#RHEOS.var_resample-Tuple{RheologyData,Symbol,Float64}",
    "page": "Preprocessing",
    "title": "RHEOS.var_resample",
    "category": "method",
    "text": "var_resample(self::RheologyData, refvar::Symbol, pcntdownsample::Float64; mapback::Bool = false)\n\nConvert a fixed sample rate array to a variable sample rate, with sampling points added according to a relative change in chosen variable refvar, 1st derivative of refvar and 2nd derivative of refvar (WRT time). Usually chosen as the measured variable, so :σ for a stress relaxation test and :ϵ for a creep test.\n\nCurrently only variable downsampling supported. pcntdown sample is approximate, works well in some cases and very poorly in others. If required, compare resampled length vs original length after processing has finished. If data is noisy, may benefit from sending smoothed signal to this algorithm and either using mapback function or interpolating onto unsmoothed data.\n\nSee help docstring for var_resample for more details on algorithm implementation.\n\n\n\n\n\n"
},

{
    "location": "preprocessing/#RHEOS.downsample-Tuple{RheologyData,Array{Float64,1},Array{Int64,1}}",
    "page": "Preprocessing",
    "title": "RHEOS.downsample",
    "category": "method",
    "text": "downsample(self::RheologyData, time_boundaries::Array{Float64,1}, elperiods::Array{Int64,1})\n\nHigh-level RheologyData interface to downsample in base.jl. Boundaries are floating point times which are then converted to the closest elements.\n\n\n\n\n\n"
},

{
    "location": "preprocessing/#RHEOS.fixed_resample-Tuple{RheologyData,Array{Float64,1},Array{Int64,1},Array{String,1}}",
    "page": "Preprocessing",
    "title": "RHEOS.fixed_resample",
    "category": "method",
    "text": "fixed_resample(self::RheologyData, time_boundaries::Array{Float64,1}, elperiods::Array{Int64,1}, direction::Array{String,1})\n\nHigh-level RheologyData interface to fixed_resample in base.jl\n\n\n\n\n\n"
},

{
    "location": "preprocessing/#RHEOS.smooth-Tuple{RheologyData,Float64}",
    "page": "Preprocessing",
    "title": "RHEOS.smooth",
    "category": "method",
    "text": "smooth(self::RheologyData, τ::Float64)\n\nSmooth data using a Gaussian Kernel to time scale τ (approximately half power).\n\nSmooths both σ and ϵ.\n\n\n\n\n\n"
},

{
    "location": "preprocessing/#",
    "page": "Preprocessing",
    "title": "Preprocessing",
    "category": "page",
    "text": "var_resample(self::RheologyData, refvar::Symbol, pcntdownsample::Float64; mapback::Bool = false)downsample(self::RheologyData, time_boundaries::Array{Float64,1}, elperiods::Array{Int64,1})\nfixed_resample(self::RheologyData, time_boundaries::Array{Float64,1}, elperiods::Array{Int64,1}, direction::Array{String,1})smooth(self::RheologyData, τ::Float64)function mapbackdata(self_new::RheologyData, self_original::RheologyData)"
},

{
    "location": "fitpredict/#",
    "page": "Fit and Predict",
    "title": "Fit and Predict",
    "category": "page",
    "text": "modelfit(data::RheologyData, modulus::Function[, p0::Array{Float64,1}, lo::Array{Float64,1}, hi::Array{Float64,1}; verbose::Bool = false])modelpredict(data::RheologyData, model::RheologyModel)"
},

{
    "location": "loadsave/#RHEOS.loaddata-Tuple{String}",
    "page": "Save and Load",
    "title": "RHEOS.loaddata",
    "category": "method",
    "text": "loaddata(filedir::String)\n\nConvenience function loads RheologyData.\n\n\n\n\n\n"
},

{
    "location": "loadsave/#",
    "page": "Save and Load",
    "title": "Save and Load",
    "category": "page",
    "text": "savedata(self::RheologyData; filedir::String = \"\", ext = \"_RheologyData.jld\")loaddata(filedir::String)savemodel(self::RheologyModel; filedir::String = \"\", ext = \"\")loadmodel(filedir::String)exportdata(self::RheologyData; filedir::String = \"\", ext = \"_mod.csv\")"
},

]}
