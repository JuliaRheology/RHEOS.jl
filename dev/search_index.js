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
    "location": "fittingdata/#",
    "page": "Fitting Data",
    "title": "Fitting Data",
    "category": "page",
    "text": ""
},

{
    "location": "fittingdata/#Fitting-Data-1",
    "page": "Fitting Data",
    "title": "Fitting Data",
    "category": "section",
    "text": "This page contains examples which demonstrate all the model fitting capabilities of RHEOS.using RHEOS\n\nfiledir = \"DataComplete.csv\"\n\ndata = fileload([\"stress\",\"strain\", \"time\"], filedir)\n\n# SLS fit\nsls_fit = modelfit(data, SLS(), :G)\n\n# Spring-pot fit: cₐ, a, kᵦ, kᵧ\nlb = [0.1, 0.01, 0.1, 0.1]\nub = [Inf, 0.99, Inf, Inf]\nfractsls_fit = modelfit(data, FractionalSLS([2.0, 0.5, 0.5, 0.7]), :G; lo=lb, hi=ub, rel_tol=1e-5, verbose=true)\n\n# # get curves based on models fitted\nsls_predicted = modelpredict(data, sls_fit, :G)\nfractsls_predicted = modelpredict(data, fractsls_fit, :G)\n\n# plot all data\nfig, ax = subplots()\nax[:plot](data.t, data.σ, label=\"Data\", color=\"black\")\nax[:plot](sls_predicted.t, sls_predicted.σ, label=\"SLS\")\nax[:plot](fractsls_predicted.t, fractsls_predicted.σ, \"--\", label=\"Fractional SLS\")\nax[:legend](loc=\"best\")\nshow()"
},

{
    "location": "predictingresponse/#",
    "page": "Predicting Responses",
    "title": "Predicting Responses",
    "category": "page",
    "text": ""
},

{
    "location": "predictingresponse/#Predicting-Responses-1",
    "page": "Predicting Responses",
    "title": "Predicting Responses",
    "category": "section",
    "text": ""
},

{
    "location": "generatingdata/#",
    "page": "Generating Loading",
    "title": "Generating Loading",
    "category": "page",
    "text": ""
},

{
    "location": "generatingdata/#Generating-Data-1",
    "page": "Generating Loading",
    "title": "Generating Data",
    "category": "section",
    "text": ""
},

{
    "location": "preprocessing/#",
    "page": "Preprocessing Tools",
    "title": "Preprocessing Tools",
    "category": "page",
    "text": ""
},

{
    "location": "preprocessing/#Preprocessing-Tools-1",
    "page": "Preprocessing Tools",
    "title": "Preprocessing Tools",
    "category": "section",
    "text": ""
},

{
    "location": "fileIO/#",
    "page": "File I/O",
    "title": "File I/O",
    "category": "page",
    "text": ""
},

{
    "location": "fileIO/#File-I/O-1",
    "page": "File I/O",
    "title": "File I/O",
    "category": "section",
    "text": ""
},

{
    "location": "models/#",
    "page": "Models",
    "title": "Models",
    "category": "page",
    "text": ""
},

{
    "location": "models/#Models-1",
    "page": "Models",
    "title": "Models",
    "category": "section",
    "text": ""
},

{
    "location": "API/#",
    "page": "API",
    "title": "API",
    "category": "page",
    "text": ""
},

{
    "location": "API/#API-1",
    "page": "API",
    "title": "API",
    "category": "section",
    "text": ""
},

{
    "location": "API/#Data-Structure-1",
    "page": "API",
    "title": "Data Structure",
    "category": "section",
    "text": "RheologyData(σ::Array{Float64,1}, ϵ::Array{Float64,1}, t::Array{Float64,1}, sampling::String, log::Array{String,1})RheologyData(colnames::Array{String,1}, data1::Array{Float64,1}, data2::Array{Float64,1}[, data3::Array{Float64,1}; filedir::String=\"none\", log::Array{String,1}=Array{String}(0)])"
},

{
    "location": "API/#Model-Structure-1",
    "page": "API",
    "title": "Model Structure",
    "category": "section",
    "text": "RheologyModel(name::Function, parameters::Array{Float64,1}[, log::Array{String,1}])"
},

{
    "location": "API/#RHEOS.var_resample-Tuple{RheologyData,Symbol,Float64}",
    "page": "API",
    "title": "RHEOS.var_resample",
    "category": "method",
    "text": "var_resample(self::RheologyData, refvar::Symbol, pcntdownsample::Float64; mapback::Bool = false)\n\nConvert a fixed sample rate array to a variable sample rate, with sampling points added according to a relative change in chosen variable refvar, 1st derivative of refvar and 2nd derivative of refvar (WRT time). Usually chosen as the measured variable, so :σ for a stress relaxation test and :ϵ for a creep test.\n\nCurrently only variable downsampling supported. pcntdown sample is approximate, works well in some cases and very poorly in others. If required, compare resampled length vs original length after processing has finished. If data is noisy, may benefit from sending smoothed signal to this algorithm and either using mapback function or interpolating onto unsmoothed data.\n\nSee help docstring for var_resample for more details on algorithm implementation.\n\n\n\n\n\n"
},

{
    "location": "API/#RHEOS.downsample-Tuple{RheologyData,Array{Float64,1},Array{Int64,1}}",
    "page": "API",
    "title": "RHEOS.downsample",
    "category": "method",
    "text": "downsample(self::RheologyData, time_boundaries::Array{Float64,1}, elperiods::Array{Int64,1})\n\nHigh-level RheologyData interface to downsample in base.jl. Boundaries are floating point times which are then converted to the closest elements.\n\n\n\n\n\n"
},

{
    "location": "API/#RHEOS.fixed_resample-Tuple{RheologyData,Array{Float64,1},Array{Int64,1},Array{String,1}}",
    "page": "API",
    "title": "RHEOS.fixed_resample",
    "category": "method",
    "text": "fixed_resample(self::RheologyData, time_boundaries::Array{Float64,1}, elperiods::Array{Int64,1}, direction::Array{String,1})\n\nHigh-level RheologyData interface to fixed_resample in base.jl\n\n\n\n\n\n"
},

{
    "location": "API/#RHEOS.smooth-Tuple{RheologyData,Float64}",
    "page": "API",
    "title": "RHEOS.smooth",
    "category": "method",
    "text": "smooth(self::RheologyData, τ::Float64)\n\nSmooth data using a Gaussian Kernel to time scale τ (approximately half power).\n\nSmooths both σ and ϵ. Essentially a low pass filter with frequencies of 1/τ being cut to approximately half power. For other pad types available see ImageFiltering documentation.\n\n\n\n\n\n"
},

{
    "location": "API/#RHEOS.loaddata-Tuple{String}",
    "page": "API",
    "title": "RHEOS.loaddata",
    "category": "method",
    "text": "loaddata(filedir::String)\n\nConvenience function loads RheologyData.\n\n\n\n\n\n"
},

{
    "location": "API/#Generated-Data-Structure-1",
    "page": "API",
    "title": "Generated Data Structure",
    "category": "section",
    "text": "RheologyArtificial(data::Array{Float64,1}, t::Array{Float64,1}, stepsize::Float64, log::Array{String,1})modelfit(data::RheologyData, modulus::Function[, p0::Array{Float64,1}, lo::Array{Float64,1}, hi::Array{Float64,1}; verbose::Bool = false])modelpredict(data::RheologyData, model::RheologyModel)modelfit(data::RheologyData, modulus::Function[, p0::Array{Float64,1}, lo::Array{Float64,1}, hi::Array{Float64,1}; verbose::Bool = false])modelpredict(data::RheologyData, model::RheologyModel)var_resample(self::RheologyData, refvar::Symbol, pcntdownsample::Float64; mapback::Bool = false)downsample(self::RheologyData, time_boundaries::Array{Float64,1}, elperiods::Array{Int64,1})\nfixed_resample(self::RheologyData, time_boundaries::Array{Float64,1}, elperiods::Array{Int64,1}, direction::Array{String,1})smooth(self::RheologyData, τ::Float64)function mapbackdata(self_new::RheologyData, self_original::RheologyData)savedata(self::RheologyData; filedir::String = \"\", ext = \"_RheologyData.jld\")loaddata(filedir::String)savemodel(self::RheologyModel; filedir::String = \"\", ext = \"\")loadmodel(filedir::String)exportdata(self::RheologyData; filedir::String = \"\", ext = \"_mod.csv\")"
},

]}
