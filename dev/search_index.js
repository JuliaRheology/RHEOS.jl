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
    "text": "This page is a tutorial on how to fit viscoelastic models to data using RHEOS. If you want to try out the code below, it can all be run from the Julia REPL but note that the importing data functions will only work if you are using the \'RHEOS/examples\' folder as your working directory as that\'s where the example data files are stored."
},

{
    "location": "fittingdata/#Stress/Strain/Time-Data-1",
    "page": "Fitting Data",
    "title": "Stress/Strain/Time Data",
    "category": "section",
    "text": "This section is for standard viscoelastic tensile or compression tests where time, stress and strain data are available.First, we need to load in RHEOSusing RHEOSRHEOS has a convenience function for importing data from CSV files. The default column delimiter is \',\' but an alternative can be specified as a keyword argument. The row delimiter is a newline character (\'\\n\'). For standard viscoelastic testing data RHEOS expects either stress, strain and time data, just stress and time, or just strain and time. The order of the columns is specified as the first argument in the function importdata. The second argument is the directory of the file, as shown below.data = importdata([\"stress\",\"strain\", \"time\"], \"DataComplete.csv\")Now we have all our data stored in the variable data which is of type RheologyData. (In this tutorial, our data file DataComplete.csv is in the same directory as our Julia script so we can just use its relative directory.)Let\'s fit a Standard Linear Solid viscoelastic model via its relaxation modulus, G, as our data is from a stress relaxation test. The first argument is our data, the second argument tells RHEOS which model to fit and the final argument tells RHEOS whether to fit the model using a relaxation modulus (:G) or creep modulus (:J).fitted_SLS_model = modelfit(data, SLS(), :G)Our first fitted model is not contained in the fitted_SLS_model variable which is an instance of the RheologyModel data type.Next, we\'ll fit a fractional Standard Linear Solid model (the only difference from the above model is that the dash-pot is replaced by a spring-pot). This time we\'ll also add upper and lower bounds on the model parameters. This is highly recommended for fractional models in particular as values less than 0 or greater than 1 for the spring-pot parameter are unphysical and can cause errors in the Mittag-Leffler function used.lb = [0.1, 0.01, 0.1, 0.1]\nub = [Inf, 0.99, Inf, Inf]\nfitted_fractSLS_model = modelfit(data, FractionalSLS(), :G; lo=lb, hi=ub)Note the two keyword arguments used – lo and hi for the lower and upper parameter boundaries respectively. The special argument Inf for the three of the parameters\' upper bounds represent a type of infinity such that the parameters can be as large as required by the optimisation algorithm.For a full list of keyword arguments and features of the modelfit function, see the relevant part of the API section. Models included in RHEOS are also listed in the API section, and discussed in more detail in the Models section."
},

{
    "location": "fittingdata/#G*/ω-Data-1",
    "page": "Fitting Data",
    "title": "G*/ω Data",
    "category": "section",
    "text": "RHEOS can also fit models to dynamic mechanical analysis data from oscillatory tests. The importdata function can again be used here but with the column names Gp for the storage modulus, Gpp for the loss modulus and frequency for the frequency column. RHEOS will detect the frequency string and know to load the data into a RheologyDynamic data type. Assuming RHEOS has already been imported, let\'s load in our data file:data = importdata([\"Gp\",\"Gpp\", \"Frequency\"], \"FrequencyData.csv\")As this is dynamic mechanical testing data it is an instance of the RheologyDynamic data type. As before, we\'ll try and fit the data to a Standard Linear Solid model – but this time we need to use the dynamicmodelfit function.fitted_SLS_model = dynamicmodelfit(data, SLS())Note that we do not have to specify whether we are fitting via the creep or relaxation modulus as RHEOS always fits dynamic data to the complex (dynamic) modulus. Let\'s also try to fit the data to a fractional Standard Linear Solid model:lb = [0.1, 0.01, 0.1, 0.1]\nub = [Inf, 0.99, Inf, Inf]\nfitted_fractSLS_model = dynamicmodelfit(data, FractionalSLS(); lo=lb, hi=ub)For more information on the dynamicmodelfit function, see the API section."
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
    "location": "moreexamples/#",
    "page": "More Examples",
    "title": "More Examples",
    "category": "page",
    "text": ""
},

{
    "location": "moreexamples/#More-Examples-1",
    "page": "More Examples",
    "title": "More Examples",
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
    "location": "API/#RHEOS.RheologyData",
    "page": "API",
    "title": "RHEOS.RheologyData",
    "category": "type",
    "text": "RheologyData(σ::Vector{T}, ϵ::Vector{T}, t::Vector{T}, sampling::String, log::Vector{String}) where T<:Real\n\nRheologyData struct contains stress, strain and time data.\n\nIf preferred, an instance can be generated manually by just providing the three data vectors in the right order, sampling type will be checked automatically. If loading partial data (either stress or strain), fill the other vector as a vector of zeros of the same length as the others.\n\nFields\n\nσ: stress\nϵ: strain\nt: time\nsampling: sampling type, either \"constant\" or \"variable\"\nlog: a log of struct\'s events, e.g. preprocessing\n\n\n\n\n\n"
},

{
    "location": "API/#RHEOS.RheologyDynamic",
    "page": "API",
    "title": "RHEOS.RheologyDynamic",
    "category": "type",
    "text": "RheologyDynamic(Gp::Vector{T}, Gpp::Vector{T}, ω::Vector{T}, log::Vector{String}) where T<:Real\n\nRheologyDynamic contains storage modulus, loss modulus and frequency data.\n\nIf preferred, an instance can be generated manually by just providing the three data vectors in the right order.\n\nFields\n\nGp: storage modulus\nGpp: loss modulus\nω: frequency\nlog: a log of struct\'s events, e.g. preprocessing\n\n\n\n\n\n"
},

{
    "location": "API/#RHEOS.RheologyModel",
    "page": "API",
    "title": "RHEOS.RheologyModel",
    "category": "type",
    "text": "RheologyModel(G::T, J::T, Gp::T, Gpp::T, parameters::Vector{S<:Real} log::Vector{String}) where T<:Function\n\nRheologyModel contains a model\'s various moduli, parameters, and log of activity.\n\nFor incomplete models, an alternative constructor is available where all arguments are keyword arguments and moduli not provided default to a null modulus which always returns [-1.0].\n\nFields\n\nG: Relaxation modulus\nJ: Creep modulus\nGp: Storage modulus\nGpp: Loss modulus\nparameters: Used for predicting and as default starting parameters in fitting\nlog: a log of struct\'s events, e.g. what file it was fitted to\n\n\n\n\n\n"
},

{
    "location": "API/#RHEOS-Types-1",
    "page": "API",
    "title": "RHEOS Types",
    "category": "section",
    "text": "RheologyData\nRheologyDynamic\nRheologyModel"
},

{
    "location": "API/#Preprocessing-Functions-1",
    "page": "API",
    "title": "Preprocessing Functions",
    "category": "section",
    "text": ""
},

]}
