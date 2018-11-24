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
    "location": "#Documentation-1",
    "page": "Home",
    "title": "Documentation",
    "category": "section",
    "text": "The sections in this documentation each aim to provide tutorials in different elements of RHEOS. The API section is a comprehensive list of RHEOS types and functions brief descriptions of their use. For corrections or further questions, please create an issue on the Github repository."
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
    "text": "This page is a tutorial on how to fit viscoelastic models to data using RHEOS. If you want to try out the code below, it can all be run from the Julia REPL but note that importing the data exactly as shown will only work if you are using the \'RHEOS/examples\' folder as your working directory as that\'s where the example data files are stored and the directories used are relative."
},

{
    "location": "fittingdata/#Stress/Strain/Time-Data-1",
    "page": "Fitting Data",
    "title": "Stress/Strain/Time Data",
    "category": "section",
    "text": "This section is for standard viscoelastic tensile or compression tests where time, stress and strain data are available.First, we need to load in RHEOSusing RHEOSRHEOS has a convenience function for importing data from CSV files. The default column delimiter is \',\' but an alternative can be specified as a keyword argument. The row delimiter is a newline character (\'\\n\'). For standard viscoelastic testing data RHEOS expects either stress, strain and time data, just stress and time, or just strain and time. The order of the columns is specified as the first argument in the function importdata. The second argument is the directory of the file, as shown below.data = importdata([\"stress\",\"strain\", \"time\"], \"DataRelaxation.csv\")Now we have all our data stored in the variable data which is of type RheologyData. (In this tutorial, our data file DataRelaxation.csv is in the same directory as our Julia script so we can just use its relative directory.)Let\'s fit a Standard Linear Solid viscoelastic model via its relaxation modulus, G, as our data is from a stress relaxation test. The first argument is our data, the second argument tells RHEOS which model to fit and the final argument tells RHEOS whether to fit the model using a relaxation modulus (:G) or creep modulus (:J).fitted_SLS_model = modelfit(data, SLS(), :G)Our first fitted model is not contained in the fitted_SLS_model variable which is an instance of the RheologyModel data type.Next, we\'ll fit a fractional Standard Linear Solid model (the only difference from the above model is that the dash-pot is replaced by a spring-pot). This time we\'ll also add upper and lower bounds on the model parameters. This is highly recommended for fractional models in particular as values less than 0 or greater than 1 for the spring-pot parameter are unphysical and can cause errors in the Mittag-Leffler function used.lb = [0.1, 0.01, 0.1, 0.1]\nub = [Inf, 0.99, Inf, Inf]\nfitted_fractSLS_model = modelfit(data, FractionalSLS(), :G; lo=lb, hi=ub)Note the two keyword arguments used – lo and hi for the lower and upper parameter boundaries respectively. The special argument Inf for the three of the parameters\' upper bounds represent a type of infinity such that the parameters can be as large as required by the optimisation algorithm.For a full list of keyword arguments and features of the modelfit function, see the relevant part of the API section. Models included in RHEOS are also listed in the API section, and discussed in more detail in the Models section."
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
    "location": "predictingresponse/#Stress/Strain/Time-Data-1",
    "page": "Predicting Responses",
    "title": "Stress/Strain/Time Data",
    "category": "section",
    "text": "Now that we know how to import datasets and fit models to them, let\'s try and make some predictions based on fitted model parameters. We will use our data imported and models fitted_SLS_model and fitted_fractSLS_model as fitted in the previous section. First we\'ll use model prediction to see how well the models were fitted to the relaxation data.SLS_predicted = modelpredict(data, fitted_SLS_model, :G)\n\nfractSLS_predicted = modelpredict(data, fitted_fractSLS_model, :G)Now we have the predicted data, we can plot. Any plotting library can be used but the example code below uses the PyPlot Julia package.using PyPlot\n\nfig, ax = subplots()\nax[:plot](data.t, data.σ, label=\"Data\", color=\"black\")\nax[:plot](SLS_predicted.t, SLS_predicted.σ, label=\"SLS\")\nax[:plot](fractSLS_predicted.t, fractSLS_predicted.σ, \"--\", label=\"Fractional SLS\")\nax[:legend](loc=\"best\")\nshow()(Image: relaxation prediction) We can see from the above plot that the fractional Standard Linear Solid model has yielded a much better fit than that regular Standard Linear Solid. Now we also have creep data generated using the same model. Let\'s see how well our fitted models perform making predictions on a completely different type of data.creep_data = importdata([\"stress\", \"strain\", \"time\"], \"DataCreep.csv\")\n\ncreep_SLS_predicted = modelpredict(creep_data, fitted_SLS_model, :J)\n\ncreep_fractSLS_predicted = modelpredict(creep_data, fitted_fractSLS_model, :J)Note the two differences in the modelpredict function calls. In the first argument we are predicting based on the creep data so we use the creep data we imported in the first line of the above code block. In the third argument we use the argument :J instead of :G to make predictions based on the creep modulus, not the relaxation modulus as before. Plotting the above but strain this time not stress as it\'s using the creep modulus:fig, ax = subplots()\nax[:plot](creep_data.t, creep_data.ϵ, label=\"Creep Data\", color=\"black\")\nax[:plot](creep_SLS_predicted.t, creep_SLS_predicted.ϵ, label=\"SLS\")\nax[:plot](creep_fractSLS_predicted.t, creep_fractSLS_predicted.ϵ, \"--\", label=\"Fractional SLS\")\nax[:legend](loc=\"best\")\nshow()(Image: creep prediction)Again we see that the fractional SLS model has yielded a far better prediction, even on data not seen during the fitting stage."
},

{
    "location": "predictingresponse/#Frequency-Data-1",
    "page": "Predicting Responses",
    "title": "Frequency Data",
    "category": "section",
    "text": "Although it is very similar to the above, for completeness this section demonstrates how to predict frequency-based rheology data.frequency_data = importdata([\"Gp\", \"Gpp\", \"frequency\"], \"FrequencyData.csv\")\n\nfreq_fractSLS_predicted = dynamicmodelpredict(frequency_data, fitted_fractSLS_model)See the dynamicmodelpredict for more information on the function itself."
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
    "text": "RHEOS has several convenience functions for generating arbitrarily complex loading patterns. These may be particularly useful for investigating the responses of viscoelastic models with which users are unfamiliar. This section will demonstrate how to use them. It should be noted at the outset that the way these functions are currently implemented, they generate the same loading in both stress and strain with the expectation that users will then use one or other to make predictions. As all the structs generated are of RheologyData type, the same addition, subtraction and multiplication overloaded methods can be used for real data. When adding two RheologyData structs and one is longer than the other (in time), the shorter one will be extended by keeping the last values of that shorter struct\'s data constant for the rest of time. Adding, subtracting and multiplying will raise an error if the data do not have the same sample rate. All plots here are generated using the PyPlot Julia package."
},

{
    "location": "generatingdata/#Step,-Ramp-and-Oscillatory-Loading-1",
    "page": "Generating Loading",
    "title": "Step, Ramp and Oscillatory Loading",
    "category": "section",
    "text": "The code below uses stepgen to create one step starting at 100 seconds (with total duration of 1000 seconds) and another step starting at 500 seconds lasting the same total duration. The first argument of stepgen determines the total length in seconds. The second step is then subtracted from the first to create a new combined loading pattern as shown in the plots below. foo = stepgen(1000, 100)\n\nbar = stepgen(1000, 500)\n\nbaz = foo - bar(Image: step gif)The above example uses an \'instantaneous\' step. However, a logistic functional transition can be used by adding a non-zero t_trans keyword argument. Next we\'ll generate ramp loading using rampgen. The arguments in order are the total time length (as before), the time to start the ramp and the time to stop the ramp.foo = rampgen(1000, 100, 200)\n\nbar = rampgen(1000, 500, 700)\n\nbaz = 2*foo - bar(Image: ramp)And finally some oscillatory loading, ramp loading, and oscillatory loading multiplied by the ramp loading. The first argument of singen is the total time length, the second is the frequency in hertz, and the keyword argument is phase in radians.foo = singen(1000, 1/50; phase = -π/2)\n\nbar = rampgen(1000, 10, 400) - rampgen(1000, 400, 800)\n\nbaz = foo*bar(Image: osci)"
},

{
    "location": "generatingdata/#Repeated-Loading,-Adding-Noise-1",
    "page": "Generating Loading",
    "title": "Repeated Loading, Adding Noise",
    "category": "section",
    "text": "For repeated loading, RHEOS provides a convenience function that can loop loading patterns a specified number of times. Similar to the step function, repeatdata also offers a t_trans keyword argument which determines the transition between the end of one loop and the start of the next. If t_trans=0 then the transition is instantaneous, if not the transition occurs by logistic function with approximate transition time t_trans. The repeatdata function currently only works if the stress and strain arrays contain the same data (as is produced by all the data generation functions). Below we use our most recently defined baz variable (the oscillatory loading multiplied by the two ramps) to demonstrate.repeatedbaz = repeatdata(baz, 3)(Image: repeatedbaz)Finally, we can add uniform white noise to the data using the noisegen function. The below example demonstrates this on a simple step.foo = stepgen(100, 50)\nbar = 0.01*noisegen(100)\nbaz = foo + bar(Image: noise)"
},

{
    "location": "generatingdata/#A-More-Complicated-Example-1",
    "page": "Generating Loading",
    "title": "A More Complicated Example",
    "category": "section",
    "text": "Below is an example which uses almost all of the RHEOS data generation functionality together in one example.## Combined Example\n\n# generate a single step at 25 seconds\nstepup = stepgen(50.0, 25.0; stepsize = 0.05, t_trans = 2.5)\n\n# generate an oscillation which starts fading in at 25.5 seconds and has faded out by 49.5 seconds\nosci = 0.1*singen(50.0, 0.2; stepsize = 0.05)\nrampup = rampgen(50.0, 25.5, 37.5; stepsize = 0.05)\nrampdown = -rampgen(50.0, 37.5, 49.5; stepsize = 0.05)\n\n# combine the step and faded oscillation\ncombined = osci*(rampup + rampdown) + stepup\n\n# repeat this three times\nrepeated = repeatdata(combined, 3)\n\n# add some white noise with amplitude of 0.01\nnoisyrepeated = repeated + 0.01*noisegen(150.0; seed = 1, stepsize = 0.05)which, when plotted, produces the following (Image: complicated)we cam zoom in to a stepped region to see the noise more clearly (Image: complicatedzoom)"
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
    "text": "RHEOS offers several preprocessing functions which may come in useful. This page is intended to be a brief tutorial of their use. For detailed descriptions of functions and their optional arguments, see the API section."
},

{
    "location": "preprocessing/#Downsampling-1",
    "page": "Preprocessing Tools",
    "title": "Downsampling",
    "category": "section",
    "text": ""
},

{
    "location": "preprocessing/#Upsampling-and-Downsampling-1",
    "page": "Preprocessing Tools",
    "title": "Upsampling and Downsampling",
    "category": "section",
    "text": ""
},

{
    "location": "preprocessing/#Variable-Resampling-1",
    "page": "Preprocessing Tools",
    "title": "Variable Resampling",
    "category": "section",
    "text": ""
},

{
    "location": "preprocessing/#Smoothing-1",
    "page": "Preprocessing Tools",
    "title": "Smoothing",
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
    "location": "fileIO/#CSV-Import/Export-1",
    "page": "File I/O",
    "title": "CSV Import/Export",
    "category": "section",
    "text": "In the Fitting Data section we saw how to import data from csv files. As a brief overview of that functionality, the 3 lines below demonstrate how to import stress/strain/time data, partial strain/time data, and G\'/G\'\'/frequency data. The columns can be in any order but need to specified by the first argument array of strings shown below. Note that if you want to try and run this code exactly as shown below, the data files used in these examples are stored in the \'RHEOS/examples\' directory.data = importdata([\"stress\",\"strain\", \"time\"], \"DataRelaxation.csv\")\n\ndata_incomplete = importdata([\"strain\", \"time\"], \"DataIncomplete.csv\")\n\ndata_dynamic = importdata([\"Gp\",\"Gpp\", \"Frequency\"], \"FrequencyData.csv\")Where the first two function calls return RheologyData objects and the last function call returns a RheologyDynamic object.If you want to analyse or plot your data in software other than Julia you will likely want to export it to a CSV file. To export RheologyData and RheologyDynamic objects to CSV files we can use the exportdata function. For the two complete data-sets we imported above, we can export them into new files in the following way.exportdata(data, \"exported_data\")\n\nexportdata(data_dynamic, \"exported_dynamic_data\")The \'.csv\' extension will automatically be added (and can be modified if needs be through use of a keyword argument). The delimiter can also be changed if necessary. See API section for information on this."
},

{
    "location": "fileIO/#Native-Julia-data-files-(JLD2)-1",
    "page": "File I/O",
    "title": "Native Julia data files (JLD2)",
    "category": "section",
    "text": "When using RHEOS and doing all subsequent analysis and plotting in Julia, it is convenient to be able to store native Julia objects to disk so that they can be loaded in to subsequent Julia sessions. RHEOS provides two pairs of convenience functions to facilitate this: savedata and loaddata, savemodel and loadmodel. Using our data imported above we can demonstrate saving and loading data (it is exactly the same process for both RHEOS data types).savedata(data, \"imported_data\")\n\nloaded_data = loaddata(\"imported_data.jld2\")The savedata function will automatically append a \'.jld2\' to the filename when saving but this can be modified if necessary through a keyword argument. savedata and loaddata work on both RheologyData and  RheologyDynamic objects.Lastly, let\'s create a Standard Linear Solid model object with default parameters, save that to disk and load it again.model = SLS()\n\nsavemodel(model, \"saved_model\")\n\nloaded_model = loadmodel(\"saved_model.jld2\")"
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
    "location": "models/#Models-included-in-RHEOS-1",
    "page": "Models",
    "title": "Models included in RHEOS",
    "category": "section",
    "text": "Severa models (i.e. their relaxation modulus, creep modulus and complex modulus) already implemented in RHEOS. Their constructors are listed below. Note that the params argument is always optional, as represented by their containment in square barckets. If left blank then the model\'s default parameters are used. Additional models are forthcoming."
},

{
    "location": "models/#RHEOS.Spring",
    "page": "Models",
    "title": "RHEOS.Spring",
    "category": "function",
    "text": "Spring([params::Vector{T}]) where T<:Real\n\nSingle spring element.\n\n\n\n\n\n"
},

{
    "location": "models/#RHEOS.DashPot",
    "page": "Models",
    "title": "RHEOS.DashPot",
    "category": "function",
    "text": "DashPot([params::Vector{T}]) where T<:Real\n\nSingle dashpot element.\n\n\n\n\n\n"
},

{
    "location": "models/#RHEOS.SpringPot",
    "page": "Models",
    "title": "RHEOS.SpringPot",
    "category": "function",
    "text": "SpringPot([params::Vector{T}]) where T<:Real\n\nSingle springpot element.\n\n\n\n\n\n"
},

{
    "location": "models/#Elements-1",
    "page": "Models",
    "title": "Elements",
    "category": "section",
    "text": "Spring\nDashPot\nSpringPot"
},

{
    "location": "models/#RHEOS.FractionalMaxwell",
    "page": "Models",
    "title": "RHEOS.FractionalMaxwell",
    "category": "function",
    "text": "FractionalMaxwell([params::Vector{T}]) where T<:Real\n\nTwo springpots in series.\n\n\n\n\n\n"
},

{
    "location": "models/#RHEOS.FractionalMaxwellSpring",
    "page": "Models",
    "title": "RHEOS.FractionalMaxwellSpring",
    "category": "function",
    "text": "FractionalMaxwellSpring([params::Vector{T}]) where T<:Real\n\nA springpot and spring in series.\n\n\n\n\n\n"
},

{
    "location": "models/#RHEOS.FractionalMaxwellDashpot",
    "page": "Models",
    "title": "RHEOS.FractionalMaxwellDashpot",
    "category": "function",
    "text": "FractionalMaxwellDashpot([params::Vector{T}]) where T<:Real\n\nA springpot and dashpot in series.\n\n\n\n\n\n"
},

{
    "location": "models/#RHEOS.Maxwell",
    "page": "Models",
    "title": "RHEOS.Maxwell",
    "category": "function",
    "text": "Maxwell([params::Vector{T}]) where T<:Real\n\nA spring and dashpot in series.\n\n\n\n\n\n"
},

{
    "location": "models/#Maxwell-Type-1",
    "page": "Models",
    "title": "Maxwell Type",
    "category": "section",
    "text": "FractionalMaxwell\nFractionalMaxwellSpring\nFractionalMaxwellDashpot\nMaxwell"
},

{
    "location": "models/#RHEOS.FractionalKelvinVoigt",
    "page": "Models",
    "title": "RHEOS.FractionalKelvinVoigt",
    "category": "function",
    "text": "FractionalKelvinVoigt([params::Vector{T}]) where T<:Real\n\nTwo springpots in parallel.\n\n\n\n\n\n"
},

{
    "location": "models/#RHEOS.FractionalKVspring",
    "page": "Models",
    "title": "RHEOS.FractionalKVspring",
    "category": "function",
    "text": "FractionalKVspring([params::Vector{T}]) where T<:Real\n\nA springpot and spring in parallel.\n\n\n\n\n\n"
},

{
    "location": "models/#RHEOS.FractionalKVdashpot",
    "page": "Models",
    "title": "RHEOS.FractionalKVdashpot",
    "category": "function",
    "text": "FractionalKVdashpot([params::Vector{T}]) where T<:Real\n\nA springpot and dashpot in parallel.\n\n\n\n\n\n"
},

{
    "location": "models/#RHEOS.KelvinVoigt",
    "page": "Models",
    "title": "RHEOS.KelvinVoigt",
    "category": "function",
    "text": "KelvinVoigt([params::Vector{T}]) where T<:Real\n\nA spring and dashpot in parallel.\n\n\n\n\n\n"
},

{
    "location": "models/#Kelvin-Voigt-Type-1",
    "page": "Models",
    "title": "Kelvin-Voigt Type",
    "category": "section",
    "text": "FractionalKelvinVoigt\nFractionalKVspring\nFractionalKVdashpot\nKelvinVoigt"
},

{
    "location": "models/#RHEOS.FractionalZener",
    "page": "Models",
    "title": "RHEOS.FractionalZener",
    "category": "function",
    "text": "FractionalZener([params::Vector{T}]) where T<:Real\n\n2 springpots in series, both in parallel with a springpot.\n\n\n\n\n\n"
},

{
    "location": "models/#RHEOS.FractionalSLS",
    "page": "Models",
    "title": "RHEOS.FractionalSLS",
    "category": "function",
    "text": "FractionalSLS([params::Vector{T}]) where T<:Real\n\nA springpot and spring in series, both in parallel with a spring.\n\n\n\n\n\n"
},

{
    "location": "models/#RHEOS.FractionalSpecial",
    "page": "Models",
    "title": "RHEOS.FractionalSpecial",
    "category": "function",
    "text": "FractionalSpecial([params::Vector{T}]) where T<:Real\n\nA springpot and dashpot in series, both in parallel with a spring.\n\n\n\n\n\n"
},

{
    "location": "models/#RHEOS.SLS",
    "page": "Models",
    "title": "RHEOS.SLS",
    "category": "function",
    "text": "SLS(params::Vector{T}) where T<:Real\n\nA spring and dashpot in series, both in parallel with a spring.\n\n\n\n\n\n"
},

{
    "location": "models/#Zener-Type-1",
    "page": "Models",
    "title": "Zener Type",
    "category": "section",
    "text": "FractionalZener\nFractionalSLS\nFractionalSpecial\nSLS"
},

{
    "location": "models/#Creating-your-own-model-1",
    "page": "Models",
    "title": "Creating your own model",
    "category": "section",
    "text": "If you know some (or all) of the moduli for a model that you would like use but hasn\'t already been implemented in RHEOS, this section will explain how to quickly import these moduli into a RheologyModel object for use with other parts of RHEOS. For the sake of example, we will use the relaxation modulus, storage modulus and loss modulus of the Standard Linear Solid model as defined in RHEOS.function G_sls(t::Vector{T}, params::Vector{T}) where T<:Real\n    η, kᵦ, kᵧ = params\n\n    G = kᵧ .+ kᵦ*exp.(-t*kᵦ/η)\nend\n\nfunction Gp_sls(ω::Vector{T}, params::Vector{T}) where T<:Real\n    η, kᵦ, kᵧ = params\n\n    τ = η/kᵦ\n\n    denominator = 1 .+ (τ^2)*(ω.^2)\n    numerator = (ω.^2)*(τ^2)*kᵦ\n\n    Gp = numerator./denominator .+ kᵧ\nend\n\nfunction Gpp_sls(ω::{T}, params::Vector{T}) where T<:Real\n    η, kᵦ, kᵧ = params\n\n    τ = η/kᵦ\n\n    denominator = 1 .+ (τ^2)*(ω.^2)\n    numerator = ω*τ*kᵦ\n\n    Gpp = numerator./denominator\nendNow we have the our moduli defined as Julia functions we can store them, along with some (optional) default parameters, in a RheologyModel struct in the following way.our_model = RheologyModel(G = G_sls, Gp = Gp_sls, Gpp = Gpp_sls, params = [1.0, 0.5, 1.0])Now we can fit this model to data and use it to make predictions. Any moduli not included in this final step will default to a null_modulus which always returns the array [-1.0]."
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
    "text": "This page is a list of all types and functions defined by RHEOS that are intended for general use."
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
    "location": "API/#RHEOS.downsample",
    "page": "API",
    "title": "RHEOS.downsample",
    "category": "function",
    "text": "downsample(self::RheologyData, time_boundaries::Vector{T} where T<:Real, elperiods::Vector{S} where S<:Integer)\n\nBoundaries are floating point times which are then converted to the closest elements. Works by just reducing on indices. For example, time_boundaries could be [0.0, 10.0, 100.0] and elperiods could be [2, 4]. So new data would take every 2nd element from 0.0 seconds to 10.0 seconds, then every 4th element from 10.0 seconds to 100.0 seconds.\n\n\n\n\n\n"
},

{
    "location": "API/#RHEOS.fixedresample",
    "page": "API",
    "title": "RHEOS.fixedresample",
    "category": "function",
    "text": "fixedresample(self::RheologyData, time_boundaries::Vector{T} where T<:Real, elperiods::Vector{K} where K<:Integer, direction::Vector{String})\n\nResample data with new sample rate(s).\n\nWhereas downsample can only reduce the sample rate by not taking every array element, fixedresample can also upsample. Whether to up or down sample for a given section is known from the direction argument.\n\n\n\n\n\n"
},

{
    "location": "API/#RHEOS.variableresample",
    "page": "API",
    "title": "RHEOS.variableresample",
    "category": "function",
    "text": "variableresample(self::RheologyData, refvar::Symbol, pcntdownsample::Real; mapback::Bool = false)\n\nConvert a fixed sample rate array to a variable sample rate, with sampling points added according to a relative change in chosen variable refvar, 1st derivative of refvar and 2nd derivative of refvar (WRT time). Usually chosen as the measured variable, so :σ for a stress relaxation test and :ϵ for a creep test.\n\nCurrently only variable downsampling supported. pcntdown sample is approximate, works well in some cases and very poorly in others. If required, compare resampled length vs original length after processing has finished. If data is noisy, may benefit from sending smoothed signal to this algorithm and either using mapback function or interpolating onto unsmoothed data.\n\nAlgorithm works as follows. 25 initial samples are generated evenly spread. After this,  the array is repeatedly sweeped, anywhere Δy, Δdy/dx, Δd2y/dx2 is greater than a threshold, α, a new sample is created at the midpoint of the two tested points. This is allowed to happen a maximum of 400 times, after which α is decreased and the process starts again. This macro process continues until the desired pcntdownsample ratio has been reached.\n\nArguments\n\nself: RheologyData instance\nrefvar: The data whose derivatives will determine sample densities\npcntdownsample: Approximate ratio of new samples to old samples\n_mapback = false: (Optional) Determines whether resampled points should \'snap\' to closest original points\n\n\n\n\n\n"
},

{
    "location": "API/#RHEOS.smooth",
    "page": "API",
    "title": "RHEOS.smooth",
    "category": "function",
    "text": "smooth(self::RheologyData, τ::Real; pad::String=\"replicate\")\n\nSmooth data using a Gaussian Kernel to time scale τ (approximately half power).\n\nSmooths both σ and ϵ. Essentially a low pass filter with frequencies of 1/τ being cut to approximately half power. For other pad types available see ImageFiltering documentation.\n\n\n\n\n\n"
},

{
    "location": "API/#RHEOS.zerotime",
    "page": "API",
    "title": "RHEOS.zerotime",
    "category": "function",
    "text": "function zerotime(self::RheologyData)\n\nConvenience function to normalize time such that the starting time is 0.0\n\n\n\n\n\n"
},

{
    "location": "API/#Preprocessing-Functions-1",
    "page": "API",
    "title": "Preprocessing Functions",
    "category": "section",
    "text": "downsample\nfixedresample\nvariableresample\nsmooth\nzerotime"
},

{
    "location": "API/#RHEOS.modelfit",
    "page": "API",
    "title": "RHEOS.modelfit",
    "category": "function",
    "text": "modelfit(data::RheologyData, model::RheologyModel, modtouse::Symbol; p0::Vector{T} = [-1.0], lo::Vector{T} = [-1.0], hi::Vector{T} = [-1.0], verbose::Bool = false, rel_tol = 1e-4, diff_method=\"BD\") where T<:Real\n\nFit RheologyData struct to model and return a fitted model as a RheologyModel object.\n\nArguments\n\ndata: RheologyData struct containing all data\nmodel: RheologyModel containing moduli and default (initial) parameters\nmodtouse: :G for relaxation modulus, :J for creep modulus\np0: Initial parameters to use in fit (uses \'model\' parameters if none given)\nlo: Lower bounds for parameters\nhi: Upper bounds for parameters\nverbose: If true, prints parameters on each optimisation iteration\nrel_tol: Relative tolerance of optimization, see NLOpt docs for more details\ndiff_method: Set finite difference formula to use for derivative, currently \"BD\" or \"CD\"\n\n\n\n\n\n"
},

{
    "location": "API/#RHEOS.modelpredict",
    "page": "API",
    "title": "RHEOS.modelpredict",
    "category": "function",
    "text": "modelpredict(data::RheologyData, model::RheologyModel, modtouse::Symbol; diff_method=\"BD\")\n\nGiven data and model, return new dataset based on model parameters and using the modulus specified by \'modtouse\'; either creep modulus (:J, only returned strain is new) or  relaxation modulus (:G, only returned stress is new). \'diff_method\' sets finite difference for  calculating the derivative used in the hereditary integral and can be either backwards difference (\"BD\") or central difference (\"CD\").\n\n\n\n\n\n"
},

{
    "location": "API/#RHEOS.modelstepfit",
    "page": "API",
    "title": "RHEOS.modelstepfit",
    "category": "function",
    "text": "modelstepfit(data::RheologyData, model::RheologyModel, modtouse::Symbol; p0::Vector{T} = [-1.0], lo::Vector{T} = [-1.0], hi::Vector{T} = [-1.0], verbose::Bool = false, rel_tol = 1e-4) where T<:Real\n\nSame as \'modelfit\' except assumes a step loading. If this assumption is appropriate for the data then fitting can be sped up greatly by use of this function. If modtouse is :G, relaxation modulus, then the first element of the strain is assumed to be the amplitude of the step. If modtouse is :j, creep modulus, then the first element of the stress is assumed to be the amplitude of the step.\n\nArguments\n\ndata: RheologyData struct containing all data\nmodel: RheologyModel containing moduli and default (initial) parameters\nmodtouse: :G for relaxation modulus, :J for creep modulus\np0: Initial parameters to use in fit (uses \'model\' parameters if none given)\nlo: Lower bounds for parameters\nhi: Upper bounds for parameters\nverbose: If true, prints parameters on each optimisation iteration\nrel_tol: Relative tolerance of optimization, see NLOpt docs for more details\n\n\n\n\n\n"
},

{
    "location": "API/#RHEOS.modelsteppredict",
    "page": "API",
    "title": "RHEOS.modelsteppredict",
    "category": "function",
    "text": "modelsteppredict(data::RheologyData, model::RheologyModel, modtouse::Symbol; step_on::Real = 0.0)\n\nSame as modelpredict but assumes a step loading with step starting at \'step_on\'. Singularities are bypassed by adding 1 to the index of the singular element. \n\n\n\n\n\n"
},

{
    "location": "API/#RHEOS.dynamicmodelfit",
    "page": "API",
    "title": "RHEOS.dynamicmodelfit",
    "category": "function",
    "text": "dynamicmodelfit(data::RheologyDynamic, model::RheologyModel; p0::Vector{T} = [-1.0], lo::Vector{T} = [-1.0], hi::Vector{T} = [-1.0], verbose::Bool = false, rel_tol = 1e-4) where T<:Real\n\nFits model to the frequency/loss+storage moduli data.\n\nAll arguments are as described below. The \'weights\' argument some more information. As this fitting procedure is fitting two functions simultaneously (the storage and loss moduli), if left untransformed the fit would tend to favour the  modulus which is larger in magnitude and not fit the other modulus well. To avoid this, RHEOS offers a number of data transforms which can be used. \n\nArguments\n\ndata: RheologyDynamic struct containing all data\nmodel: RheologyModel containing moduli and default (initial) parameters\np0: Initial parameters to use in fit (uses \'model\' parameters if none given)\nlo: Lower bounds for parameters\nhi: Upper bounds for parameters\nverbose: If true, prints parameters on each optimisation iteration\nrel_tol: Relative tolerance of optimization, see NLOpt docs for more details\nweights: Weighting mode for storage and loss modulus (see above)\n\n\n\n\n\n"
},

{
    "location": "API/#RHEOS.dynamicmodelpredict",
    "page": "API",
    "title": "RHEOS.dynamicmodelpredict",
    "category": "function",
    "text": "dynamicmodelpredict(data::RheologyDynamic, model::RheologyModel)\n\nGiven dynamic rheology data and model, return new dataset based on model parameters. Returns another RheologyDynamic instance with the predicted Gp and Gpp based on the frequencies and model parameters.\n\n\n\n\n\n"
},

{
    "location": "API/#Fitting-and-Predicting-Functions-1",
    "page": "API",
    "title": "Fitting and Predicting Functions",
    "category": "section",
    "text": "modelfit\nmodelpredict\nmodelstepfit\nmodelsteppredict\ndynamicmodelfit\ndynamicmodelpredict"
},

{
    "location": "API/#RHEOS.linegen",
    "page": "API",
    "title": "RHEOS.linegen",
    "category": "function",
    "text": "linegen(t_total::Real; stepsize::Real = 1.0)\n\nGenerate RheologyData struct with a simple line loading of height 1.0.\n\nArguments\n\nt_total: Total time length of data\nstepsize: Time sampling period\n\n\n\n\n\n"
},

{
    "location": "API/#RHEOS.stepgen",
    "page": "API",
    "title": "RHEOS.stepgen",
    "category": "function",
    "text": "stepgen(t_total::Real, t_on::Real; t_trans::Real = 0.0, stepsize::Real = 1.0)\n\nGenerate RheologyData struct with a step loading of height 1.0. If t_trans is 0.0 then  the step is instantaneous, otherwise the step is approximated by a logistic function  approximately centered at t_on.\n\nArguments\n\nt_total: Total time length of data\nt_on: Step on time\nt_trans: Step transition time\nstepsize: Time sampling period\n\n\n\n\n\n"
},

{
    "location": "API/#RHEOS.rampgen",
    "page": "API",
    "title": "RHEOS.rampgen",
    "category": "function",
    "text": "rampgen(t_total::Real, t_start::Real, t_stop::Real; amplitude::Real = 1.0, baseval::Real = 0.0, stepsize::Real = 1.0)\n\nGenerate RheologyData struct with a ramp function. Reaches amplitude of 1.0 at t_stop.\n\nArguments\n\nt_total: Total time length of data\nt_start: Time for starting ramp\nt_stop: Time of stopping ramp\nstepsize: Time sampling period\n\n\n\n\n\n"
},

{
    "location": "API/#RHEOS.singen",
    "page": "API",
    "title": "RHEOS.singen",
    "category": "function",
    "text": "singen(t_total::Real, frequency::Real; t_start::Real = 0.0, phase::Real = 0.0, stepsize::Real = 1.0)\n\nGenerate RheologyData struct with a sinusoidal loading of amplitude 1.0.\n\nArguments\n\nt_total: Total time length of data\nfrequency: Frequency of oscillation (Hz)\nt_start: Time for oscillation to begin\nphase: Phase of oscillation (radians)\nstepsize: Time sampling period\n\n\n\n\n\n"
},

{
    "location": "API/#RHEOS.noisegen",
    "page": "API",
    "title": "RHEOS.noisegen",
    "category": "function",
    "text": "noisegen(t_total::Real; seed::Union{Int, Nothing} = nothing, stepsize::Real = 1.0)\n\nGenerate uniform random noise of maximum amplitude +/- 1.0. If reproducibility is required,  always use the same number in the seed keyword argument with the same non-negative integer.\n\nArguments\n\nt_total: Total time length of data\nseed: Seed used for random number generation\nbaseval: Initial amplitude before oscillation started\nstepsize: Time sampling period\n\n\n\n\n\n"
},

{
    "location": "API/#RHEOS.repeatdata",
    "page": "API",
    "title": "RHEOS.repeatdata",
    "category": "function",
    "text": "repeatdata(self::RheologyData, n::Integer)\n\nRepeat a given RheologyData data set n times.\n\n\n\n\n\n"
},

{
    "location": "API/#Data-Generation-Functions-1",
    "page": "API",
    "title": "Data Generation Functions",
    "category": "section",
    "text": "linegen\nstepgen\nrampgen\nsingen\nnoisegen\nrepeatdata"
},

{
    "location": "API/#RHEOS.importdata",
    "page": "API",
    "title": "RHEOS.importdata",
    "category": "function",
    "text": "importdata(colnames::Array{String,1}, filedir::String; delimiter=\',\')\n\nLoad data from a CSV file (two/three columns, comma seperated by default but delimiter can be specified in the delimiter keyword argument). Columns must be identified by providing an array of strings which tell the function which data is contained in each column. \n\nCan be used to construct either a RheologyData instance or a RheologyDynamic instance. Function detects whether \"time\" or \"frequency\" has been included and proceeds accordingly. For oscillatory data, all three columns (Gp, Gpp, Frequency) must be provided. For regular viscoelastic data two or three columns can be provided but time must always be provided. I.e. either stress/strain/time, stress/time or  strain/time.\n\n\n\n\n\n"
},

{
    "location": "API/#RHEOS.exportdata",
    "page": "API",
    "title": "RHEOS.exportdata",
    "category": "function",
    "text": "exportdata(self::Union{RheologyData, RheologyDynamic}, filedir::String; ext = \".csv\", delimiter=\',\')\n\nExport RheologyData or RheologyDynamic type to csv format. Exports three columns in order:  stress, strain, time for standard viscoelastic data. Or: Gp, Gpp, frequency for oscillatory data. File extension can be modified using the ext keyword argument. As with importdata, the delimiter can also be set by keyword argument.\n\nUseful for plotting/analysis in other software.\n\n\n\n\n\n"
},

{
    "location": "API/#RHEOS.savedata",
    "page": "API",
    "title": "RHEOS.savedata",
    "category": "function",
    "text": "savedata(self::Union{RheologyData, RheologyDynamic}, filedir::String; ext = \".jld2\")\n\nSave RheologyData or RheologyDynamic object using JLD2 format reuse in  a later Julia session.\n\n\n\n\n\n"
},

{
    "location": "API/#RHEOS.loaddata",
    "page": "API",
    "title": "RHEOS.loaddata",
    "category": "function",
    "text": "loaddata(filedir::String)\n\nLoads RheologyData or RheologyDynamic object from a jld2 file.\n\n\n\n\n\n"
},

{
    "location": "API/#RHEOS.savemodel",
    "page": "API",
    "title": "RHEOS.savemodel",
    "category": "function",
    "text": "savemodel(self::RheologyModel, filedir::String; ext = \".jld2\")\n\nSave RheologyModel object using JLD2 format reuse in a later Julia session.\n\n\n\n\n\n"
},

{
    "location": "API/#RHEOS.loadmodel",
    "page": "API",
    "title": "RHEOS.loadmodel",
    "category": "function",
    "text": "loaddata(filedir::String)\n\nLoads RheologyModel from a JLD2 file.\n\n\n\n\n\n"
},

{
    "location": "API/#File-I/O-1",
    "page": "API",
    "title": "File I/O",
    "category": "section",
    "text": "importdata\nexportdata\nsavedata\nloaddata\nsavemodel\nloadmodel"
},

]}
