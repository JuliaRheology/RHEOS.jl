println("===============================================")
println("Testing rheodata.jl")
println("===============================================")


function _RheoTimeData_explicit_nolog()
    a = Vector{RheoFloat}(1.0:1.0:3.0)
    b = Vector{RheoFloat}(4.0:1.0:6.0)
    c = Vector{RheoFloat}(7.0:1.0:9.0)

    data = RheoTimeData(a, b, c, nothing)

    data.σ==a && data.ϵ==b && data.t==c && isnothing(data.log)
end
@test _RheoTimeData_explicit_nolog()

function _RheoTimeData_const_nolog()
    eps = Vector{RheoFloat}(1.0:1.0:3.0)
    sig = Vector{RheoFloat}(4.0:1.0:6.0)
    t = Vector{RheoFloat}(7.0:1.0:9.0)

    data = RheoTimeData(strain=eps, stress=sig, t=t, savelog=false)
    showlog(data)

    getstress(data)==sig && getstrain(data)==eps && gettime(data)==t && isnothing(data.log)
end
@test _RheoTimeData_const_nolog()

function _RheoTimeData_const_nostrain()
    sig = Vector{RheoFloat}(4.0:1.0:6.0)
    t = Vector{RheoFloat}(7.0:1.0:9.0)

    data = RheoTimeData(σ=sig, t=t)
    showlog(data)

    data.σ==sig && data.ϵ==RheoFloat[] && data.t==t && data.log[1].info.type == stress_only::TimeDataType
end
@test _RheoTimeData_const_nostrain()

function _operators_logs()
    d1=strainfunction(timeline(0:0.5:10),t->exp(-t))
    d2=strainfunction(timeline(0:0.5:10),t->1-exp(-t))
    d=2*d1 - (-d2) + d2

    d3=rheologrun(d.log)

    (d3.ϵ == d.ϵ) && all([ abs(e-2.)<=eps(RheoFloat) for e in d.ϵ ])
end
@test _operators_logs()



function _RheoFreqData_explicit_nolog()
    a = Vector{RheoFloat}(1.0:1.0:3.0)
    b = Vector{RheoFloat}(4.0:1.0:6.0)
    c = Vector{RheoFloat}(7.0:1.0:9.0)

    data = RheoFreqData(a, b, c, nothing)

    data.Gp==a && data.Gpp==b && data.ω==c && isnothing(data.log)
end
@test _RheoFreqData_explicit_nolog()

function _RheoFreqData_const_nolog()
    a = Vector{RheoFloat}(1.0:1.0:3.0)
    b = Vector{RheoFloat}(4.0:1.0:6.0)
    c = Vector{RheoFloat}(7.0:1.0:9.0)

    data = RheoFreqData(omega=a, Gp=b, Gpp=c , savelog=false)
    showlog(data)

    getfreq(data)==a && getstorage(data)==b && getloss(data)==c && isnothing(data.log)
end
@test _RheoFreqData_const_nolog()



function _rheoconvert()
    vi64=Int64(1)
    ai64=[vi64,vi64]
    arf=Vector{RheoFloat}([1,2,3])
    rheoconvert(RheoFloat(1.0))===rheoconvert(vi64) && typeof(rheoconvert(ai64)) == Vector{RheoFloat} &&
        arf === rheoconvert(arf)  &&  ai64 !== rheoconvert(ai64)
end
@test _rheoconvert()


function _union_strain1stress2()
    time_instance = timeline(t_start=0.0, t_end=15.0, step=0.2)
    imposed_strain = strainfunction(time_instance, sawtooth(offset=5.0, amp=2, period=5))
    imposed_stress = stressfunction(time_instance, stairs(offset=3.0, amp=0.5, width=4.0))

    combined = imposed_strain|imposed_stress

    (combined.σ==imposed_stress.σ) && (combined.ϵ==imposed_strain.ϵ) && (combined.t==imposed_stress.t)
end
@test _union_strain1stress2()

function _union_stress1strain2()
    time_instance = timeline(t_start=0.0, t_end=15.0, step=0.2)
    imposed_strain = strainfunction(time_instance, sawtooth(offset=5.0, amp=2, period=5))
    imposed_stress = stressfunction(time_instance, stairs(offset=3.0, amp=0.5, width=4.0))

    combined = imposed_stress|imposed_strain

    (combined.σ==imposed_stress.σ) && (combined.ϵ==imposed_strain.ϵ) && (combined.t==imposed_stress.t)
end
@test _union_stress1strain2()



function __setdata!()
  t=timeline(0:0.1:1)
  RHEOS._setdata!(t.ϵ,t.t)
  b = t.ϵ==t.t
  RHEOS._setdata!(t.ϵ,2 .* t.t)
  b && (t.ϵ== 2 .* t.t)
end
@test __setdata!()


function __mapdata!()
  t=timeline(0:0.1:1)
  RHEOS._mapdata!(x->x,t.ϵ,t.t)
  b = t.ϵ==t.t
  RHEOS._mapdata!(x->2*x, t.ϵ, t.t)
  b && (t.ϵ== 2 .* t.t)
end
@test __mapdata!()


using RHEOS

# Step 1: Generate Timeline
datat = timeline(t_start = 0, t_end = 20.0, step = 0.02)  # Create a timeline from 0 to 20 seconds with a step size of 0.02 seconds
rheotimedatatype(datat)  # Ensure the timeline data type is correct for RHEOS

# Step 2: Generate Strain Data (Ramp & hold)
dramp_strain = strainfunction(datat, ramp(offset = 1.0, gradient = 0.8))  # Generate a ramp strain function with offset 1.0 and gradient 0.8
dhold_strain = dramp_strain - strainfunction(datat, ramp(offset = 3.0, gradient = 0.8))  # Generate a hold strain function by subtracting a shifted ramp

# Step 3: Generate Stress Data (Ramp & hold)
dramp_stress = stressfunction(datat, ramp(offset = 4.0, gradient = 0.8))  # Generate a ramp stress function with offset 4.0 and gradient 0.8
dhold_stress = dramp_stress - stressfunction(datat, ramp(offset = 5.0, gradient = 0.8))  # Generate a hold stress function by subtracting a shifted ramp

# Define the rheological model
model = RheoModel(SLS_Zener, (η = 1, kᵦ = 1, kᵧ = 1))
data_ext = dhold_stress
rheotimedatatype(data_ext)
SLS_predict = modelpredict(data_ext, model)
data = SLS_predict

# Fit the model to the data
SLS_Zener_model = modelfit(data, SLS_Zener, strain_imposed)

# Define the test function
function _test_extractfitdata()
    # Call the extractfitdata function with data.log
    extracted_data = extractfitdata(data.log)
    all_tests_passed = true
    
    # Iterate through data.log to dynamically verify modelfit entries
    for (index, log_entry) in enumerate(data.log)
        if log_entry.action.funct == :modelfit
            model_name = log_entry.info.model_name
            expected_params = log_entry.info.model_params
            expected_info = log_entry.info

            if haskey(extracted_data, model_name)
                model_entries = extracted_data[model_name]

                # Find the corresponding entry in the extracted data
                matching_entries = filter(x -> x.index == index, model_entries)

                if length(matching_entries) == 1
                    extracted_entry = matching_entries[1]

                    if extracted_entry.params == expected_params && extracted_entry.info == expected_info && extracted_entry.index == index
                        println("Test passed for entry $index: Extracted data matches expected data.")
                    else
                        println("Test failed for entry $index: Extracted data does not match expected data.")
                        println("Extracted params: $(extracted_entry.params)")
                        println("Expected params: $expected_params")
                        println("Extracted info: $(extracted_entry.info)")
                        println("Expected info: $expected_info")
                        println("Extracted index: $(extracted_entry.index)")
                        println("Expected index: $index")
                        all_tests_passed = false
                    end
                else
                    println("Test failed for entry $index: Number of matching entries does not match expected.")
                    println("Matching entries: $matching_entries")
                    all_tests_passed = false
                end
            else
                println("Test failed for entry $index: $model_name model not found in extracted data.")
                println("Extracted data: $extracted_data")
                all_tests_passed = false
            end
        end
    end
    
    return all_tests_passed
end

# Run the test function
@test _test_extractfitdata()
