tol = (eps(RHEOS.RheoFloat))^(0.125) 

function _timeline()
    time_instance = timeline(t_start=0.0, t_end=15.0, step=0.1)

    typeof(time_instance)==RheoTimeData && time_instance.t==(0.0:0.1:15.0)
end
@test _timeline()

function _strainfunction()
    time_instance = timeline(t_start=0.0, t_end=15.0, step=0.1)

    strain_imposed = strainfunction(time_instance, t -> t)

    typeof(strain_imposed)==RheoTimeData && strain_imposed.t==(0.0:0.1:15.0) && strain_imposed.ϵ==(0.0:0.1:15.0)
end
@test _strainfunction()

function _stressfunction()
    time_instance = timeline(t_start=0.0, t_end=15.0, step=0.1)

    stress_imposed = stressfunction(time_instance, t -> t)

    typeof(stress_imposed)==RheoTimeData && stress_imposed.t==(0.0:0.1:15.0) && stress_imposed.σ==(0.0:0.1:15.0)
end
@test _stressfunction()

function _hstep()
    time_instance = timeline(t_start=0.0, t_end=15.0, step=0.1)
    step_stress_imposed = stressfunction(time_instance, hstep(offset=1.0, amp=2.0))

    step_stress_imposed.σ==(0.0:0.1:15.0)
end
@test _hstep()