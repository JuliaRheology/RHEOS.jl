function _timeline()
    time_instance = timeline(t_start=0.0, t_end=15.0, step=0.1)

    typeof(time_instance)==RheoTimeData && time_instance.t==(0.0:0.1:15.0)
end
@test _timeline()

function _frequencyspec()
    freq_instance = frequencyspec(ω_start=0.0, ω_end=15.0, step=0.1)

    typeof(freq_instance)==RheoFreqData && freq_instance.ω==(0.0:0.1:15.0)
end
@test _frequencyspec()

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
    imposed = stressfunction(time_instance, hstep(offset=1.0, amp=2.0))

    all(v -> v==0, imposed.σ[1:10]) && all(v -> v==2.0, imposed.σ[11:end])
end
@test _hstep()

function _ramp()
    dt = 0.5
    _offset = 1.0
    time_instance = timeline(t_start=0.0, t_end=15.0, step=dt)
    imposed = stressfunction(time_instance, ramp(offset=_offset, gradient=1.0))

    all(v -> v==0, imposed.σ[imposed.t.<_offset]) && imposed.σ[imposed.t.==(_offset+dt)]≈[dt] && imposed.σ[end]≈14.0
end
@test _ramp()

function _stairs()
    time_instance = timeline(t_start=0.0, t_end=15.0, step=0.2)
    imposed = stressfunction(time_instance, stairs(offset=3.0, amp=0.5, width=4.0))

    test1 = all(v -> v==0.0, imposed.σ[imposed.t.<3])
    test2 = all(v -> v==0.5, imposed.σ[(imposed.t.>=3) .& (imposed.t.<7)])
    test3 = all(v -> v==1.0, imposed.σ[(imposed.t.>=7) .& (imposed.t.<11)])
    test4 = all(v -> v==1.5, imposed.σ[(imposed.t.>=11) .& (imposed.t.<15)])
    test5 = all(v -> v==2.0, imposed.σ[imposed.t.>=15])
    
    test1 && test2 && test3 && test4 && test5
end
@test _stairs()

function _square()
    dt = 0.1
    time_instance = timeline(t_start=0.0, t_end=3.0, step=dt)
    imposed = stressfunction(time_instance, square(offset=1.0, amp=1.0, period=1.0, width=0.5))

    testresult = all(v -> v==0, imposed.σ[imposed.t.<1.0])
    ons = [1.0, 1.1, 1.2, 1.3, 1.4, 2.0, 2.1, 2.2, 2.3, 2.4, 3.0]
    offs = [1.5, 1.6, 1.7, 1.8, 1.9, 2.5, 2.6, 2.7, 2.8, 2.9]
    for (i,t) in enumerate(imposed.t)
        if t in ons && imposed.σ==0.0
            testresult = testresult && false
        elseif t in offs && imposed.σ==1.0
            testresult = testresult && false
        end
    end
    testresult
end
@test _square()

function _sawtooth()
    time_instance = timeline(t_start=0.0, t_end=5.0, step=0.1)
    _offset = 1.0
    imposed = stressfunction(time_instance, sawtooth(offset=_offset, amp=1.0, period=1.0))
    # this test could be improved
    all(v -> v==0, imposed.σ[imposed.t.<_offset]) && all(v -> v==0.0, imposed.σ[imposed.t.%1.0.==0.0])
end
@test _sawtooth()

function _triangle()
    time_instance = timeline(t_start=0.0, t_end=5.0, step=0.1)
    _amp = 1.0
    _offset = 1.0
    imposed = stressfunction(time_instance, triangle(offset=_offset, amp=1.0, period=1.0, width=0.5))

    all(v -> v==0, imposed.σ[imposed.t.<_offset]) && all(v -> v==0.0, imposed.σ[imposed.t.%1.0.==0.0]) && all(v -> v==_amp, imposed.σ[((imposed.t.%1.0).==0.5) .& (imposed.t.>_offset)])
end
@test _triangle()
