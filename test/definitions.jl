function _operators_logs()
    d1=strainfunction(timeline(),t->exp(-t))
    d2=strainfunction(timeline(),t->1-exp(-t))
    d=2*d1 - (-d2) + d2

    d3=rheologrun(d.log)

    (d3.ϵ == d.ϵ) && all([ abs(e-2.)<=eps(RheoFloat) for e in d.ϵ ])
end
@test _operators_logs()
