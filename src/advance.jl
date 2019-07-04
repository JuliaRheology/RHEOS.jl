"""
    function zerotime(self::RheologyData)

Convenience function to normalize time such that the starting time is 0.0
"""
function zerotime(self::RheologyData)

    return RheologyData(self.σ, self.ϵ, self.t .- minimum(self.t), self.sampling, vcat(self.log, ["Normalized time to start at 0.0"]))

end

function sampleratecompare(t1::Vector{RheoFloat}, t2::Vector{RheoFloat})

    @assert constantcheck(t1) "Sample-rate of both arguments must be constant, first argument is non-constant"
    @assert constantcheck(t2) "Sample-rate of both arguments must be constant, second argument is non-constant"

    diff1 = getsampleperiod(t1)
    diff2 = getsampleperiod(t2)

    diff1 == diff2

end