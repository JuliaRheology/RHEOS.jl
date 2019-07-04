# """
#     function zerotime(self::RheologyData)

# Convenience function to normalize time such that the starting time is 0.0
# """
# function zerotime(self::RheologyData)

#     return RheologyData(self.σ, self.ϵ, self.t .- minimum(self.t), self.sampling, vcat(self.log, ["Normalized time to start at 0.0"]))

# end
