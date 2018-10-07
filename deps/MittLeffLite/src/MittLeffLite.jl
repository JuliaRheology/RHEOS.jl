__precompile__()

module MittLeffLite

import SpecialFunctions.gamma

export mittleff
export mittlefferr
export mittleffderiv

include("mittag_leffler.jl")

end # module
