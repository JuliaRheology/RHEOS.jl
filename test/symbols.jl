println("Testing symbol.jl")
println("===============================================")


function _symbol_to_unicode()
    RHEOS.symbol_to_unicode(:ε) == :ϵ
end
@test _symbol_to_unicode()



function _symbol_to_unicode_tuple()
    RHEOS.symbol_to_unicode((ε = [1,2,3], sigma = [3,4,5], time = [0.1, 0.2, 0.3], eta = 2.5, k = 6)) == (ϵ = [1, 2, 3], σ = [3, 4, 5], t = [0.1, 0.2, 0.3], η = 2.5, k = 6) 
end
@test _symbol_to_unicode_tuple()

function _symbol_to_unicode_dict()
    ds=Dict(:beta=>0.05, :c_beta=>0.05, :eta=>10.0)
    nt=RHEOS.symbol_to_unicode(ds) 
    nt.η==10.0
end
@test _symbol_to_unicode_dict()


println("===============================================")





