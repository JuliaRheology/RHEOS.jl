println("===============================================")
println("Testing rheomodel.jl")
println("===============================================")




function _freezeparams()
    # SLS2_mod = freezeparams( SLS2, G₀=2, η₂=3.5)
    SLS2_mod = freezeparams( SLS2,(G₀=2,η₂=3.5))

    relaxmod(SLS2, 1, [2,1,2,3,3.5]) == relaxmod(SLS2_mod, 1, [1,2,3])
end
@test _freezeparams()


function _getparams()
    m = RheoModel(Maxwell, k=1, eta=2)
    p1 = getparams(m)
    p2 = getparams(m,unicode=false)
    (p1.k==1.0) && (p1.η==2.0) &&  (p2.k==1.0) && (p2.eta==2.0) 
end
@test _getparams()

function _getparams_dict()
    m = RheoModel(Maxwell, k=1, eta=2)
    p = dict(getparams(m,unicode=false))
    (p[:k]==1.0) && (p[:eta]==2.0) 
end
@test _getparams_dict()



#
#   The moduli functions on arrays are extensively tested as part of the model tests
#


function _scalar_moduli()
    m=RheoModel(Spring, k=2)
    (relaxmod(m))(1) == relaxmod(m, 1) == relaxmod(Spring, k=2, 1) == 2. &&
        (creepcomp(m))(1) == creepcomp(m, 1) == creepcomp(Spring, k=2, 1) == 0.5 &&
        (storagemod(m))(1) == storagemod(m, 1) == storagemod(Spring, k=2, 1) == 2. &&
        (lossmod(m))(1) == lossmod(m, 1) == lossmod(Spring, k=2, 1) == 0.0 &&
        (dynamicmod(m))(1) == dynamicmod(m, 1) == dynamicmod(Spring, k=2, 1) == 2.0 + 0.0*im 
end
@test _scalar_moduli()


function _builddiffequation()
    equation = (ϵ =((:cₐ,:a),), σ =((1.0,0.0),(:(cₐ/cᵦ),:(a-β))))
    p = (:cₐ, :a, :cᵦ, :β)
    built= RHEOS.builddiffequation(equation,p)
    @test isa(built, RHEOS.DiffEqu)
    @test built.leftvar == :ϵ
    @test built.rightvar == :σ

    @test built.leftde isa Set{RHEOS.DETerm{RHEOS.DiffScaFree}}
    @test built.rightde isa Set{RHEOS.DETerm{RHEOS.DiffScaFree}}
    @test length(built.leftde) == 1
    @test length(built.rightde) == 2

    p_fixed = [0.5, 0.6, 0.5, 0.35]

    built_fixed = RHEOS._builddiffequation(built, p_fixed)

    left = Set{RHEOS.DETerm{RheoFloat}}()
    push!(left, RHEOS.DETerm{RheoFloat}(0.5,0.6))

    right = Set{RHEOS.DETerm{RheoFloat}}()
    push!(right, RHEOS.DETerm{RheoFloat}(1.0,0.0))
    push!(right, RHEOS.DETerm{RheoFloat}(1.0,0.25))
    
    


    @test left == built_fixed.leftde
    @test right == built_fixed.rightde
end

_builddiffequation()
