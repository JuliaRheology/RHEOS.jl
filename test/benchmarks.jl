using RHEOS
using BenchmarkTools
using InverseLaplace
import MittagLeffler
using FunctionWrappers: FunctionWrapper

t=RheoFloat.(0.0:0.001:1.0);
t0=1.0

println("=========================")
println("Test on maths functions ")
println("=========================")

println()
println("Test on the MittalgLeffler function")
println()
println("Version from Julia repo")
@btime MittagLeffler.mittleff(0.5, t0)
println("Version from RHEOS (should be the same as above)")
@btime RHEOS.mittleff(0.5, t0)

println()
println()

println()
println("Test on the inverse Laplace transform")
println()
cₐ = 1.0; a=0.5; kᵦ = 2.0; kᵧ=3.0
println("Talbot algorithm on FractSLS_Zener creep function")
Ĵ(s) = (1/s)*(cₐ*s^a + kᵦ)/(cₐ*s^a*kᵦ + kᵧ*(cₐ*s^a + kᵦ))
print("  scalar value: ")
@btime InverseLaplace.talbot(s->(1/s)*(cₐ*s^a + kᵦ)/(cₐ*s^a*kᵦ + kᵧ*(cₐ*s^a + kᵦ)), t0)
print("  vector using talbot.: ")
@btime InverseLaplace.talbot.(s->(1/s)*(cₐ*s^a + kᵦ)/(cₐ*s^a*kᵦ + kᵧ*(cₐ*s^a + kᵦ)), t)
print("  vector using talbotarr: ")
@btime InverseLaplace.talbotarr(s -> Ĵ(s), t)
println()

println("Weeks (80) algorithm on FractSLS_Zener creep function")
print("  Initialisation: ")
@btime InverseLaplace.Weeks(Ĵ,80)
ft = InverseLaplace.Weeks(Ĵ,80)
print("  Scalar value: ")
@btime ft(t0)
print("  Vector values: ")
@btime ft.(t)
println()
println()

println("=========================")
println("Test on Moduli functions ")
println("=========================")

println()
println("Test on the Maxwell model")
println()
@btime relaxmod(Maxwell, t0, [2.0,3.0])
@btime relaxmod(Maxwell, t, [2.0,3.0]);

m = RheoModel(Maxwell,  (η=2.0, k=3.0))
@btime relaxmod(m, t0)
@btime relaxmod(m, t);


println()
println("Test on the FractSLS_Zener model relaxation modulus - with mittleff")
println()
cₐ = 1.0; a=0.5; kᵦ = 2.0; kᵧ=3.0
println("Direct function evaluation")
@btime kᵦ*MittagLeffler.mittleff(a, -(kᵦ/cₐ)*t0^a) + kᵧ
@btime kᵦ.*MittagLeffler.mittleff.(a, -(kᵦ./cₐ).*t.^a) .+ kᵧ;
println("From the RheoModelClass")
@btime relaxmod(FractSLS_Zener, t0, [1.0,0.5,2.0,3.0])
@btime relaxmod(FractSLS_Zener, t, [1.0,0.5,2.0,3.0])
m = RheoModel(FractSLS_Zener,  (cₐ = 1.0, a=0.5, kᵦ = 2.0, kᵧ=3.0))
println("From the RheoModel")
@btime relaxmod(m, t0)
@btime relaxmod(m, t)

println()
println("Test on the FractSLS_Zener model creep compliance - with inverse laplace")
println()

fwp = ((t,p)->(InverseLaplace.Weeks(s->(1/s)*(p[1]*s^p[2] + p[3])/(p[1]*s^p[2]*p[3] + p[4]*(p[1]*s^p[2] + p[3])),80))(t)) |> FunctionWrapper{RheoFloat,Tuple{RheoFloat,Vector{RheoFloat}}}
fwpa = ((t,p)->(InverseLaplace.Weeks(s->(1/s)*(p[1]*s^p[2] + p[3])/(p[1]*s^p[2]*p[3] + p[4]*(p[1]*s^p[2] + p[3])),80)).(t)) |> FunctionWrapper{Vector{RheoFloat},Tuple{Vector{RheoFloat},Vector{RheoFloat}}}
fw = InverseLaplace.Weeks(s->(1/s)*(1.0*s^0.5 + 2.0)/(1.0*s^0.5*2.0 + 3.0*(1.0*s^0.5 + 2.0)),80) |> FunctionWrapper{RheoFloat,Tuple{RheoFloat}}
fwa = t->fw.(t) |> FunctionWrapper{Vector{RheoFloat},Tuple{Vector{RheoFloat}}}
println()
println("function wrappers using Weeks")
@btime fwp(t0,[1.0,0.5,2.0,3.0])
@btime fwpa(t,[1.0,0.5,2.0,3.0])
@btime fw(t0)
@btime fwa(t)

println()
println("From the RheoModelClass")
@btime creepcomp(FractSLS_Zener, t0, [1.0,0.5,2.0,3.0])
@btime creepcomp(FractSLS_Zener, t, [1.0,0.5,2.0,3.0])
m = RheoModel(FractSLS_Zener,  (cₐ = 1.0, a=0.5, kᵦ = 2.0, kᵧ=3.0))
println("From the RheoModel")
@btime creepcomp(m, t0)
@btime creepcomp(m, t)
nothing




#=

=========================
Test on maths functions 
=========================

Test on the MittalgLeffler function

Version from Julia repo
  58.545 ns (1 allocation: 16 bytes)
Version from RHEOS (should be the same as above)
  58.453 ns (1 allocation: 16 bytes)



Test on the inverse Laplace transform

Talbot algorithm on FractSLS_Zener creep function
  scalar value:   20.937 μs (642 allocations: 18.75 KiB)
  vector using talbot.:   21.988 ms (642143 allocations: 18.33 MiB)
  vector using talbotarr:   3.649 ms (173413 allocations: 3.61 MiB)

Weeks (80) algorithm on FractSLS_Zener creep function
  Initialisation:   106.554 μs (2634 allocations: 110.75 KiB)
  Scalar value:   1.720 μs (1 allocation: 16 bytes)
  Vector values:   1.641 ms (2 allocations: 8.03 KiB)


=========================
Test on Moduli functions 
=========================

Test on the Maxwell model

  67.570 ns (2 allocations: 112 bytes)
  7.409 μs (3 allocations: 8.12 KiB)
  33.446 ns (1 allocation: 16 bytes)
  7.451 μs (1 allocation: 8.00 KiB)

Test on the FractSLS_Zener model relaxation modulus - with mittleff

Direct function evaluation
  209.077 ns (7 allocations: 112 bytes)
  68.211 μs (15 allocations: 8.53 KiB)
From the RheoModelClass
  168.912 ns (4 allocations: 160 bytes)
  64.546 μs (3 allocations: 8.16 KiB)
From the RheoModel
  105.382 ns (1 allocation: 16 bytes)
  59.494 μs (1 allocation: 8.00 KiB)

Test on the FractSLS_Zener model creep compliance - with inverse laplace


function wrappers using Weeks
  62.336 μs (71 allocations: 30.80 KiB)
  1.708 ms (72 allocations: 38.81 KiB)
  1.720 μs (1 allocation: 16 bytes)
  1.654 ms (5 allocations: 8.12 KiB)

From the RheoModelClass
  10.712 μs (101 allocations: 3.67 KiB)
  10.987 ms (99102 allocations: 3.48 MiB)
From the RheoModel
  9.578 μs (1 allocation: 16 bytes)
  9.682 ms (1 allocation: 8.00 KiB)
=#