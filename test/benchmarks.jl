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

Hardware:
Intel Core i5-7200U Processor (3MB Cache, up to 3.1Ghz)
8GB LPDDR3 1866

=========================
Test on maths functions 
=========================

Test on the MittalgLeffler function

Version from Julia repo
  57.567 ns (1 allocation: 16 bytes)
Version from RHEOS (should be the same as above)
  58.195 ns (1 allocation: 16 bytes)



Test on the inverse Laplace transform

Talbot algorithm on FractSLS_Zener creep function
  scalar value:   21.206 μs (642 allocations: 18.75 KiB)
  vector using talbot.:   21.383 ms (642143 allocations: 18.33 MiB)
  vector using talbotarr:   3.418 ms (173413 allocations: 3.61 MiB)

Weeks (80) algorithm on FractSLS_Zener creep function
  Initialisation:   105.563 μs (2634 allocations: 110.75 KiB)
  Scalar value:   1.616 μs (1 allocation: 16 bytes)
  Vector values:   1.596 ms (2 allocations: 8.03 KiB)


=========================
Test on Moduli functions 
=========================

Test on the Maxwell model

  67.185 ns (2 allocations: 112 bytes)
  7.559 μs (2 allocations: 8.09 KiB)
  33.592 ns (1 allocation: 16 bytes)
  8.153 μs (1 allocation: 8.00 KiB)

Test on the FractSLS_Zener model relaxation modulus - with mittleff

Direct function evaluation
  207.000 ns (7 allocations: 112 bytes)
  68.707 μs (15 allocations: 8.53 KiB)
From the RheoModelClass
  166.950 ns (4 allocations: 160 bytes)
  64.459 μs (2 allocations: 8.11 KiB)
From the RheoModel
  104.457 ns (1 allocation: 16 bytes)
  59.557 μs (1 allocation: 8.00 KiB)

Test on the FractSLS_Zener model creep compliance - with inverse laplace


function wrappers using Weeks
  62.121 μs (71 allocations: 30.80 KiB)
  1.660 ms (72 allocations: 38.81 KiB)
  1.636 μs (1 allocation: 16 bytes)
  1.610 ms (5 allocations: 8.12 KiB)

From the RheoModelClass
  10.674 μs (101 allocations: 2.64 KiB)
  10.721 ms (99102 allocations: 2.47 MiB)
From the RheoModel
  9.535 μs (1 allocation: 16 bytes)
  9.580 ms (1 allocation: 8.00 KiB)

=#
