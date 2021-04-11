using RHEOS
using BenchmarkTools
import MittagLeffler

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


println("=========================")
println("Test on Moduli functions ")
println("=========================")

println()
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
println("From the RheoModel(Class)")
@btime relaxmod(FractSLS_Zener, t0, [1.0,0.5,2.0,3.0])
@btime relaxmod(FractSLS_Zener, t, [1.0,0.5,2.0,3.0])
m = RheoModel(FractSLS_Zener,  (cₐ = 1.0, a=0.5, kᵦ = 2.0, kᵧ=3.0))
@btime relaxmod(m, t0)
@btime relaxmod(m, t)
nothing


#=
=========================
Test on maths functions 
=========================

Test on the MittalgLeffler function

Version from Julia repo
  56.957 ns (1 allocation: 16 bytes)
Version from RHEOS (should be the same as above)
  56.930 ns (1 allocation: 16 bytes)
=========================
Test on Moduli functions 
=========================


Test on the Maxwell model

  67.689 ns (2 allocations: 112 bytes)
  7.500 μs (3 allocations: 8.12 KiB)
  35.061 ns (1 allocation: 16 bytes)
  6.811 μs (1 allocation: 8.00 KiB)

Test on the FractSLS_Zener model relaxation modulus - with mittleff

Direct function evaluation
  207.834 ns (7 allocations: 112 bytes)
  67.608 μs (15 allocations: 8.53 KiB)
From the RheoModel(Class)
  169.915 ns (4 allocations: 160 bytes)
  65.957 μs (3 allocations: 8.16 KiB)
  103.534 ns (1 allocation: 16 bytes)
  59.479 μs (1 allocation: 8.00 KiB)


=#