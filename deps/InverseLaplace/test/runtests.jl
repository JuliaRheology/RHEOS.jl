using InverseLaplace
import InverseLaplace: talbot, gwr

using Base.Test

include("weeks_test.jl")
include("interface_test.jl")

#### plain function interface

@test isapprox( talbot(s -> 1/s,1.0) , 1.000000000007737)
@test isapprox( talbot(s -> 1/s,1) , 1.0)

# ilt wraps talbot
@test isapprox( talbot(s -> 1/s,1), ilt(s -> 1/s,1))

@test isapprox( gwr(s -> 1/s,1) , 1.0)
@test isapprox( gwr(s -> 1/s,1.0) , 1.0000000002465594)
@test isapprox( talbot(s -> s/(1+s^2),1) , cos(1))
@test isapprox( talbot(s -> 1/(1+s),1) , exp(-1))
@test isapprox( gwr(s -> 1/(1+s),1) , exp(-1))

if Int != Int32 && VERSION >= v"0.5.0"
    @test isapprox( talbot(s -> 1/s^4,1//10) * 6 , 1e-3)
    @test isapprox( gwr(s -> 1/s^4,1//10) * 6 , 0.0010000012595365085; atol = 1e-9)

    # At one point, talbot returned rational with rational input. Now we convert to BigFloat
    @test isapprox( talbot.(s -> 1/s^3,[1//10,2//10], 64) , [5e-3, 2e-2]; atol = 1e-18)
    @test (talbot.(s -> 1/s^3,[[.1,.2] ; [.4,.5]]); true)
    @test isapprox( talbot.(s -> 1/s^3, [.1 .2; .3 .4] ),  [0.005 0.02; 0.045 0.08], atol = 1e-13)
else
    @test isapprox( talbot(s -> 1/s^4,1//10) * 6 , 1e-3; atol = 1e-9)
    @test isapprox( gwr(s -> 1/s^4,1//10) * 6 , 0.0010000012595365085; atol = 1e-9)
    @test isapprox( talbot.(s -> 1/s^3, [.1 .2; .3 .4] ),  [0.005 0.02; 0.045 0.08]; atol = 1e-9)
end

@test isapprox( talbot.(s -> 1/s^3,[.1,.2]) , [0.005,0.02])

# gwr is often less accurate.

if VERSION >= v"0.5.0"
    @test isapprox( gwr.(s -> 1/s^3,[.1,.2]) , [0.005,0.02]; atol = 1e-10)
else
    @test isapprox( gwr.(s -> 1/s^3,[.1,.2]) , [0.005,0.02]; atol = 1e-7)
end
   
@test isapprox( InverseLaplace.talbotarr(s -> 1/s^3,[.1,.2]) , [0.005,0.02])

@test isapprox( map(x -> convert(Float64,x) , talbot.(s -> 1/s^3,[BigFloat(1//10) , BigFloat(2//10)])), [0.005,0.02])

@test typeof(talbot(s -> 1/s^3,1//10)) == BigFloat
@test typeof(talbot(s -> 1/s^3,1//10, 64)) == BigFloat

