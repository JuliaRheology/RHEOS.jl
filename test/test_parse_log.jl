using Test

include("parse_log.jl")

# Mock log data for testing
const mock_log = [
    (action = (type = :source, funct = :importcsv), info = ()),
    (action = (type = :analysis, funct = :modelfit), info = (model_name = "maxwell", model_params = Dict("η" => 76.30267813573391, "k" => 0.5451389878296595), error = 4.234991209005228)),
    (action = (type = :analysis, funct = :modelfit), info = (model_name = "SLS_Zener", model_params = Dict("η" => 0.2995472252368927, "kᵦ" => 1.0333390861749647, "kᵧ" => 0.5), error = 3.134312236497204e-11)),
]

# Define the test function
function test_parse_log()
    expected_output = Dict(
        "maxwell" => [
            Dict("params" => Dict("η" => 76.30267813573391, "k" => 0.5451389878296595), "error" => 4.234991209005228, "index" => 2),
        ],
        "SLS_Zener" => [
            Dict("params" => Dict("η" => 0.2995472252368927, "kᵦ" => 1.0333390861749647, "kᵧ" => 0.5), "error" => 3.134312236497204e-11, "index" => 3),
        ]
    )

    actual_output = parse_log(mock_log)

    @test actual_output == expected_output
end

# Run tests
@testset "Tests for parse_log function" begin
    @testset "Basic functionality" begin
        test_parse_log()
    end
end
