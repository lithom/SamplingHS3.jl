using SamplingHS3
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

# write your own tests here
@testset "HR Tests" begin include("test_hr.jl") end
@testset "MH Tests" begin include("test_mh.jl") end
@testset "NS Integration" begin include("test_ns_integration_01.jl") end
@testset "SamplingHS3 Utils 01" begin include("test_utils_01.jl") end

#@testset "HR Tests" begin include("test_hr.jl") end
