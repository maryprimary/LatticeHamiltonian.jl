#=
测试费米子算符
=#

using Test
using LatticeHamiltonian


@testset "算符" begin
    a = DenseFermiOpt{Float64}("dense", [1. 2.; 3. 4.])
    println(a)
    @test 1 == 1
    b = SparseFermiOpt{Float64}("sparse", 4, [1, 2], [1, 2], [0 1; 1 0], 0)
    println(b)
    c = sparse2dense(b)
    println(c)
    println(exp(a))
    println(exp(b))
    println(exp(c))
end

