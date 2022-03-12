#=
测试slater相关的功能
=#

using Test
using LatticeHamiltonian.SlaterDeterminant
using LinearAlgebra


@testset "Slater" begin
    sla1 = Slater{ComplexF64}("ϕ", rand(ComplexF64, 4, 4))
    @show sla1
    @test det(sla1) == det(sla1.V)
    sla2 = Slater{ComplexF64}("ψ", rand(ComplexF64, 4, 4))
    sla3 = Slater{ComplexF64}("τ", copy(sla2.V))
    @test all(isapprox.(sla2.V, sla3.V))
    multiply_left!(sla1.V, sla3)
    sla4 = sla1*sla2
    @test all(isapprox.(sla3.V, sla4.V))
    mat = rand(Float64, 4, 4)
    mat = mat + mat'
    sla1 = Slater{Float64}("V", eigvecs(mat))
    sla2 = sla1'
    sla3 = inv(sla1)
    @test all(isapprox.(sla2.V, sla3.V))
end

