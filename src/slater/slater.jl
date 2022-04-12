#=
Slater行列式的模块
=#

module SlaterDeterminant

export Slater, multiply_left!
using LinearAlgebra
using ..LatticeHamiltonian

"""
Slater行列式
"""
struct Slater{T}
    name :: String
    V :: Matrix{T}
end


"""
给Slater行列式左乘
"""
function multiply_left!(mat::Matrix{T}, sla::Slater{T}) where T
    sla.V .= mat * sla.V
end

"""
算符左作用
"""
function multiply_left!(op::DenseFermiOpt{T}, sla::Slater{T}) where T
    sla.V .= op.V * sla.V
end


"""
转置
"""
function Base.:adjoint(sla::Slater{T}) :: Slater{T} where T
    mat = adjoint(sla.V)
    return Slater{T}(sla.name*"'", mat)
end


"""
矩阵乘法 
"""
function Base.:*(sla1::Slater{T}, sla2::Slater{T}) :: Slater{T} where T
    mat = sla1.V * sla2.V
    return Slater{T}(sla1.name*"*"*sla2.name, mat)
end


"""
求逆矩阵
"""
function Base.:inv(sla::Slater{T}) :: Slater{T} where T
    mat = inv(sla.V)
    return Slater{T}("inv("*sla.name*")", mat)
end


"""
求行列式
"""
function LinearAlgebra.:det(sla::Slater{T}) :: T where T
    return det(sla.V)
end

end
