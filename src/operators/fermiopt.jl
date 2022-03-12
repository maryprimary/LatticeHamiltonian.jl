#=
费米子算符
=#


abstract type AbstractFermiOpt{T} <: AbstractQtmOpt{T} end

"""
包含完整矩阵的算符
"""
struct DenseFermiOpt{T} <: AbstractFermiOpt{T}
    name :: String
    V :: Matrix{T}
end


"""
只标记了指标的算符
"""
struct SparseFermiOpt{T} <: AbstractFermiOpt{T}
    name :: String
    dim :: Int64
    rowi :: Vector{Int64}
    coli :: Vector{Int64}
    V :: Matrix{T}
    diag :: T
end



function Base.:show(::IO, op::DenseFermiOpt)
    println(op.name)
    println(op.V)
end


function Base.:show(::IO, op::SparseFermiOpt)
    println(op.name)
    println(op.dim)
    println(op.rowi)
    println(op.coli)
    println(op.V)
    println(op.diag)
end



"""
将sparse转成dense
"""
function sparse2dense(op::SparseFermiOpt{T}) :: DenseFermiOpt{T} where T
    mat = zeros(T, op.dim, op.dim)
    for didx in 1:1:op.dim
        mat[didx, didx] = op.diag
    end
    for idxs in CartesianIndices(op.V)
        ri, ci = Tuple(idxs)
        mat[op.rowi[ri], op.coli[ci]] = op.V[ri, ci]
    end
    return DenseFermiOpt{T}(op.name, mat)
end


"""
给算符exp
"""
function Base.:exp(op::DenseFermiOpt{T}) :: DenseFermiOpt{T} where T
    mat = exp(op.V)
    return DenseFermiOpt{T}("exp^"*op.name, mat)
end


"""
给算符exp
"""
function Base.:exp(op::SparseFermiOpt{T}) :: SparseFermiOpt{T} where T
    mat = exp(op.V)
    diag = exp(op.diag)
    return SparseFermiOpt{T}(
        "exp^"*op.name, op.dim,
        op.rowi, op.coli,
        mat, diag
    )
end

