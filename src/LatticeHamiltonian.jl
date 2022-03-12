module LatticeHamiltonian


include("operators/operators.jl")
export FermiHamiltonian, DenseFermiOpt, SparseFermiOpt
export sparse2dense
export AuxiliaryField, FourComponentGamma, FourComponentEta

#所有的Vector都要有明确的类型，不能是抽象类型
struct FermiHamiltonian{LOPT, SOPT}
    LinearOps :: Vector{LOPT}
    LinearCoe :: Vector{Float64}
    SquareOps :: Vector{SOPT}
    SquareCoe :: Vector{Float64}
end



include("slater/slater.jl")
#using .SlaterDeterminant
#export Slater



end # module
