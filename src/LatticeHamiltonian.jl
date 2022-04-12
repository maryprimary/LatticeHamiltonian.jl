module LatticeHamiltonian


include("operators/operators.jl")
export FermiHamiltonian, DenseFermiOpt, SparseFermiOpt
export sparse2dense
export AuxiliaryField4, AuxiliaryField2#, FourComponentGamma, FourComponentEta
export fourier_matrix, rlmat2ftmat, ftmat2rlmat
export rlvec2ftvec, ftvec2rlvec, fourier_matrix_pts
export ftint2rlint, rlint2ftint

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


include("interactions/bogoliubov.jl")
export create_triangluar_kpair, triangle_ground_state_construct
export sort_triangluar_kpair


include("interactions/hartreefock.jl")
export triangle_k4tab, mean_field_slater, mean_field_parameters, mean_field_hamiltonian
export mean_field_hamiltonian_rl


end # module
