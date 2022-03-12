#=
HS变换
=#


"""
辅助场
"""
struct AuxiliaryField
    name :: String
    g :: ComplexF64
    type :: Symbol
    #AuxiliaryField(gamma, eta) = begin
    #    @assert length(gamma) ==  length(eta)
    #    new(gamma, eta)
    #end
end


"""
四分量辅助场
"""
FourComponentGamma(::Val{-2}) = 1 - 6/√3
FourComponentGamma(::Val{-1}) = 1 + 6/√3
FourComponentGamma(::Val{1}) = 1 + 6/√3
FourComponentGamma(::Val{2}) = 1 - 6/√3

"""
四分量辅助场
"""
FourComponentEta(::Val{-2}) = -√(2(3+√6))
FourComponentEta(::Val{-1}) = -√(2(3-√6))
FourComponentEta(::Val{1}) = √(2(3-√6))
FourComponentEta(::Val{2}) = √(2(3+√6))


