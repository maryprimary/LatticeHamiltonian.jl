#=
HS变换
=#


"""
四个分量的辅助场
"""
struct AuxiliaryField4
    name :: String
    g :: ComplexF64
    alpha :: ComplexF64
    ΔV :: Vector{ComplexF64}
    Coef :: Vector{ComplexF64}
end


"""
四个分量的
"""
function AuxiliaryField4(name::String, g::ComplexF64, alpha::ComplexF64)
    delv = Vector{ComplexF64}(undef, 4)
    coef = Vector{ComplexF64}(undef, 4)
    for (idx, comp) in enumerate([-2, -1, 1, 2])
        delv[idx] = exp(g*FourComponentEta(Val{comp}())) - 1.0
        coef[idx] = FourComponentGamma(Val{comp}())*exp(g*alpha)
    end
    return AuxiliaryField4(
        name, g, alpha, delv, coef
    )
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



"""
四个分量的辅助场
"""
struct AuxiliaryField2
    name :: String
    λ :: Float64
    ΔV :: Matrix{Float64}
    Coef :: Vector{Float64}
end



"""
两个分量的
"""
function AuxiliaryField2(name::String, U::Float64, dtau::Float64)
    delv = Matrix{Float64}(undef, 2, 2)
    coef = Vector{Float64}(undef, 2)
    exph = exp(U*dtau/2)
    lamb = log(exph + sqrt(exph^2 - 1))
    delv[1, 1] = exp(-lamb-U*dtau/2) - 1
    delv[1, 2] = exp(lamb-U*dtau/2) - 1
    delv[2, 1] = exp(lamb-U*dtau/2) - 1
    delv[2, 2] = exp(-lamb-U*dtau/2) - 1
    coef .= 1
    return AuxiliaryField2(
        name, lamb, delv, coef
    )
end



