#=
利用Hatree-Fock整理哈密顿量
=#

using ..SlaterDeterminant
using LinearAlgebra


"""
获得平均场的解
"""
function mean_field_slater(hamup::Matrix{Float64}, hamdn::Matrix{Float64}, N)
    #println(ham)
    eigv = eigvecs(hamup)
    #println(eigv)
    slaup = Slater{Float64}("mfsla_up", real(copy(eigv[:, 1:N])))
    #
    eigv = eigvecs(hamdn)
    sladn = Slater{Float64}("mfsla_dn", real(copy(eigv[:, 1:N])))
    return slaup, sladn
end



"""
获取平均场的参数
"""
function mean_field_parameters(slaup::Slater{Float64}, sladn::Slater{Float64})
    eqgrup = slaup.V * adjoint(slaup.V)
    eqgrdn = sladn.V * adjoint(sladn.V)
    return eqgrup, eqgrdn
end


"""
重新构造平均场哈密顿量
"""
function mean_field_hamiltonian(hamup, hamdn, Γ4, k4tab, eqgrup, eqgrdn)
    hsize = size(Γ4)[1]
    hintup = zeros(Float64, hsize, hsize)
    hintdn = zeros(Float64, hsize, hsize)
    for k1=1:1:hsize[1]; for k2=1:1:hsize[1]; for k3=1:1:hsize[1]
        k4 = k4tab[k1, k2, k3]
        #c+_1up c+_2up c_3up c_4up
        hintup[k1, k4] += 0.5*eqgrup[k2, k3] * Γ4[k1, k2, k3]
        hintup[k2, k3] += 0.5*eqgrup[k1, k4] * Γ4[k1, k2, k3]
        hintup[k1, k3] -= 0.5*eqgrup[k2, k4] * Γ4[k1, k2, k3]
        hintup[k2, k4] -= 0.5*eqgrup[k1, k3] * Γ4[k1, k2, k3]
        #c+_1up c+_2dn c_3dn c_4up
        hintup[k1, k4] += 0.5*eqgrdn[k2, k3] * Γ4[k1, k2, k3]
        hintdn[k2, k3] += 0.5*eqgrup[k1, k4] * Γ4[k1, k2, k3]
        #c+_1dn c+_2up c_3up c_4dn
        hintdn[k1, k4] += 0.5*eqgrup[k2, k3] * Γ4[k1, k2, k3]
        hintup[k2, k3] += 0.5*eqgrdn[k1, k4] * Γ4[k1, k2, k3]
        #c+_1dn c+_2dn c_3dn c_4dn
        hintdn[k1, k4] += 0.5*eqgrdn[k2, k3] * Γ4[k1, k2, k3]
        hintdn[k2, k3] += 0.5*eqgrdn[k1, k4] * Γ4[k1, k2, k3]
        hintdn[k1, k3] -= 0.5*eqgrdn[k2, k4] * Γ4[k1, k2, k3]
        hintdn[k2, k4] -= 0.5*eqgrdn[k1, k3] * Γ4[k1, k2, k3]
    end; end; end
    hintup = hintup / (hsize)^3
    hintdn = hintdn / (hsize)^3
    return hamup + hintup, hamdn + hintdn
end


"""
重新构造平均场哈密顿量
"""
function mean_field_hamiltonian_rl(hamup, hamdn, Γ4, eqgrup, eqgrdn)
    hsize = size(Γ4)[1]
    hintup = zeros(Float64, hsize, hsize)
    hintdn = zeros(Float64, hsize, hsize)
    for cidx in CartesianIndices(Γ4)
        k1, k2, k3, k4 = Tuple(cidx)
        #c+_1up c+_2up c_3up c_4up
        hintup[k1, k4] += 0.5*eqgrup[k2, k3] * Γ4[k1, k2, k3, k4]
        hintup[k2, k3] += 0.5*eqgrup[k1, k4] * Γ4[k1, k2, k3, k4]
        hintup[k1, k3] -= 0.5*eqgrup[k2, k4] * Γ4[k1, k2, k3, k4]
        hintup[k2, k4] -= 0.5*eqgrup[k1, k3] * Γ4[k1, k2, k3, k4]
        #c+_1up c+_2dn c_3dn c_4up
        hintup[k1, k4] += 0.5*eqgrdn[k2, k3] * Γ4[k1, k2, k3, k4]
        hintdn[k2, k3] += 0.5*eqgrup[k1, k4] * Γ4[k1, k2, k3, k4]
        #c+_1dn c+_2up c_3up c_4dn
        hintdn[k1, k4] += 0.5*eqgrup[k2, k3] * Γ4[k1, k2, k3, k4]
        hintup[k2, k3] += 0.5*eqgrdn[k1, k4] * Γ4[k1, k2, k3, k4]
        #c+_1dn c+_2dn c_3dn c_4dn
        hintdn[k1, k4] += 0.5*eqgrdn[k2, k3] * Γ4[k1, k2, k3, k4]
        hintdn[k2, k3] += 0.5*eqgrdn[k1, k4] * Γ4[k1, k2, k3, k4]
        hintdn[k1, k3] -= 0.5*eqgrdn[k2, k4] * Γ4[k1, k2, k3, k4]
        hintdn[k2, k4] -= 0.5*eqgrdn[k1, k3] * Γ4[k1, k2, k3, k4]
    end
    hintup = hintup
    hintdn = hintdn
    return hamup + hintup, hamdn + hintdn
end


"""
创建实空间位置
a1 = (-0.5, 0.5*√3)
a2 = (0.5, 0.5*√3)
"""
function triangle_points(L1, L2)
    a1 = (-0.5, 0.5*√3)
    a2 = (0.5, 0.5*√3)
    latt_arr = Tuple{Float64, Float64}[]
    for yi = 0:1:L2-1; for xi = 0:1:L1-1
        rx = xi*a1[1] + yi*a2[1]
        ry = xi*a1[2] + yi*a2[2]
        push!(latt_arr, (rx, ry))
    end; end
    return latt_arr
end

"""
创建k4tab
b1 = (−2π, 2π/√3)
b2 = (2π, 2π/√3)
"""
function triangle_k4tab(L1, L2)
    b1 = (-2π/L1, 2π/(L1*√3))
    b2 = (2π/L2, 2π/(L2*√3))
    latt_arr = Tuple{Float64, Float64}[]
    for yi = 0:1:L2-1; for xi = 0:1:L1-1
        kx = xi*b1[1] + yi*b2[1]
        ky = xi*b1[2] + yi*b2[2]
        kx, ky = map_to_fbz_triangluar(kx, ky)
        push!(latt_arr, (kx, ky))
    end; end
    Lsq = L1*L2
    k4tab = zeros(Int64, Lsq, Lsq, Lsq)
    for k1=1:1:Lsq; for k2 = 1:1:Lsq; for k3 = 1:1:Lsq
        k4x = latt_arr[k1][1] + latt_arr[k2][1] - latt_arr[k3][1]
        k4y = latt_arr[k1][2] + latt_arr[k2][2] - latt_arr[k3][2]
        #k4x, k4y = map_to_fbz_triangluar(k4x, k4y)
        for pidx = 1:1:length(latt_arr)
            kdiffx = latt_arr[pidx][1] - k4x
            kdiffy = latt_arr[pidx][2] - k4y
            kdiffx, kdiffy = map_to_fbz_triangluar(kdiffx, kdiffy)
            if abs(kdiffx) < 1e-6 && abs(kdiffy) < 1e-6
                k4tab[k1, k2, k3] = pidx
                break
            end
        end
    end; end; end
    return latt_arr, k4tab
end



