#=
测试hf
=#

using LatticeHamiltonian
using Plots
using LinearAlgebra
using Test
using JLD

"""
平均场哈密顿量下的自旋关联
"""
function spin_corr(i, j, grup, grdn)
    if i == j
        throw(error("not implement"))
    end
    corr = 0.
    #n_i↑ n_j↑ / c+_i c_i c+_j c_j
    corr += grup[i, i] * grup[j, j]
    corr += grup[i, j] * (-grup[j, i])
    #-n_i↑ n_j↓
    corr -= grup[i, i] * grdn[j, j]
    #-n_i↓ n_j↑
    corr -= grdn[i, i] * grup[j, j]
    #n_i↓ n_j↓
    corr += grdn[i, i] * grdn[j, j]
    corr += grdn[i, j] * (-grdn[j, i])
    return corr
end



L = 6
N = 19

rlpts = Tuple{Float64, Float64}[]
for yidx = 1:1:L; for xidx = 1:1:L
    push!(rlpts, (xidx-1, yidx-1))
end; end

kpts, k4tab = triangle_k4tab(L)

ftpts, fmat = fourier_matrix(L, (-0.5, 0.5*sqrt(3)), (0.5, 0.5*sqrt(3)), rlpts)
#ftop = rlmat2ftmat(h0, fmat)

fmat2 = fourier_matrix_pts(kpts, (-0.5, 0.5*sqrt(3)), (0.5, 0.5*sqrt(3)), rlpts)


#for kp in zip(kpts, ftpts)
#    println(kp[1], " ", kp[2])
#end

@test all(isapprox.(fmat, fmat2, atol=1e-10))


#latplt = plot()
#scatter!(latplt, kpts)
#savefig(latplt, "lat.png")
#
##println(k4tab)
#
#Γ4 = zeros(Float64, L^2, L^2, L^2)
##
#for k1=1:1:L^2; for k2=1:1:L^2; for k3=1:1:L^2
#    global Γ4
#    if k1 == k3
#        Γ4[k1, k2, k3] += 2
#    end
#    if k2 == k3
#        Γ4[k1, k2, k3] += 1
#    end
#end; end; end
#
#htp = heatmap(Γ4[:, :, 1])
#
#savefig(htp, "htp.png")


h0 = zeros(Float64, L^2, L^2)
for (kid, kpt) in enumerate(kpts)
    kx, ky = kpt[1], kpt[2]
    ϵ = -2cos(kx) - 4cos(sqrt(3)*ky/2.0)*cos(kx/2.0)
    h0[kid, kid] = ϵ / L^2
    #println(ϵ)
end

h0up = copy(h0)
h0dn = copy(h0)

slatup, slatdn = mean_field_slater(h0up, h0dn, N)
eqgrtup, eqgrtdn = mean_field_parameters(slatup, slatdn)


hamup_rl = ftmat2rlmat(h0up, fmat2)
hamdn_rl = ftmat2rlmat(h0dn, fmat2)

slaup_rl, sladn_rl = mean_field_slater(hamup_rl, hamdn_rl, N)
grup_rl, grdn_rl = mean_field_parameters(slaup_rl, sladn_rl)

println(spin_corr(1, 2, grup_rl, grdn_rl))


#h0_rl = ftmat2rlmat(h0, fmat2)
#println(h0_rl - h0_rl')
#
#slat_rl = mean_field_slater(h0_rl, N)
#eqgrt_rl = mean_field_parameters(slat_rl)
#
#sc = spin_corr(1, 2, eqgrt_rl)
#println(sc)

Γ4 = zeros(Float64, L^2, L^2, L^2)
for k1=1:1:L^2; for k2=1:1:L^2; for k3=1:1:L^2
    global Γ4
    if k1 == k3
        Γ4[k1, k2, k3] += 20
    end
    if k2 == k3
        Γ4[k1, k2, k3] += 10
    end
end; end; end

#println(k4tab[2, 1, 2])

hamup, hamdn = mean_field_hamiltonian(h0up, h0dn, Γ4, k4tab, eqgrtup, eqgrtdn)

#println(hamup)

#exit()

hamup_rl = ftmat2rlmat(hamup, fmat2)
hamdn_rl = ftmat2rlmat(hamup, fmat2)

slaup_rl, sladn_rl = mean_field_slater(hamup_rl, hamdn_rl, N)
grup_rl, grdn_rl = mean_field_parameters(slaup_rl, sladn_rl)

println(spin_corr(1, 2, grup_rl, grdn_rl))
#println(hamup)

for iidx = 1:1:10
    global hamup, hamdn
    slaup, sladn = mean_field_slater(hamup, hamdn, N)
    eqgrup, eqgrdn = mean_field_parameters(slaup, sladn)
    hamup, hamdn = mean_field_hamiltonian(h0up, h0dn, Γ4, k4tab, eqgrup, eqgrdn)
end

#println(hamup)


hamup_rl = ftmat2rlmat(hamup, fmat2)
hamdn_rl = ftmat2rlmat(hamup, fmat2)
slaup_rl, sladn_rl = mean_field_slater(hamup_rl, hamdn_rl, N)
grup_rl, grdn_rl = mean_field_parameters(slaup_rl, sladn_rl)

println(spin_corr(1, 2, grup_rl, grdn_rl))


#实空间

#Γ4 = Γ4 / L^2 / L^2 / L^2
#Γ4rl = real(ftint2rlint(Γ4, fmat2, k4tab))
#save("Gamma4.jld", "Gamma4rl", Γ4rl)

fjld = load("Gamma4.jld")
Γ4rl = fjld["Gamma4rl"]
#zeros(Float64, L^2, L^2, L^2, L^2)#fjld["Gamma4rl"]
#Γ4rl[1, 2, 2, 1] = 1.
#Γ4rl[1, 2, 1, 2] = 2.

println(Γ4rl[1, 2, 2, 1])
println(Γ4rl[1, 2, 1, 2])

#exit()

#println(Γ4rl[2, 1, 2, 1])

#exit()
#Γ4 = zeros(Float64, L^2, L^2, L^2, L^2)
#for k1=1:1:L^2
#    global Γ4
#    Γ4[k1, k1, k1, k1] += 1
#end



h0up_rl = ftmat2rlmat(h0up, fmat2)
h0dn_rl = ftmat2rlmat(h0dn, fmat2)

#println((Int64 ∘ round).(h0up_rl))

slaup_rl, sladn_rl = mean_field_slater(h0up_rl, h0dn_rl, N)
grup_rl, grdn_rl = mean_field_parameters(slaup_rl, sladn_rl)

println(spin_corr(1, 3, grup_rl, grdn_rl))

hamup_rl, hamdn_rl = mean_field_hamiltonian_rl(h0up_rl, h0dn_rl, Γ4rl, grup_rl, grdn_rl)

println(hamdn_rl[1:6, 1:6])

slaup_rl, sladn_rl = mean_field_slater(hamup_rl, hamdn_rl, N)
grup_rl, grdn_rl = mean_field_parameters(slaup_rl, sladn_rl)
println(spin_corr(1, 3, grup_rl, grdn_rl))

#exit()


for iidx = 1:1:10
    global hamup_rl, hamdn_rl, slaup_rl, sladn_rl, grup_rl, grdn_rl
    slaup_rl, sladn_rl = mean_field_slater(hamup_rl, hamdn_rl, N)
    grup_rl, grdn_rl = mean_field_parameters(slaup_rl, sladn_rl)
    println(spin_corr(1, 3, grup_rl, grdn_rl))
    hamup_rl, hamdn_rl = mean_field_hamiltonian_rl(h0up_rl, h0dn_rl, Γ4rl, grup_rl, grdn_rl)
end
slaup_rl, sladn_rl = mean_field_slater(hamup_rl, hamdn_rl, N)
grup_rl, grdn_rl = mean_field_parameters(slaup_rl, sladn_rl)
println(spin_corr(1, 3, grup_rl, grdn_rl))

