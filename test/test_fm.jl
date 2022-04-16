#=
测试铁磁mf
=#

using LinearAlgebra

"""铁磁的哈密顿量"""
function fmham(grup, grdn)
    hsize = size(grup)
    hamup = zeros(hsize[1], hsize[2])
    hamdn = zeros(hsize[1], hsize[2])
    for i=1:1:hsize[1]; for j=1:1:hsize[1]
        if i == j
            #+ c+_i u c_i u c+_j d c_j d
            #hamup[i, i] += 0.5*grdn[j, j]
            #hamdn[j, j] += 0.5*grup[i, i]
            ##+ c+_i d c_i d c+_j u c_j u
            #hamdn[i, i] += 0.5*grup[j, j]
            #hamup[j, j] += 0.5*grdn[i, i]
            continue
        end
        #- c+_i u c_i u c+_j u c_j u
        #hamup[i, i] += -0.5*grup[j, j]
        #hamup[j, j] += -0.5*grup[i, i]
        #hamup[i, j] += 0.5*grup[j, i]
        #hamup[j, i] += 0.5*grup[i, j]
        ##- c+_i d c_i d c+_j d c_j d
        #hamdn[i, i] += -0.5*grdn[j, j]
        #hamdn[j, j] += -0.5*grdn[i, i]
        #hamdn[i, j] += 0.5*grdn[j, i]
        #hamdn[j, i] += 0.5*grdn[i, j]
        #+ c+_i u c_i u c+_j d c_j d
        hamup[i, i] += 0.5*grdn[j, j]
        hamdn[j, j] += 0.5*grup[i, i]
        #+ c+_i d c_i d c+_j u c_j u
        hamdn[i, i] += 0.5*grup[j, j]
        hamup[j, j] += 0.5*grdn[i, i]
    end; end
    return hamup, hamdn
end


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



slaup = zeros(2, 1)
sladn = zeros(2, 1)


slaup[1, 1] = rand()
slaup[2, 1] = sqrt(1 - slaup[1, 1]^2)
sladn[1, 1] = rand()
sladn[2, 1] = sqrt(1 - sladn[1, 1]^2)
println(slaup)
println(sladn)

grup = slaup * adjoint(slaup)
grdn = sladn * adjoint(sladn)


println(spin_corr(1, 2, grup, grdn))

mfhup, mfhdn = fmham(grup, grdn)

println(mfhup)
println(mfhdn)

slaup = eigvecs(mfhup)[:, 1]
sladn = eigvecs(mfhdn)[:, 1]

println(slaup)
println(sladn)
grup = slaup * adjoint(slaup)
grdn = sladn * adjoint(sladn)
println(spin_corr(1, 2, grup, grdn))

mfhup, mfhdn = fmham(grup, grdn)
slaup = eigvecs(mfhup)[:, 1]
sladn = eigvecs(mfhdn)[:, 1]

println(slaup)
println(sladn)
grup = slaup * adjoint(slaup)
grdn = sladn * adjoint(sladn)
println(spin_corr(1, 2, grup, grdn))


mfhup, mfhdn = fmham(grup, grdn)
slaup = eigvecs(mfhup)[:, 1]
sladn = eigvecs(mfhdn)[:, 1]

println(slaup)
println(sladn)
grup = slaup * adjoint(slaup)
grdn = sladn * adjoint(sladn)
println(spin_corr(1, 2, grup, grdn))


println(slaup)
println(sladn)
grup = slaup * adjoint(slaup)
grdn = sladn * adjoint(sladn)
println(spin_corr(1, 2, grup, grdn))
