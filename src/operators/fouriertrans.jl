#=
傅立叶变换的矩阵
=#


"""
傅立叶变换的矩阵
最后的pts按照编号的顺序给出点的坐标
"""
function fourier_matrix(L, a1, a2, pts)
    #计算b1,b2
    area = a1[1] * a2[2] - a1[2] * a2[1]
    b1 = (a2[2]/area, -a2[1]/area)
    b2 = (-a1[2]/area, a1[1]/area)
    println(b1)
    println(b2)
    #
    #println(b1)
    #println(b2)
    #println(a1[1]*b1[1] + a1[2]*b1[2])
    #println(a1[1]*b2[1] + a1[2]*b2[2])
    #println(a2[1]*b1[1] + a2[2]*b1[2])
    #println(a2[1]*b2[1] + a2[2]*b2[2])
    ftspc = Tuple{Float64, Float64}[]
    for yidx = 1:1:L; for xidx = 1:1:L
        push!(ftspc, (
            2π*((xidx-1)*b1[1]+(yidx-1)*b2[1])/L,
            2π*((xidx-1)*b1[2]+(yidx-1)*b2[2])/L
        ))
    end; end
    rlspc = Tuple{Float64, Float64}[]
    for pidx = 1:1:length(pts)
        pt = pts[pidx]
        push!(rlspc, (pt[1]*a1[1]+pt[2]*a2[1], pt[1]*a1[2]+pt[2]*a2[2]))
    end
    tmat = Matrix{ComplexF64}(undef, length(ftspc), length(rlspc))
    for kidx = 1:1:length(ftspc); for ridx = 1:1:length(rlspc)
        kx = ftspc[kidx][1]; ky = ftspc[kidx][2]
        rx = rlspc[ridx][1]; ry = rlspc[ridx][2]
        tmat[kidx, ridx] = exp(-im*(kx*rx+ky*ry))
    end; end
    return ftspc, tmat
end



"""
傅立叶变换的矩阵
最后的pts按照编号的顺序给出点的坐标
"""
function fourier_matrix_pts(ftpts, a1, a2, rlpts)
    ftspc = ftpts
    rlspc = Tuple{Float64, Float64}[]
    for pidx = 1:1:length(rlpts)
        pt = rlpts[pidx]
        push!(rlspc, (pt[1]*a1[1]+pt[2]*a2[1], pt[1]*a1[2]+pt[2]*a2[2]))
    end
    tmat = Matrix{ComplexF64}(undef, length(ftspc), length(rlspc))
    for kidx = 1:1:length(ftspc); for ridx = 1:1:length(rlspc)
        kx = ftspc[kidx][1]; ky = ftspc[kidx][2]
        rx = rlspc[ridx][1]; ry = rlspc[ridx][2]
        tmat[kidx, ridx] = exp(-im*(kx*rx+ky*ry))
    end; end
    return tmat
end



"""
将实空间的算符转移到动量空间
c_r = 1/n ∑_{k} adj(tmat_{k, r}) c_k
c^+_r = 1/n ∑_{k} tmat_{k, r} c^+_k
"""
function rlmat2ftmat(rlmat, tmat)
    rlsize = size(rlmat)
    ftmatc = zeros(ComplexF64, rlsize[1], rlsize[2])
    for rc in 1:1:rlsize[1]; for ra in 1:1:rlsize[2]
        rcoef = rlmat[rc, ra]
        for fc in 1:1:rlsize[1]; for fa in 1:1:rlsize[2]
            ftmatc[fc, fa] += rcoef * tmat[fc, rc] * adjoint(tmat[fa, ra])
        end; end
    end; end
    #注意哈密顿量就是n*n的，n就是格点数
    ftmatc = ftmatc / rlsize[1] / rlsize[1]
    return real(ftmatc)
end


"""
将动量空间的算符转移到实空间
c_k = ∑_{r} tmat_{k, r} c_r
c^+_k = ∑_{r} adj(tmat_{k, r}) c^+_r
"""
function ftmat2rlmat(ftmat, tmat)
    ftsize = size(ftmat)
    rlmatc = zeros(ComplexF64, ftsize[1], ftsize[2])
    for fc = 1:1:ftsize[1]; for fa = 1:1:ftsize[2]
        fcoef = ftmat[fc, fa]
        for rc = 1:1:ftsize[1]; for ra = 1:1:ftsize[2]
            rlmatc[rc, ra] += fcoef * adjoint(tmat[fc, rc]) * tmat[fa, ra]
        end; end
    end; end
    return real(rlmatc)
end


"""
将实空间的向量转移到动量空间
c_r = 1/n ∑_{k} adj(tmat_{k, r}) c_k
c^+_r = 1/n ∑_{k} tmat_{k, r} c^+_k
"""
function rlvec2ftvec(rlvec, tmat)
    rlsize = size(tmat)
    ftv = tmat * rlvec
    return ftv / rlsize[1]
end


"""
将动量空间的向量转移到实空间
c_k = ∑_{r} tmat_{k, r} c_r
c^+_k = ∑_{r} adj(tmat_{k, r}) c^+_r
"""
function ftvec2rlvec(ftvec, tmat)
    nsite = length(ftvec)
    rlv = zeros(ComplexF64, nsite)
    for nidx = 1:1:nsite
        fcoef = ftvec[nidx]
        for idx = 1:1:nsite
            rlv[idx] += fcoef * adjoint(tmat[nidx, idx])
        end
    end
    return rlv
end



"""
将相互作用从动量空间变换
c_k = ∑_{r} tmat_{k, r} c_r
c^+_k = ∑_{r} adj(tmat_{k, r}) c^+_r
"""
function ftint2rlint(ftint, tmat, k4tab)
    nsite = size(ftint)[1]
    rlint = zeros(ComplexF64, nsite, nsite, nsite, nsite)
    for k1=1:1:nsite; for k2=1:1:nsite; for k3=1:1:nsite
        k4 = k4tab[k1, k2, k3]
        fcoef = ftint[k1, k2, k3]
        for i1=1:1:nsite; for i2=1:1:nsite; for i3=1:1:nsite; for i4=1:1:nsite
            rlint[i1, i2, i3, i4] += fcoef * adjoint(tmat[k1, i1]) *
            adjoint(tmat[k2, i2]) * tmat[k3, i3] * tmat[k4, i4]
        end; end; end; end
    end; end; end
    return rlint
end


"""
将相互作用从实空间变换
c_r = 1/n ∑_{k} adj(tmat_{k, r}) c_k
c^+_r = 1/n ∑_{k} tmat_{k, r} c^+_k
"""
function rlint2ftint(rlint, tmat, k4tab)
    nsite = size(rlint)[1]
    ftint = zeros(ComplexF64, nsite, nsite, nsite)
    for i1=1:1:nsite; for i2=1:1:nsite; for i3=1:1:nsite; for i4=1:1:nsite
        rcoef = rlint[i1, i2, i3, i4]
        for k1=1:1:nsite; for k2=1:1:nsite; for k3=1:1:nsite
            k4 = k4tab[k1, k2, k3]
            ftint[k1, k2, k3] += rcoef * tmat[k1, i1] *
            tmat[k2, i2] * adjoint(tmat[k3, i3]) * adjoint(tmat[k4, i4])
        end; end; end
    end; end; end; end
    return ftint / (nsite^4)
end

