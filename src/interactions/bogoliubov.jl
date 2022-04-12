#=
进行Bogoliubov变换
=#
using LinearAlgebra
#=平均场哈密顿量
H_int = ∑_{σ,k,p} V(k,p) c^+_kσ1 c^+_-kσ2 c_-pσ2 c_pσ1
= ∑_{σ,k,p} V(k,p) ⟨ c^+_kσ1 c^+_-kσ2 ⟩ c_-pσ2 c_pσ1 + V(k,p) c^+_kσ1 c^+_-kσ2 ⟨ c_-pσ2 c_pσ1 ⟩
重新整理左侧的求和,交换kp,而根据对称性的要求，V(k，p)是对称的矩阵
= ∑_{σ,k,p} V(p,k) ⟨ c^+_pσ1 c^+_-pσ2 ⟩ c_-kσ2 c_kσ1 + V(k,p) ⟨ c_-pσ2 c_pσ1 ⟩ c^+_kσ1 c^+_-kσ2
= ∑_{σ,k,p} V(k,p) D*_{σ1σ2}(p) c_-kσ2 c_kσ1 + D_{σ1σ2}(p) c^+_kσ1 c^+_-kσ2
整理指标
= ∑_{σ,k} (∑_{p} V(k,p) D*_{σ1σ2}(p)) c_-kσ2 c_kσ1 + (∑_{p} V(k,p) D_{σ1σ2}(p)) c^+_kσ1 c^+_-kσ2
= ∑_{σ,k} Δ*(k)_{σ1σ2} c_-kσ2 c_kσ1 + Δ(k)_{σ1σ2} c^+_kσ1 c^+_-kσ2
#
#
#Γ(1,2,3,4)= -Γ(2,1,3,4) 如果 <c^+1c^+2> = <c^+2c^+1>
#那么会导致Γ(1,2,3,4)<c^+1c^+2>c_3c_4 和 Γ(1,2,3,4)<c^+2c^+1>c_3c_4相互抵消
#<c^+1c^+2>和<c^+2c^+1>中对称的成分都会抵消
#利用 ⟨ c^+_pσ1 c^+_-pσ2 ⟩ = -⟨ c^+_-pσ2 c^+_pσ1 ⟩
# D*_{σ1σ2}(p) = -D*_{σ2σ1}(-p)
重新整理上面的哈密顿量
= ∑_{σ,k} Δ*(k)_{σ2σ1} c_-kσ1 c_kσ2 + Δ(k)_{σ1σ2} c^+_kσ1 c^+_-kσ2
= ∑_{σ,k} -Δ*(-k)_{σ1σ2} c_-kσ1 c_kσ2 + Δ(k)_{σ1σ2} c^+_kσ1 c^+_-kσ2
#其中Δ(k)_{σ1σ2} = (∑_{p} V(k,p) ⟨ c_-pσ2 c_pσ1 ⟩)
#下面不用这个整理了自旋方向的哈密顿量，直接用一个4×4的
=#

#能带之间是不会混合在一起的
#基(c_k↑ c_k↓ c^+_-k↑ c^+_-k↓)^T
#=
  ϵk↑     0       0      Δk↑↓
   0     ϵk↓     Δk↓↑     0
       -Δ*-k↑↓   -ϵ-k↑    0
-Δ*-k↓↑   0       0      -ϵ-k↓
=#
#=假设知道一个Δ对角化上面的矩阵，得到一个新的本正值和本征失量(ν_ij)
新的本正值为 ξ1 ξ2 ξ3 ξ4 
本征向量 d_j = ∑ ν_{i, j} c_i   c_i = ∑ ν_{i, j} d_j
新的哈密顿量就是 H = ∑_{k} ξ_i d_i，确定基态为ξ<0的态
重新计算 ⟨ c_-pσ2 c_pσ1 ⟩ = ⟨ ∑ ν_{i,j} d_j ∑ ν_{i1,j1} d_j1  ⟩
只有⟨ d_i d_i ⟩可能有数值，⟨ c_-pσ2 c_pσ1 ⟩ = ⟨ ∑_{i, i1, j} ν_{i,j} ν_{i1,j} d_j  d_j ⟩
=#

"""
计算
"""
function bogoliubov_selfconsist(ϵk, ϵ_k, Δkud, Δkdu, Δ_kud, Δ_kdu)
    mat = zeros(Float64, 4, 4)
    mat[1, 1] = ϵk
    mat[2, 2] = ϵk
    mat[3, 3] = -ϵ_k
    mat[4, 4] = -ϵ_k
    mat[1, 4] = Δkud
    mat[2, 3] = Δkdu
    mat[3, 2] = -Δ_kud
    mat[4, 1] = -Δ_kdu
    #只能得出bogoliubov变换后空间中的Slater，没有用
end


"""
构造一个三角格子基态
"""
function triangle_ground_state_construct(N, kpts, kpairs, fmat)
    disps = Tuple{Int64, Float64}[]
    for kp in kpairs
        kx, ky = kpts[kp[1]][1], kpts[kp[1]][2]
        ϵk = -2cos(kx) - 4cos(sqrt(3)*ky/2.0)*cos(kx/2.0)
        push!(disps, (kp[1], ϵk))
        if kp[1] != kp[2]
            nkx, nky = kpts[kp[2]][1], kpts[kp[2]][2]
            ϵ_k = -2cos(nkx) - 4cos(sqrt(3)*nky/2.0)*cos(nkx/2.0)
            push!(disps, (kp[2], ϵ_k))
        end
    end
    sort!(disps, by=x->x[2], rev=true)
    #Slater行列式
    nsite = length(kpts)
    slate = zeros(Float64, nsite, N)
    #能量前几的态
    println(disps, length(disps))
    for eidx = 1:1:N
        (kidx, ϵk) = disps[eidx]
        println(ϵk)
        slate[kidx, eidx] = 1 / sqrt(length(kpts))
    #throw(error("a"))
    end
    println(slate)
    #将动量空间转到实空间
    rlslate = zeros(Float64, nsite, N)
    for eidx = 1:1:N
        println(slate[:, eidx])
        println(ftvec2rlvec(slate[:, eidx], fmat))
        #rlslate[:, eidx] = ftvec2rlvec(slate[:, eidx], fmat)
        #println(rlslate[:, eidx])
    end
    return ftvec2rlvec(slate[:, 1], fmat),ftvec2rlvec(slate[:, 2], fmat)#rlslate[:, 1], rlslate[:, 2]
    #println(rlslate)
end


"""
将三角格子的k值进行配对
"""
function sort_triangluar_kpair(kpts)
    latt_arr = Tuple{Float64, Float64}[]
    for pt in kpts
        kx, ky = map_to_fbz_triangluar(pt[1], pt[2])
        push!(latt_arr, (kx, ky))
    end
    #latt_arr中和kpts中是一样的
    #格点的配对数组
    pair_arr = zeros(Int64, length(latt_arr))
    for kidx = 1:1:length(latt_arr)
        nkv = map_to_fbz_triangluar(-latt_arr[kidx][1], -latt_arr[kidx][2])
        tpidx = -1
        for pidx = 1:1:length(latt_arr)
            if all(isapprox.(latt_arr[pidx], nkv, atol=1e-6))
                tpidx = pidx
                break
            end
        end
        if tpidx == -1
            tpidx = kidx
        end
        #println(latt_arr[kidx], __map_to_fbz_triangluar(-latt_arr[kidx][1], -latt_arr[kidx][2]))
        pair_arr[kidx] = tpidx
        #println(kidx, " ",tpidx)
        if pair_arr[tpidx] != 0 && pair_arr[tpidx] != kidx
            throw(error("duplicate"))
        end
        pair_arr[tpidx] = kidx
    end
    kpairs = NTuple{2, Int64}[]
    for pidx = 1:1:length(pair_arr)
        if pair_arr[pidx] != -1
            push!(kpairs, (pidx, pair_arr[pidx]))
            pair_arr[pair_arr[pidx]] = -1
        end
    end
    return kpairs
end



"""
创建成对的k
b1 = (−2π, 2π/√3)
b2 = (2π, 2π/√3)
"""
function create_triangluar_kpair(L)
    #在惯用晶胞中设置k点
    b1 = (-2π/L, 2π/(L*√3))
    b2 = (2π/L, 2π/(L*√3))
    #给动量空间的格点一个编号
    latt_arr = Tuple{Float64, Float64}[]
    for xi = 0:1:L-1; for yi = 0:1:L-1
        kx = xi*b1[1] + yi*b2[1]
        ky = xi*b1[2] + yi*b2[2]
        #kx, ky = __map_to_fbz_triangluar(kx, ky)
        push!(latt_arr, (kx, ky))
    end; end
    println(latt_arr)
    return latt_arr, sort_triangluar_kpair(latt_arr)
end


