#=
测试kpair
=#

using LatticeHamiltonian
using Plots
using LinearAlgebra

L = 6

function brillouin()
    edgelen = 4pi/3.0
    xlim = edgelen / 2.0
    ylim = sqrt(3) * xlim
    vxs = [2xlim, xlim, -xlim, -2xlim, -xlim, +xlim, 2xlim]
    vys = [0, ylim, ylim, 0, -ylim, -ylim, 0]
    plt = plot(vxs, vys, label="")
    return plt
end

"""
从x，y获得编号
"""
function get_index_from_posi(L, x, y)
    x_ = x
    while x_ > L || x_ < 1
        x_ = x_ - sign(x_-1) * L
    end
    y_ = y
    while y_ > L || y_ < 1
        y_ = y_ - sign(y_-1) * L
    end
    return x_ + (y_ - 1)*L
end


"""
从编号获得xy
"""
function get_posi_from_index(L, idx)
    x = mod(idx, L)
    y = (idx - x) // L
    return x, y
end


function triangle_bonds(L)
    bonds = NTuple{3, Int64}[]
    for x = 1:1:L
        for y = 1:1:L
            sidx = get_index_from_posi(L, x, y)
            #
            eidx = get_index_from_posi(L, x, y+1)
            push!(bonds, (sidx, eidx, -1))
            #
            eidx = get_index_from_posi(L, x+1, y-1)
            push!(bonds, (sidx, eidx, -1))
            #
            eidx = get_index_from_posi(L, x-1, y)
            push!(bonds, (sidx, eidx, -1))
        end
    end
    return bonds
end



h0 = zeros(Float64, L^2, L^2)
bonds = triangle_bonds(L)
for bnd in bonds
    h0[bnd[1], bnd[2]] = bnd[3]
    h0[bnd[2], bnd[1]] = bnd[3]
end

latplt = plot()
for bnd in bonds
    global latplt
    x1, y1 = get_posi_from_index(L, bnd[1])
    x2, y2 = get_posi_from_index(L, bnd[2])
    plot!(latplt, [x1, x2], [y1, y2], label="")
end
savefig(latplt, "lat.png")
#
println(eigvals(h0))
eigve = eigvecs(h0)
n1 = eigve[:, 36]
n2 = eigve[:, 35]
#println(eigve[:, 36])
#println(eigve[:, 35])
#
#
#plt = brillouin()
#
rlpts = Tuple{Float64, Float64}[]
for yidx = 1:1:L; for xidx = 1:1:L
    push!(rlpts, (xidx-1, yidx-1))
end; end

ftpts, fmat = fourier_matrix(L, (-0.5, 0.5*sqrt(3)), (0.5, 0.5*sqrt(3)), rlpts)
ftop = rlmat2ftmat(h0, fmat)

for idx = 1:1:36
    println(ftop[idx, idx])
end
#
println(ftpts[17], ftpts[32])
#println(fmat[17, :], fmat[31, :])
#vec1 = zeros(36)
#vec1[17] = 1/6
#rvec1 = ftvec2rlvec(vec1, fmat)
#println(rvec1)
#println("========")
#println(h0*rvec1)
#println("========")
#println(h0*rvec1 ./ rvec1)
#println("========")
#
#vec1 = zeros(36)
#vec1[17] = 1/6
#rvec1 = ftvec2rlvec(vec1, fmat)
#println(norm(rvec1))
#println(rvec1)
#println("========")
#println(h0*rvec1)
#println("========")
#println(h0*rvec1 ./ rvec1)
#println("========")
#
##println(fmat)
##println(ftpts)
tkp = sort_triangluar_kpair(ftpts)
##println(tkp)
#tidx = 1
#for kpi in tkp
#    scatter!(plt,
#    [ftpts[kpi[1]][1], ftpts[kpi[2]][1]],
#    [ftpts[kpi[1]][2], ftpts[kpi[2]][2]],
#    text=string(tidx), label="")
#    global tidx += 1
#end
#savefig(plt, "test_kpair.png")
#
#
#
g1, g2 = triangle_ground_state_construct(2, ftpts, tkp, fmat)
#
#
c1 = dot(g1, n1)
c2 = dot(g2, n1)
println(c1, c2)
println(c1 * g1 + c2 * g2)
println("=======")
println(n1)
##
##println(g2)
##println("=========")
##println(h0*g2)
##println("=========")
##println(h0*g2 ./ g2)
#
#ftop = rlmat2ftmat(h0, fmat)
#
#for idx = 1:1:36
#    println(ftop[idx, idx])
#end
#println(ftop[17, 17])
#
#println(ftpts[32])
##println(tkp[32])
