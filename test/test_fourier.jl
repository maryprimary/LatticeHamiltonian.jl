#=
测试傅立业变换
=#

using Test
using LatticeHamiltonian
using LinearAlgebra


"""
从x，y获得编号
"""
function get_index_from_posi(L, x, y)
    x_ = x
    while x_ > L || x_ < 1
        x_ = x_ - sign(x_) * L
    end
    y_ = y
    while y_ > L || y_ < 1
        y_ = y_ - sign(y_) * L
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


"""
返回正方格子的hopping
"""
function square_bonds(L)
    bonds = NTuple{3, Int64}[]
    for x = 1:1:L
        for y = 1:1:L
            sidx = get_index_from_posi(L, x, y)
            #
            eidx = get_index_from_posi(L, x+1, y)
            push!(bonds, (sidx, eidx, -1))
            #
            eidx = get_index_from_posi(L, x, y+1)
            push!(bonds, (sidx, eidx, -1))
        end
    end
    return bonds
end


@testset "傅立业" begin
    L = 4
    NWALKERS = 2
    
    pts = Tuple{Float64, Float64}[]
    for yidx = 1:1:L; for xidx = 1:1:L
        push!(pts, (xidx-1, yidx-1))
    end; end

    h0 = zeros(Float64, L^2, L^2)
    bonds = square_bonds(L)
    for bnd in bonds
        h0[bnd[1], bnd[2]] = bnd[3]
        h0[bnd[2], bnd[1]] = bnd[3]
    end
    println(h0)

    kpts, fmat = fourier_matrix(L, (1.0, 0.0), (0.0, 1.0), pts)

    println(kpts)
    println(fmat)
    ftspcmat = rlmat2ftmat(h0, fmat)
    for idx1 = 1:1:L^2; for idx2 = 1:1:L^2
        if idx1 != idx2
        @test abs(ftspcmat[idx1, idx2]) < 1e-10
        end
    end; end
    rlspcmat = ftmat2rlmat(ftspcmat, fmat)
    println("====================")
    println(rlspcmat)
    @test all(isapprox.(h0, rlspcmat, atol=1e-10))
    #
    orbs = eigvecs(rlspcmat)
    #
    for oidx = 1:1:L^2
        orb = orbs[:, oidx:oidx]
        ftvec = rlvec2ftvec(orb, fmat)
        println(ftvec)
        rlvec = ftvec2rlvec(ftvec, fmat)
        @test all(isapprox.(orb, rlvec, atol=1e-10))
        #@test isapprox(sum(abs.(ftvec)), 1/L, rtol=1e-8)
    end
end

