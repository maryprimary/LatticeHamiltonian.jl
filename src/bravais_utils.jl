#=
bravais lattice
=#


"""
从x，y获得编号
x, y从0开始
"""
function bravais_index_from_posi(L, x, y)
    x_ = x+1
    while x_ > L || x_ < 1
        x_ = x_ - sign(x_-1) * L
    end
    y_ = y+1
    while y_ > L || y_ < 1
        y_ = y_ - sign(y_-1) * L
    end
    return x_ + (y_ - 1)*L
end


"""
从编号获得xy
"""
function bravais_posi_from_index(L, idx)
    x = mod(idx-1, L)+1
    y = (idx - x) // L
    return x, Int64(y+1)
end



"""
将点对应回第一布里渊区
a1 = (-0.5, 0.5√3)
a2 = (0.5, 0.5√3)
b1 = (−2π, 2π/√3)
b2 = (2π, 2π/√3)
"""
function map_to_fbz_triangluar(ptx, pty)
    destx, desty = ptx, pty
    #用一个长方型的框包起来
    while abs(desty) > π*2/√3
        desty = desty - sign(desty) * π*4/√3
    end
    while abs(destx) > 2*π
        destx = destx - sign(destx) * π*4
    end
    #如果在竖线上，肯定已经在第一布里渊区
    if isapprox(destx, 0.)
        return (destx, desty)
    end
    #如果在横线上
    if isapprox(desty, 0.)
        if destx > pi*4/3
            return (destx-π*2, π*2/√3)
        end
        if desty < -pi*4/3
            return (destx+π*2, π*2/√3)
        end
        return (destx, desty)
    end
    #如果在四个角落上
    lslope = -√3 * sign(destx) * sign(desty)
    y0 = sign(desty)*π*4/√3
    #println(lslope)
    #println(y0)
    #println(lslope*destx + y0)
    if desty*sign(desty) > (lslope*destx + y0)*sign(desty)
        destx = destx - sign(destx)*π*2
        desty = desty - sign(desty)*π*2/√3
    end
    return (destx, desty)
end
