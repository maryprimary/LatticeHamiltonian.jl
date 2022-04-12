#=
实空间的两粒子
=#


struct FourFermionInteraction
    name :: String
    type :: Symbol
    V :: Array{Float64, 4}
end


struct TwoOrderChannel
    name :: String
    type :: Symbol
    V :: Matrix{Float64}
end



"""
将结果打印成容易观察
"""
function pretty(ffi::FourFermionInteraction) :: Union{Missing, String}
    missing
end


"""
将结果打印成容易观察的格式
"""
function pretty(toc::TwoOrderChannel) :: String
    str = name
end



