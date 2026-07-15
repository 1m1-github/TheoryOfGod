module TOGMatrix

# using Adapt
using TOGColor: rgba2scalar

const MATRIX_CACHE = Dict{UInt,AbstractMatrix}()

struct ϕMatrix{M}
    matrix::M
    H::UInt32
    W::UInt32
end
function (f::ϕMatrix)(x)
    # @info "(f::ϕMatrix)(x)", x
    px = unsafe_trunc(UInt32, x[2] * (f.W - 1)) + 1 # todo x[2] and x[3], is that general?
    py = unsafe_trunc(UInt32, (1 - x[3]) * (f.H - 1)) + 1
    y = f.matrix[py, px]
    # @info "px,py, f.matrix[py, px]", px,py, y
    y
end
# function Adapt.adapt_structure(to, f::ϕMatrix)
#     ϕMatrix(Adapt.adapt(to, f.mat), f.H, f.W)
# end
function rgbamatrix(matrix)
    H, W = size(matrix)
    # scalar = map(rgba -> rgba2scalar(rgba, one()), matrix)
    scalar = rgba2scalar.(matrix)
    ϕMatrix(scalar, UInt32(H), UInt32(W))
end

end
