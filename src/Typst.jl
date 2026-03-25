# const DPI = 300
const DPI = 144
const TYPST_TEMPLATE(content) = """
#set page(width: auto, height: auto, margin: (top: 5pt, bottom: 5pt, left: 5pt, right: 5pt))
#set text(size: 10pt)
$content
"""

const TYPST_CACHE = Dict{UInt,SMatrix}()

function typst_to_matrix(typst_code)
    h = hash(typst_code)
    haskey(TYPST_CACHE, h) && return TYPST_CACHE[h]
    cmd = `typst compile - --format png --ppi $DPI -`
    rgba = pipeline(IOBuffer(TYPST_TEMPLATE(typst_code)), cmd) |> read |> IOBuffer |> PNGFiles.load
    TYPST_CACHE[h] = SMatrix{size(rgba)...}(rgb2c.(rgba))
end
function Φ_typst(mat)
    H, W = size(mat)
    x -> begin
        # px = clamp(unsafe_trunc(UInt32, x[2] * (W - 1)) + 1, 1, W)
        # py = clamp(unsafe_trunc(UInt32, x[3] * (H - 1)) + 1, 1, H)
        px = unsafe_trunc(UInt32, x[2] * (W - 1)) + 1
        py = unsafe_trunc(UInt32, x[3] * (H - 1)) + 1
        mat[py, px]
    end
end
typst(code) = Φ_typst(typst_to_matrix(code))
mat = typst_to_matrix("iii")
# φ_hi = Φ_typst(mat)
# gpu_safe(φ_hi, 4)

# @inline function (Φ::ΦTex)(x, y)
#     px = clamp(unsafe_trunc(Int32, x * Φ.W) + Int32(1), Int32(1), Φ.W)
#     py = clamp(unsafe_trunc(Int32, y * Φ.H) + Int32(1), Int32(1), Φ.H)
#     Φ.tex[py, px]
# end
# Φ_typst(code::String, backend=GPU_BACKEND) = 
#     ΦTex(adapt(backend, typst_to_matrix(code))..., Int32.(size(typst_to_matrix(code)))...)
# function Φ_typst(code::String, backend=GPU_BACKEND)
#     mat = typst_to_matrix(code)
#     H, W = size(mat)
#     ΦTex(adapt(backend, mat), Int32(H), Int32(W))
# end

# struct ΦTex
#     tex::Matrix{T}
#     H::Int32
#     W::Int32
# end

# @inline function (Φ::ΦTex)(x, y)
#     px = clamp(unsafe_trunc(Int32, x * Φ.W) + Int32(1), Int32(1), Φ.W)
#     py = clamp(unsafe_trunc(Int32, y * Φ.H) + Int32(1), Int32(1), Φ.H)
#     Φ.tex[py, px]
# end

# # function Φ_typst(code::String)
# #     mat = typst_to_matrix(code)
# #     H, W = size(mat)
# #     ΦTex(mat, Int32(H), Int32(W))
# # end
# function Φ_typst(code::String; offset::Int32=Int32(0))
#     mat = typst_to_matrix(code)
#     H, W = size(mat)
#     ΦTex3(Int32(H), Int32(W), offset), mat
# end