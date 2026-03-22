import Pkg
Pkg.activate(".")

"""
TheoryOfGod

I = [ZERO < ○ < ONE] denotes a unit 1-dim space of information with origin ○ (no information) in its center including the corners ZERO and ONE.
Ω = I^I an ∞-dim metric and smooth vector space.
We have a Pretopology 𝕋 on Ω such that ϵᵢ ∈ 𝕋:
* ϵᵢ ⊆ Ω
* ϵ₂ ∈ ϵ₁.ϵ̃ => ϵ₂|ϵ₁ ⊆ ϵ₁ <=> ϵ₂ ⫉ ϵ₁ ⩓ ϵ₂ ∈ ϵ₃.ϵ̃ => ϵ₁ = ϵ₃
* ϵ₁ ≠ ϵ₂ => ϵ₁ ∩ ϵ₂ = ∅
* x ∈ ϵᵢ ⊊ Ω: x.ρ = 0 => ϵᵢ.Φ(x) ∈ I is arbitrary, computable and smooth fuzzy existence potential towards ONE=true xor ZERO=false.

ϵ ⊊ Ω defines its existence inside a subset of Ω using an origin (μ), a radius (ρ) and a closed vs. open in each direction (∂) vector. These vectors are finite and all other dimensional coordinates of ϵ follow from linear interpolation.
If we use a horizontal axis for dimension and a vertical axis for coordinate in the dimension, for any ϵ, the chart looks like a stepwise linear function with finite non-zero radius intervals (active dimensions) and zero interval points within the interpolated regions.
Each child ϵ is a subset of its parent in the active dimensions declared by the parent.

god ⊊ God = Ω = I^I = I^(.) = [ZERO < ○ < ONE]^(.)
god observes or creates, God iterates.
"""
# module tog

# export ∃, ∃̇, ∃!

# using DataStructures
using KernelAbstractions, StaticArrays, LinearAlgebra, MiniFB
using Metal
const GPU_BACKEND = MetalBackend()
const GPU_BACKEND_WORKGROUPSIZE = 2^2^3

const T = Float32

include("00101_TheoryOfGod∃.jl")
const Ω = 𝕋()
# const name = Dict{∃, String}()
include("00103_TheoryOfGodgod.jl")
include("Octahedron.jl")

const invϕ = one(T) / MathConstants.golden
♯space = 10^2
g = god(
    d=sort(SA[zero(T), invϕ, invϕ^2, one(T)]),
    μ=SA[t(), ○, ○, ○],
    ρ=SA[zero(T), zero(T), zero(T), zero(T)],
    ♯=(♯space, ♯space))

include("00102_TheoryOfGodMiniFB.jl")

function start(g::god)
    t = time()
    while true
        # sleep(10)
        yield()
        t̂ = time()
        dt = t̂ - t
        t = t̂
        step(g, dt)
        # p̂ixel = ∃̇(g)
        # @show p̂ixel
        # δ = Δ(pixel, p̂ixel)
        # @show δ
        # @show isempty(δ)
        # isempty(δ) && continue
        # todo
        # global buffer = p̂ixel
        global buffer = floor.(UInt32, reshape(∃̇(g, 10), prod(g.♯)) .* MAX_RGB)
    end
end
# const godTASK = @async start(g)

# include("00104_TheoryOfGodTypst.jl")
# Φ_hi = Φ_typst(typst_to_matrix("hi"))
# Φ_hi = Φ_typst("hi")
# gpu_safe(○̂, 2)
# gpu_safe((x,y) -> ○, 2)
# φ_hi, mat_hi = Φ_typst("hi")
# gpu_safe(φ_hi, 2)
# ∃!(g, φ_hi)

# println("idx=2 (tex): ", Array(out)[1])

# # μ(::𝕋) = SVector(ntuple(_ -> ○, length(d)))
# # μ(ϵ::∃) = ϵ.μ
# # ρ(Ω::𝕋) = μ(Ω)
# # ρ(ϵ::∃) = ϵ.ρ


# # include("00090_BroadcastBrowser2Module.jl")
# # import Main.BroadcastBrowserModule: BroadcastBrowser, start
# # include("00105_TheoryOfGodgodBrowser.jl")
# # const BROWSERTASK = Threads.@spawn start(b -> godBrowser(b))
# # g=collect(values(godBROWSER[]))[1].g



# ∃!(g, (x...)->begin
#     # return value in I=[0,1] given x in Ω=I^I
#     # sqrt(x[2]^2+x[3]^2+x[4]^2) ≤ T(0.1) && return one(T)
#     # zero(T)
#     x[1]
# end)
# ∃!(g, (x...)->T(0.2))
# ∃!(g, (x...)->x[2])
# ∃!(g, x -> x[2])
# ∃!(g, x -> prod(x))
# ∃!(g, x -> begin
#     sqrt(x[2]^2+x[3]^2+x[4]^2) ≤ T(0.5) && return prod(x)
#     zero(T)
# end)
# ∃!(g, x -> begin
#     sum(i -> (x[i] - ○)^2, 2:4) ≤ ○^2 && return prod(x)
#     zero(T)
# end)
# ∃!(g, x -> begin
#     sqrt((x[2]-T(0.5))^2 + (x[3]-T(0.5))^2 + (x[4]-T(0.5))^2) ≤ T(0.5) && return prod(x)
#     zero(T)
# end)
# ∃!(g, x -> begin
#     r2 = (x[2]-T(0.5))^2 + (x[3]-T(0.5))^2 + (x[4]-T(0.5))^2
#     r2 ≤ T(4.0) && return prod(x)
#     zero(T)
# end)
# ∃!(g, x -> begin
#     sin(prod(x))
# end)
# ∃!(g, x -> begin
#     prod(x)
# end)
# ∃!(g, x -> begin
#     sum(x)
# end)
# ∃!(g, x -> begin
#     s = T(10)
#     v = sin(s*x[2])*cos(s*x[3]) + sin(s*x[3])*cos(s*x[4]) + sin(s*x[4])*cos(s*x[2])
#     abs(v) < T(0.3) ? one(T) : zero(T)
# end)
# ∃!(g, x -> (one(T) + sin(T(30) * (x[2] + x[3] + x[4]))) * ○)
# ∃!(g, x -> abs(sin(T(7)*x[2]) * cos(T(11)*x[3]) * sin(T(13)*x[4])))
# ∃!(g, x -> T(50) * abs(sin(T(10)*x[2]) * cos(T(10)*x[3])))
# ∃!(g, x -> begin
#     s = T(5)
#     v = sin(s*x[2])*cos(s*x[3]) + sin(s*x[3])*cos(s*x[4]) + sin(s*x[4])*cos(s*x[2])
#     abs(v) < T(0.5) ? T(100) : zero(T)
# end)
# ∃!(g, x -> T(10))
# ∃!(g, x -> T(100) * abs(sin(T(5)*x[2]) * sin(T(5)*x[3]) * sin(T(5)*x[4])))

# ∃!(g, x -> x[2]) # shows pattern
# ∃!(g, x -> x[2] * x[3]) # shows pattern
# ∃!(g, x -> x[2] * x[2]) # shows pattern

# Fast GPU-safe sin approximation (Bhaskara-style or Taylor)
@inline function gpu_sin(x::T) where T
    # Normalize to [-π, π] range
    x = x - round(x / T(2π)) * T(2π)
    # Taylor series: sin(x) ≈ x - x³/6 + x⁵/120
    x3 = x * x * x
    x5 = x3 * x * x
    x - x3 / T(6) + x5 / T(120)
end
@inline function gpu_cos(x::T) where T
    gpu_sin(x + T(π/2))
end
# ∃!(g, x -> abs(gpu_sin(T(10)*x[2]) * gpu_sin(T(10)*x[3]) * gpu_sin(T(10)*x[4]))) # all black
# ∃!(g, x -> abs(x[2] * T(10) - round(x[2] * T(10))) * abs(x[3] * T(10) - round(x[3] * T(10)))) # works
# ∃!(g, x -> abs(x[2] * x[3] - x[3] * x[4] + x[4] * x[2])) # works

# Sawtooth wave: x - round(x), maps to [-0.5, 0.5]
@inline saw(x) = x - round(x)
# Triangle wave: abs(sawtooth), maps to [0, 0.5]
@inline tri(x) = abs(saw(x))
# Pseudo-sin: triangle wave scaled, roughly approximates sin shape
@inline psin(x) = tri(x * T(0.5)) * T(2)

# # Gyroid approximation
# ∃!(g, x -> begin
#     s = T(3)
#     v = psin(s*x[2])*psin(s*x[3]) + psin(s*x[3])*psin(s*x[4]) + psin(s*x[4])*psin(s*x[2])
#     v > T(0.4) ? one(T) - v : zero(T)
# end) # all black

# # Interference moiré
# ∃!(g, x -> tri(T(8)*x[2]) * tri(T(11)*x[3]) * tri(T(13)*x[4])) # all black

# # Lattice with voids
# ∃!(g, x -> begin
#     a = tri(T(6) * x[2])
#     b = tri(T(6) * x[3])
#     c = tri(T(6) * x[4])
#     min(a, b) * min(b, c)
# end) # all black

# # Diagonal planes intersecting
# ∃!(g, x -> tri(T(5)*(x[2]+x[3])) * tri(T(5)*(x[3]+x[4]))) # works

# # Woven fabric / basket weave
# ∃!(g, x -> begin
#     a = tri(T(8)*x[2] + T(4)*x[4])
#     b = tri(T(8)*x[3] - T(4)*x[4])
#     a * b
# end) # works

# # Turbulent clouds
# ∃!(g, x -> begin
#     a = saw(T(3)*x[2]) * saw(T(5)*x[3])
#     b = saw(T(7)*x[3]) * saw(T(11)*x[4])
#     c = saw(T(13)*x[4]) * saw(T(3)*x[2])
#     abs(a + b + c)
# end) # all black

# Crystal lattice
# ∃!(g, x -> tri(T(6)*x[2]) + tri(T(6)*x[3]) + tri(T(6)*x[4]))

# Egg crate surface
# ∃!(g, x -> tri(T(5)*x[2]) * tri(T(5)*x[3]) + tri(T(5)*x[3]) * tri(T(5)*x[4]))

# Diamond weave
# ∃!(g, x -> tri(T(4)*(x[2]+x[3]+x[4])) * tri(T(4)*(x[2]-x[3])))

# Plaid
# ∃!(g, x -> tri(T(7)*x[2]) + tri(T(7)*x[3]))

# Spiral columns
# ∃!(g, x -> tri(T(6)*x[2] + T(3)*x[4]) + tri(T(6)*x[3] + T(3)*x[4]))

# Warped grid
# ∃!(g, x -> begin
#     a = tri(T(5)*x[2] + x[3]*x[4]*T(10))
#     b = tri(T(5)*x[3] + x[2]*x[4]*T(10))
#     a + b
# end)

# Interference fringes — depth-dependent frequency
# ∃!(g, x -> tri(x[4] * T(20) * x[2]) + tri(x[4] * T(20) * x[3]))

# Nested shells via sum of different frequencies
# ∃!(g, x -> tri(T(3)*(x[2]+x[3])) + tri(T(7)*(x[3]+x[4])) + tri(T(11)*(x[2]+x[4])))

# ∃!(g, x -> begin
#     s = T(1)  # or T(0.001), T(1), etc.
#     px = x[2] * s
#     py = x[3] * s
#     pz = x[4] * s
    
#     iterations = 5
#     scale = one(T)
    
#     for _ in 1:iterations
#         # Fold into positive octant relative to current scale
#         if px + py < zero(T)
#             px, py = -py, -px
#         end
#         if px + pz < zero(T)
#             px, pz = -pz, -px
#         end
#         if py + pz < zero(T)
#             py, pz = -pz, -py
#         end
#         # Scale and translate toward corner
#         px = T(2) * px - scale
#         py = T(2) * py - scale
#         pz = T(2) * pz - scale
#         scale *= T(2)
#     end
    
#     # Distance estimate — inside if close to origin relative to scale
#     d = (abs(px) + abs(py) + abs(pz)) / scale
#     d < T(0.5) ? one(T) : zero(T)
# end)

# ∃!(g, x -> one(T)) # all pink
# ∃!(g, x -> x[2] > zero(T) ? one(T) : zero(T)) # all pink
# ∃!(g, x -> x[2] > one(T) ? one(T) : zero(T)) # boxes with diff colors
# ∃!(g, x -> x[2] > T(10) ? one(T) : zero(T)) # top rect green, bottom rect black
# ∃!(g, x -> x[2] > T(100) ? one(T) : zero(T)) # top rect green, bottom rect black

# Sphere centered at (50, 50, 50) with large radius
# ∃!(g, x -> begin
#     d2 = (x[2]-T(50))*(x[2]-T(50)) + (x[3]-T(50))*(x[3]-T(50)) + (x[4]-T(50))*(x[4]-T(50))
#     d2 < T(2500) ? one(T) : zero(T)
# end) # all black

# Normalize assuming range is roughly [0, 100]
# ∃!(g, x -> begin
#     nx = x[2] * T(0.01)
#     ny = x[3] * T(0.01)
#     nz = x[4] * T(0.01)
#     d2 = (nx - ○)*(nx - ○) + (ny - ○)*(ny - ○) + (nz - ○)*(nz - ○)
#     d2 < T(0.25) ? one(T) : zero(T)
# end) # all black

# (wx, wy) = (0.1715175f0, 0.30113247f0)

# ∃!(g, x -> begin
#     d2 = (x[2]-○)*(x[2]-○) + (x[3]-○)*(x[3]-○) + (x[4]-○)*(x[4]-○)
#     d2 < T(0.1) ? one(T) : zero(T)
# end)

# ∃!(g, x -> begin
#     d2 = (x[2]-one(T))*(x[2]-one(T)) + (x[3]-one(T))*(x[3]-one(T)) + (x[4]-one(T))*(x[4]-one(T))
#     d2 < T(0.09) ? one(T) : zero(T)
# end) # all brownish
# ∃!(g, x -> begin
#     d2 = (x[2]-one(T))*(x[2]-one(T)) + (x[3]-one(T))*(x[3]-one(T)) + (x[4]-one(T))*(x[4]-one(T))
#     d2 < T(0.25) ? one(T) : zero(T)
# end) # all pinkish

# ∃!(g, x -> begin
#     d2 = (x[2]-one(T))*(x[2]-one(T)) + (x[3]-one(T))*(x[3]-one(T)) + (x[4]-one(T))*(x[4]-one(T))
#     d2 < T(0.01) ? one(T) : zero(T)
# end) # moving out enough shows sphere inside another sphere
# ∃!(g, x -> begin
#     d2 = (x[2]-one(T))*(x[2]-one(T)) + (x[3]-one(T))*(x[3]-one(T)) + (x[4]-one(T))*(x[4]-one(T))
#     d2 < T(0.001) ? one(T) : zero(T)
# end) # all black
# ∃!(g, x -> begin
#     d2 = (x[2]-one(T))*(x[2]-one(T)) + (x[3]-one(T))*(x[3]-one(T)) + (x[4]-one(T))*(x[4]-one(T))
#     r2 = T(0.005)
#     d2 < r2 ? one(T) - d2 / r2 : zero(T)
# end) # many concentric circles

# # Torus (ring) — distance from a circle in xz plane
# ∃!(g, x -> begin
#     # Distance from ring of radius R in x-z plane
#     ring = T(0.07)  # major radius
#     tube = T(0.015) # minor radius squared
#     dx = x[2] - one(T)
#     dy = x[3] - one(T)
#     dz = x[4] - one(T)
#     # Distance from ring axis
#     qx = dx * dx + dz * dz
#     # This is rough but GPU safe — avoid sqrt
#     # Approximate: (sqrt(qx) - ring)^2 ≈ qx - 2*ring*sqrt(qx) + ring^2
#     # Instead just use the squared torus implicit:
#     d = (qx + dy*dy + ring*ring - tube) 
#     a = qx
#     # Proper torus: (x²+y²+z²+R²-r²)² - 4R²(x²+z²) = 0
#     v = d*d - T(4)*ring*ring*a
#     v < zero(T) ? T(5) : zero(T)
# end)

# # Two spheres
# ∃!(g, x -> begin
#     ox = T(0.06)
#     d1 = (x[2]-one(T)-ox)*(x[2]-one(T)-ox) + (x[3]-one(T))*(x[3]-one(T)) + (x[4]-one(T))*(x[4]-one(T))
#     d2 = (x[2]-one(T)+ox)*(x[2]-one(T)+ox) + (x[3]-one(T))*(x[3]-one(T)) + (x[4]-one(T))*(x[4]-one(T))
#     r2 = T(0.005)
#     v = zero(T)
#     d1 < r2 && (v += T(5) * (one(T) - d1/r2))
#     d2 < r2 && (v += T(5) * (one(T) - d2/r2))
#     v
# end)

# Sierpinski tetrahedron — now with correct coordinates!
# ∃!(g, x -> begin
#     px = x[2] - one(T)
#     py = x[3] - one(T)
#     pz = x[4] - one(T)
    
#     for _ in 1:5
#         if px + py < zero(T); px, py = -py, -px end
#         if px + pz < zero(T); px, pz = -pz, -px end
#         if py + pz < zero(T); py, pz = -pz, -py end
#         px = T(2) * px - T(0.1)
#         py = T(2) * py - T(0.1)
#         pz = T(2) * pz - T(0.1)
#     end
    
#     d = (abs(px) + abs(py) + abs(pz)) / T(3.2)
#     d < one(T) ? T(5) : zero(T)
# end) 

# Sierpinski — scaled to fit within r=0.1 of center
# ∃!(g, x -> begin
#     s = T(10)  # scale into unit cube
#     px = (x[2] - one(T)) * s + ○
#     py = (x[3] - one(T)) * s + ○
#     pz = (x[4] - one(T)) * s + ○
    
#     for _ in 1:4
#         if px + py < one(T); px, py = one(T) - py, one(T) - px end
#         if px + pz < one(T); px, pz = one(T) - pz, one(T) - px end
#         if py + pz < one(T); py, pz = one(T) - pz, one(T) - py end
#         px = T(2) * px - one(T)
#         py = T(2) * py - one(T)
#         pz = T(2) * pz - one(T)
#     end
    
#     d = abs(px) + abs(py) + abs(pz)
#     d < T(3) ? T(5) : zero(T)
# end)

# Two spheres — bigger, closer together
# ∃!(g, x -> begin
#     ox = T(0.03)
#     r2 = T(0.01)
#     dx1 = x[2] - one(T) - ox
#     dx2 = x[2] - one(T) + ox
#     dy = x[3] - one(T)
#     dz = x[4] - one(T)
#     d1 = dx1*dx1 + dy*dy + dz*dz
#     d2 = dx2*dx2 + dy*dy + dz*dz
#     v = zero(T)
#     d1 < r2 && (v += T(5) * (one(T) - d1/r2))
#     d2 < r2 && (v += T(5) * (one(T) - d2/r2))
#     v
# end)

∃!(g, x -> prod(x))

# const Ω = 𝕋()
# Ω.Ο[Ω]
# Ω.Ο[Ω.ϵ̃[Ω][1]]
# Ω.Ο[ϵ]
# t()
# t(Ω.ϵ̃[Ω][1])
# t(Ω.ϵ̃[Ω.ϵ̃[Ω][1]][1])
# Ω.ϵ̃
g.ẑero.μ
g.ẑero.ρ
g.ône.μ
g.ône.ρ
# g = step(g);
# ∃!(g, (x...)->T(0.1))
# g = move(g, SA[t(Ω.ϵ̃[Ω][1]), T(0.75), T(0.75), T(0.75)])
# ∃!(g, (x...)->T(0.2))
# g = move(g, SA[t(Ω.ϵ̃[Ω][1]), T(0.5), T(0.5), T(0.5)])
