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
using KernelAbstractions, StaticArrays, LinearAlgebra
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
♯space = 10
g = god(
    d=sort(SA[zero(T), invϕ, invϕ^2, one(T)]),
    μ=SA[t(), ○, ○, ○],
    ρ=SA[zero(T), zero(T), zero(T), zero(T)],
    ♯=(♯space, ♯space))

Ω.Ο[Ω]
Ω.Ο[Ω.ϵ̃[Ω][1]]
Ω.Ο[ϵ]
t()
t(Ω.ϵ̃[Ω][1])
t(Ω.ϵ̃[Ω.ϵ̃[Ω][1]][1])
Ω.ϵ̃
g.ẑero.μ
g.ône.μ
g = step(g)
g = move(g, SA[t(Ω.ϵ̃[Ω][1]), T(0.75), T(0.75), T(0.75)])
g = move(g, SA[t(Ω.ϵ̃[Ω][1]), T(0.5), T(0.5), T(0.5)])

∃!(g, (x...)->T(0.1))
∃!(g, (x...)->T(0.2))

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
        global buffer = ∃̇(g)
    end
end
const godTASK = @async start(g)

# include("00104_TheoryOfGodTypst.jl")
# Φ_hi = Φ_typst(typst_to_matrix("hi"))
# Φ_hi = Φ_typst("hi")
# gpu_safe(○̂, 2)
# gpu_safe((x,y) -> ○, 2)
# φ_hi, mat_hi = Φ_typst("hi")
# gpu_safe(φ_hi, 2)
# ∃!(g, φ_hi)

# println("idx=2 (tex): ", Array(out)[1])
# # const MAX_RGB = T(mfb_rgb(255, 255, 255))
# # rgb2c(r, g, b) = T(mfb_rgb(r * 255, g * 255, b * 255)) / MAX_RGB
# # c2rgb(c2) = begin
# #     c = floor(UInt32, c2 * MAX_RGB)
# #     ((c >> 16) & 0xFF, (c >> 8) & 0xFF, c & 0xFF) ./ 255
# # end

# # μ(::𝕋) = SVector(ntuple(_ -> ○, length(d)))
# # μ(ϵ::∃) = ϵ.μ
# # ρ(Ω::𝕋) = μ(Ω)
# # ρ(ϵ::∃) = ϵ.ρ



# # include("00090_BroadcastBrowser2Module.jl")
# # import Main.BroadcastBrowserModule: BroadcastBrowser, start
# # include("00105_TheoryOfGodgodBrowser.jl")
# # const BROWSERTASK = Threads.@spawn start(b -> godBrowser(b))
# # g=collect(values(godBROWSER[]))[1].g
