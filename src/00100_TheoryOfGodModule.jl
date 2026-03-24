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

using KernelAbstractions, StaticArrays, LinearAlgebra, MiniFB
using Metal
const GPU_BACKEND = MetalBackend()
const GPU_BACKEND_WORKGROUPSIZE = 2^2^3

const T = Float32

include("00101_TheoryOfGod∃.jl")
const Ω = Ref(𝕋())
# const name = Dict{∃, String}()
include("00103_TheoryOfGodgod.jl")
include("Octahedron.jl")

const invϕ = one(T) / MathConstants.golden
♯space = 10^2
t() = t(Ω[])
const G = Ref(god(
    d=sort(SA[zero(T), invϕ, invϕ^2, one(T)]), # t, x, y, z
    ẑeroμ=SA[t(), ○, ○, ○],
    f̂ocusμ=SA[t(), ○*exp(T(0.1)), ○*exp(T(0.1)), ○*exp(T(0.1))],
    ρ=(T(0.1), T(0.1), one(T)),
    ♯=(♯space, ♯space)))

include("00102_TheoryOfGodMiniFB.jl")

updatebuffer() = begin
    try # todo rm
        Base.invokelatest() do
            global MINIFB_BUFFER[] = floor.(UInt32, reshape(∃̇(G[], Ω[]), prod(G[].♯)) .* MAX_RGB)
        end
    catch e
        showerror(stderr, e, catch_backtrace())
    end
end
const UPDATE_MINIFB_BUFFER_TASK = @async while true
    yield()
    # sleep(1) # todo rm
    updatebuffer()
end

const god_TASK = @async while true
    yield()
    while isready(PENDING_ACTIONS)
        Α = take!(PENDING_ACTIONS)
        global G[] = Α(G[])
    end
end

const TIME_TASK = @async begin
    t = time()
    while true
        yield()
        # sleep(1) # todo rm
        t̃ = time()
        dt = t̃ - t
        t = t̃
        g, δ = step(G[], dt)
        δ || continue
        global G[] = g
    end
end

∃!(G[], x -> prod(x), Ω[])

# g[], δ = step(g[], zero(T))
# ω = Ω[]
# g=G[]
# ∃̇(G[], Ω[])

# const Ω = 𝕋()
# Ω[].Ο[Ω[]]
# Ω.Ο[Ω.ϵ̃[Ω][1]]
# Ω.Ο[ϵ]
# t()
# t(Ω.ϵ̃[Ω][1])
# t(Ω.ϵ̃[Ω.ϵ̃[Ω][1]][1])
# Ω.ϵ̃
# G[].ẑero.μ
# G[].ẑero.ρ
# G[].f̂ocus.μ
# G[].f̂ocus.ρ
# G[].ρ
# G[], δ = step(G[])
# g = G[]
# ∃!(g, (x...)->T(0.1))
# g = move(g, SA[t(Ω.ϵ̃[Ω][1]), T(0.75), T(0.75), T(0.75)])
# ∃!(g, (x...)->T(0.2))
# g = move(g, SA[t(Ω.ϵ̃[Ω][1]), T(0.5), T(0.5), T(0.5)])

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
