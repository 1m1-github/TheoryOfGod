struct god
    ẑero::∃
    ône::∃
    ∂t₀::Bool
    v::T
    ρ::T
    Ω::𝕋
    ⚷::UInt
    ♯::NTuple
    ∇::UInt
    d::Function
end
function god(; d, μ, ρ, ⚷=zero(UInt), Φ=○̂, ♯=SA[1], ∇=typemax(UInt))
    N = length(d)
    ∂₀ = SVector(ntuple(_ -> (true, false), N))
    ẑero = ∃(Ω, d, μ, ρ, ∂₀, Φ)
    # μ₁ = SA[μ[1], ones(T, N - 1)...]
    ones = @SVector ones(T, N)
    ∂₁ = SVector(ntuple(_ -> (false, true), N))
    zeros = @SVector zeros(T, N)
    ône = ∃(Ω, d, ones, zeros, ∂₁, ○̂)
    god(ẑero, ône, true, zero(T), zero(T), 𝕋(), ⚷, ♯, ∇, (x, y) -> sqrt.(x .^ 2 .+ y .^ 2))
end
dh(⚷, g, n) = powermod(g, ⚷, n)
⚷⚷(⚷, dh) = powermod(dh, ⚷, DH_N)
⚷i(i, ⚷, ♯) =
    ntuple(length(♯)) do d
        mod1(i[d] + ⚷, ♯[d])
    end
i⚷(i, ⚷, ♯) =
    ntuple(length(♯)) do d
        î = mod(⚷, ♯[d])
        mod1(i[d] + ♯[d] - î, ♯[d])
    end
function ∃!(g::god, Φ, ω=Ω)
    ϵ = -(g.ône, g.ẑero, ω)
    ṫ = t(ω.Ο[ω] + 1)
    ρt = (one(T) - ṫ) * ○
    μt = ṫ + ρt
    μ = SA[μt, ϵ.μ[2:end]...]
    ρ = SA[ρt, ϵ.ρ[2:end]...]
    ϵ = ∃(ϵ, ϵ.d, μ, ρ, ϵ.∂, Φ)
    ∃!(ϵ, ω)
end

trivial(ϵ) = ϵ isa 𝕋 || ϵ.Φ === ○̂
function ∃̇(g::god, ∇̄=1, ω=Ω)
    ϵ = β(-(g.ône, g.ẑero, ω), ω)
    # todo dim==2
    ẽx, ẽy, wx, wy = calc_ew(g.ẑero.μ[end-2:end], g.ône.μ[end-2:end])
    istrivial = trivial(ϵ)
    N = length(g.ẑero.μ)
    ex = SVector(zeros(T, N - 3)..., ẽx...)
    ey = SVector(zeros(T, N - 3)..., ẽy...)
    c = (g.ẑero.μ .+ g.ône.μ) * ○
    d = g.ône.μ .- g.ẑero.μ
    @show c d ẽx ẽy wx wy
    i = fill(istrivial ? 0 : 1, g.♯..., GL_N)
    Φ̃Φ̃ = []
    !istrivial && push!(Φ̃Φ̃, ϵ.Φ)
    owners!(g, ϵ, i, Φ̃Φ̃, 0, ∇̄, ex, ey, wx, wy, ω)
    ΦΦ = ΦTuple(ntuple(i -> Φ̃Φ̃[i], length(Φ̃Φ̃)))
    project(g, ΦΦ, i, c, d, ex, ey, wx, wy)
end
function owners!(g, ϵ, i, ΦΦ, ∇, ∇̄, ex, ey, wx, wy, ω)
    if 0 < ∇ && ϵ isa ∃ && !trivial(ϵ)
        intersects = pyramid_box_intersection(
            i, length(ΦΦ) + 1,
            g.ẑero.μ, g.ône.μ,
            ex, ey, wx, wy,
            ϵ.μ .- ϵ.ρ, ϵ.μ .+ ϵ.ρ,
            g.♯..., GL_N)
        intersects || return
        push!(ΦΦ, ϵ.Φ)
    end
    ∇ == ∇̄ && return
    for ϵ̃ = ω.ϵ̃[ϵ]
        owners!(g, ϵ̃, i, ΦΦ, ∇ + 1, ∇̄, ex, ey, wx, wy, ω)
    end
end

struct ΦTuple{ΦT}
    ϕ::ΦT
end
@generated function Φ̇(Φ::ΦTuple{ΦT}, i, x) where ΦT
    N = length(ΦT.parameters)
    branches = []
    for ĩ = 1:N
        push!(branches, quote
            if i == $ĩ
                return Φ.ϕ[$ĩ](x)
            end
        end)
    end
    quote
        $(branches...)
        return zero(T)
    end
end
function project(g, ΦΦ, i, c, d, ex, ey, wx, wy)
    out = KernelAbstractions.zeros(GPU_BACKEND, T, g.♯...)
    i̇ = KernelAbstractions.allocate(GPU_BACKEND, UInt16, size(i))
    copyto!(i̇, i)
    Base.invokelatest() do
        project!(GPU_BACKEND, GPU_BACKEND_WORKGROUPSIZE)(
            out,
            ΦΦ, i̇, c, d, ex, ey, wx, wy, g.♯...,
            ndrange=g.♯
        )
    end
    KernelAbstractions.synchronize(GPU_BACKEND)
    Array(out)
end
@kernel function project!(out, ΦΦ, i, c, d, ex, ey, wx, wy, nx, ny)
    o = zero(T)
    ĩ = @index(Global, NTuple)
    for zi = 1:GL_N
        ϕi = i[ĩ..., zi]
        if iszero(ϕi)
            o += ○ * GL_WEIGHTS[zi]
            continue
        end
        t = zi / (GL_N + one(T))
        omt = one(T) - t
        wxk, wyk = wx * omt, wy * omt
        si = (2 * ĩ[1] - one(T) - nx) / (nx - one(T))
        sj = (2 * ĩ[2] - one(T) - ny) / (ny - one(T))
        ck = c .+ t * d
        x = ck .+ si * wxk * ex .+ sj * wyk * ey
        o += Φ̇(ΦΦ, ϕi, x) * GL_WEIGHTS[zi]
    end
    out[ĩ...] = one(T) - exp(-o)
end
function step(g::god, dt̂=one(T))
    if g.∂t₀
        ṫ = t()
        μ = SVector(ntuple(i -> i == 1 ? ṫ : g.ẑero.μ[i], length(g.ẑero.μ)))
    else
        δμ = g.ône.μ .- g.ẑero.μ
        all(d -> iszero(d), δμ) && return
        α = clamp(g.v * dt̂, zero(T), one(T))
        μ = g.ẑero.μ .+ α .* δμ
    end
    g = move(g, μ)
    # g.ẑero.Φ !== ○̂ && ∃!(g.ẑero)
    g
end
jerk(g::god, δ) = accelerate(g, g.v * exp(δ))
accelerate(g::god, δ) = speed(g, iszero(g.v) ? δ : g.v * exp(δ))
speed(g::god, v) = god(g.ẑero, g.ône, g.∂t₀, clamp(T(v), zero(T), one(T)), g.ρ, g.Ω, g.⚷, g.♯, g.∇, g.d)
# stop(g::god) = speed(g, zero(T))
# stoptime(g::god) = god(g.ẑero, g.ône, g.∂t₀, SA[zero(T), g.v[2:end]...], g.ρ, g.Ω, g.⚷, g.♯, g.∇, g.d)
scale(g::god, δ, ω=Ω) = begin
    ϵ = -(g.ône, g.ẑero, ω)
    ρ = iszero(ϵ.ρ) ? δ : min.(ϵ.ρ * exp(δ), ○)
    N = length(ϵ.μ)
    ône = SVector(min.(ϵ.μ .+ ρ, ones(T, N)))
    ẑero = SVector(max.(zeros(T, N), ϵ.μ .- ρ))
    move(g, ône) # todo could be parallel
    move(g, ẑero) # todo could be parallel
end
move(g::god, ẑeroμ) =
    god(
        ∃(g.ẑero.ϵ̂, g.ẑero.d, ẑeroμ, g.ẑero.ρ, g.ẑero.∂, g.ẑero.Φ),
        # ∃(g.ône.ϵ̂, g.ône.d, SA[ẑeroμ[1], g.ône.μ[2:end]...], g.ône.ρ, g.ône.∂, g.ône.Φ),
        g.ône,
        g.∂t₀, g.v, g.ρ, g.Ω, g.⚷, g.♯, g.∇, g.d
    )
focus(g::god, ôneμ) =
    god(
        # ∃(g.ẑero.ϵ̂, g.ẑero.d, SA[ôneμ[1], g.ẑero.μ[2:end]...], g.ẑero.ρ, g.ẑero.∂, g.ẑero.Φ),
        g.ẑero,
        ∃(g.ône.ϵ̂, g.ône.d, ôneμ, g.ône.ρ, g.ône.∂, g.ône.Φ),
        g.∂t₀, g.v, g.ρ, g.Ω, g.⚷, g.♯, g.∇, g.d
    )
move(g::god, d, δ) = move(g, SVector(ntuple(i -> begin
        g.ẑero.d[i] == d && return δ(g.ẑero.μ[i], T(0.01))
        g.ẑero.μ[i]
    end, length(g.ẑero.μ))))
focus(g::god, d, δ) = focus(g, SVector(ntuple(i -> begin
        g.ône.d[i] == d && return δ(g.ône.μ[i], T(0.01))
        g.ône.μ[i]
    end, length(g.ône.μ))))
focusup(g, d) = focus(g, d, +)
focusdown(g, d) = focus(g, d, -)
moveup(g, d) = move(g, d, +)
movedown(g, d) = move(g, d, -)
jerkdown(g) = jerk(g, T(-0.01))
jerkup(g) = jerk(g, T(0.01))
scaledown(g) = scale(g, T(-0.01))
scaleup(g) = scale(g, T(0.01))

const GL_N = 8
const GL_NODES_RAW = SVector{GL_N,T}(
    -0.9602898564975363,
    -0.7966664774136267,
    -0.5255324099163290,
    -0.1834346424956498,
    0.1834346424956498,
    0.5255324099163290,
    0.7966664774136267,
    0.9602898564975363
)
const GL_NODES = ○ .+ GL_NODES_RAW ./ (GL_NODES_RAW[end] - GL_NODES_RAW[1])
const GL_WEIGHTS = SVector{GL_N,T}(
    0.1012285362903763,
    0.2223810344533745,
    0.3137066458778873,
    0.3626837833783620,
    0.3626837833783620,
    0.3137066458778873,
    0.2223810344533745,
    0.1012285362903763,
)
