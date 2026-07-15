module TOGOctahedron

export Octahedron

# using StaticArrays
using KernelAbstractions, LinearAlgebra
using TOGGPU: GPU_BACKEND, GPU_BACKEND_WORKGROUPSIZE
using TOG∃: ○, 𝕋, Ο, ∃, ○̂
using TOGRay: GL_NODES, GL_N
import TOGRay.∃̇
using TOGPrivacy
using LoopOS: Peripheral
import Base: ∩
import Base.take!

mutable struct Octahedron{T} <: Peripheral
    Ο::UInt
    ∂Ο::Bool
    vΟ::Int
    v::T
    d::AbstractVector{T}
    ẑeroμ::AbstractVector{T}
    ôneμ::AbstractVector{T}
    ρ::AbstractVector{T}
    θ::T
    ♯::NTuple
    norm::Function
    ⚷::UInt # todo use
end
function Octahedron(; t, d, ẑeroμ, ôneμ, ρ=○(eltype(d))*(ôneμ .- ẑeroμ), θ=zero(typeof(t)), ⚷=zero(UInt), ♯=(2, 2), n̂orm=x -> sqrt(sum(x̃ -> x̃^2, x)))
    @assert length(d) == length(ẑeroμ) == length(ôneμ)
    Octahedron{typeof(t)}(Ο(t), true, zero(UInt), zero(t), d, ẑeroμ, ôneμ, ρ, θ, ♯, n̂orm, ⚷)
end

# function trivial(ϵ)
#     isempty(ϵ) && return true
#     1 < length(ϵ) && return false
#     only(ϵ).ϕ === ○̂
# end
function take!(oo::AbstractVector{Tuple{Octahedron,𝕋,T}}) where T
    ψ, α = take!(oo[1][1], oo[1][2])
    for i = 2:length(oo)
        ψ̃, α̃ = take!(oo[i][1], oo[i][2])
        ψ, α = take!(ψ, α, ψ̃, α̃, oo[i][3])
    end
    ψ, α
end
take!(oo::AbstractVector{Tuple{Octahedron,T}}, ω::𝕋) where T = take!([(o[1], ω, o[2]) for o = oo])
take!(oo::AbstractVector{Octahedron}, ω::𝕋, β::Number) = take!([(o, ω, β) for o = oo])
take!(ψ₁, α₁, ψ₂, α₂, β) = (one(β) - β) * ψ₁ .+ β * ψ₂, ○(β) * (α₁ .+ α₂)
function take!(o::Octahedron)
    # @info "∃̇Octahedron(o::Octahedron)"
    # try
    N, z, dx, dy, c, a, za, ca, _, ϵμ, ϵρ = pyramid(o)
    # @info "∃̇Octahedron", N, z, dx, dy, c, a, za, ca, ϵμ, ϵρ
    ∂ = fill(true, length(o.d))
    ϵ = ∃(o.d, ϵμ, ϵρ, ∂, ∂, ○̂)
    ϵϵϵ = ∩(ϵ, o.Ο)
    # ϵϵϵ = β̂(∃(o.d, ϵμ, ϵρ, ∂, ∂, ○̂), ω)
    # @info "∃̇Octahedron, length(ϵϵϵ), map(ϵ->ϵ.μ,collect(ϵϵϵ))", length(ϵϵϵ), map(ϵ->ϵ.μ, collect(ϵϵϵ))
    # istrivial = trivial(ϵϵϵ)
    # @info "∃̇Octahedron", istrivial
    hasdepth = !iszero(o.ρ[end])
    # @info "∃̇Octahedron, hasdepth", hasdepth
    nz = hasdepth ? GL_N - 1 : 1
    # i = fill(istrivial ? 0 : 1, o.♯..., nz)
    i = zeros(Int, o.♯[1], o.♯[2], nz)
    ϵϵ = owners!(i, o, c, ca, dx, dy, ϵϵϵ, nz, N)
    # @info "∃̇Octahedron, unique(i), length(ϵϵ)", unique(i), length(ϵϵ)
    # @info "∃̇Octahedron, i", i
    # @info "∃̇Octahedron", i
    isempty(ϵϵ) && return fill(○(c), o.♯[1], o.♯[2]), ones(eltype(c), o.♯[1], o.♯[2])
    ϵϵϵϵ = ϵTuple(ntuple(i -> ϵϵ[i], length(ϵϵ)))
    ôneϕ, ôneα = if hasdepth
        # z = @SVector zeros(c, N)
        # ôneϵ = β(o.d, a, z, ∂, ω)
        # aϵϵϵ = β̂(∃(o.d, a, zeros(typeof(a), N), ∂, ∂, ○̂), ω)
        ix = a .≠ ○(eltype(a))
        ôneϵ = ∃(o.d[ix], a[ix], zeros(eltype(a), sum(ix)), ∂[ix], ∂[ix], ○̂)
        aϵϵϵ = ∩(ôneϵ, o.Ο)
        ϕ = isempty(aϵϵϵ) ? ○(eltype(a)) : only(aϵϵϵ).ϕ(a)
        ϕ, ι(ϕ)
    else
        zero(eltype(a)), one(eltype(a))
    end
    # @info "∃̇Octahedron, ôneϕ", ôneϕ
    Φ, Α = ΦΨ(i, ϵϵϵϵ, o.d, z, za, dx, dy, N, o.♯[1], o.♯[2], nz)
    # @info "∃̇Octahedron", Φ, Α
    hasdepth ? ∃̇(Φ, Α, ôneϕ, ôneα, o.♯[1], o.♯[2], nz) : Φ[:, :, 1], Α[:, :, 1]
    # catch e
    #     bt = catch_backtrace()
    #     showerror(stderr, e, bt)
    #     fill(○(c), o.♯...)
    # end
end
function owners!(ϵ̂, o, c, ca, dx, dy, ϵϵϵ, nz, N)
    ϵϵ = ∃[]
    for ϵ = ϵϵϵ
        ∩(
            ϵ̂, length(ϵϵ) + 1,
            c, ca,
            dx, dy,
            ϵ.μ .- ϵ.ρ, ϵ.μ .+ ϵ.ρ, ϵ.∂₀, ϵ.∂₁,
            o.♯[1], o.♯[2], nz, N) && push!(ϵϵ, ϵ)
    end
    ϵϵ
end

@inline function tighten(lo, hi, coeff, abs_other, L, U)
    if zero(lo) < coeff
        lo = max(lo, (L - abs_other) / coeff)
        hi = min(hi, (U + abs_other) / coeff)
    elseif coeff < zero(lo)
        lo = max(lo, (U + abs_other) / coeff)
        hi = min(hi, (L - abs_other) / coeff)
    else
        (abs_other < L || U < -abs_other) && return (lo, hi, true)
    end
    return (lo, hi, hi < lo)
end
@inline function s_to_indices(slo, shi, ns)
    h = (ns - 1) * ○(slo)
    mid = (ns + 1) * ○(slo)
    ilo = clamp(ceil(Int, slo * h + mid), 1, ns)
    ihi = clamp(floor(Int, shi * h + mid), 1, ns)
    ihi < ilo && return (0, 0)
    return (ilo, ihi)
end
function ∩̇(i, i̇, c, ca, dx, dy, ϵlo, ϵhi, nx, ny, N, k)
    t = GL_NODES(eltype(c))[k]
    sx_lo, sx_hi = -one(t), one(t)
    sy_lo, sy_hi = -one(t), one(t)
    empty = false
    for m = 1:N
        ct = c[m] + t * ca[m]
        dxm = dx[m]
        dym = dy[m]
        L = ϵlo[m] - ct
        U = ϵhi[m] - ct
        sx_lo, sx_hi, empty = tighten(sx_lo, sx_hi, dxm, abs(dym), L, U)
        empty && break
        sy_lo, sy_hi, empty = tighten(sy_lo, sy_hi, dym, abs(dxm), L, U)
        empty && break
    end
    empty && return false
    xlo, xhi = s_to_indices(sx_lo, sx_hi, nx)
    ylo, yhi = s_to_indices(sy_lo, sy_hi, ny)
    (iszero(xlo) || iszero(ylo)) && return false
    # @info "∩̇, xlo:xhi, ylo:yhi, k", xlo:xhi, ylo:yhi, k
    i[xlo:xhi, ylo:yhi, k] .= i̇
    true
end
function ∩(i, i̇, c, ca, dx, dy, ϵlo, ϵhi, ϵ∂₀, ϵ∂₁, nx, ny, nz, N)
    intersects = false
    @inbounds for k = 1:nz
        intersects |= ∩̇(i, i̇, c, ca, dx, dy, ϵlo, ϵhi, nx, ny, N, k)
    end
    intersects
end

pyramid(o::Octahedron) = pyramid(o.ẑeroμ, o.ôneμ, o.ρ, o.θ, o.norm)
# using StaticArrays, LinearAlgebra
# ẑeroμ, ôneμ, gρ, θ = [0.0, 0.5, 0.5, 0.5], [0.0, 0.5, 0.5, 0.6], [0.1, 0.1, 0.1, 1.0], 0.0
# n̂orm=x -> sqrt(sum(x̃ -> x̃^2, x))
function pyramid(ẑeroμ, ôneμ, gρ, θ, n̂orm)
    # @show "pyramid", ẑeroμ, ôneμ, gρ, θ
    N = length(ẑeroμ)
    zo = ôneμ .- ẑeroμ
    D = n̂orm(zo)^2
    # @show "pyramid", zo, D
    perm = sortperm(abs.(zo[(end-1):end]))
    i1, i2 = N - 3 + perm[1], N - 3 + perm[2]
    # e1 = zeros(SVector{N,eltype(ôneμ)})
    e1 = zeros(eltype(ôneμ), N)
    e1[i1] = one(eltype(ôneμ))
    # e1 = setindex(e1, one(eltype(ôneμ)), i1)
    d1 = e1 - dot(e1, zo) / D * zo
    d1 /= n̂orm(d1)
    # e2 = zeros(SVector{N,eltype(ôneμ)})
    e2 = zeros(eltype(ôneμ), N)
    e2[i2] = one(eltype(ôneμ))
    # e2 = setindex(e2, one(eltype(ôneμ)), i2)
    d2 = e2 - dot(e2, zo) / D * zo
    d2 = d2 - dot(d2, d1) * d1
    d2 /= n̂orm(d2)
    # @show "pyramid", d1, d2
    dx = cos(θ) * d1 + sin(θ) * d2
    dy = -sin(θ) * d1 + cos(θ) * d2
    dx *= gρ[end-2]
    dy *= gρ[end-1]
    # @show "pyramid", dx, dy
    c = (ôneμ .+ ẑeroμ) * ○(ôneμ)
    ca = gρ[end] * zo * ○(ôneμ)
    # ca = gρ[end] / sqrt(D) * zo
    a = c .+ ca
    # @show "pyramid", c, ca, a
    μ, ρ = pyramid_aabb(c, a, dx, dy)
    # @show "pyramid", μ, ρ
    z = c .- dx .- dy
    za = a .- z
    # @show "pyramid", z, za
    N, z, dx, dy, c, a, za, ca, zo, μ, ρ
end

zo2μρ(z, o) = (o .+ z)*○(z), (o .- z)*○(z)
function pyramid_aabb(c, a, dx, dy)
    half = abs.(dx) .+ abs.(dy)
    z = min.(a, c .- half)
    o = max.(a, c .+ half)
    zo2μρ(z, o)
end
function box_aabb(c, v)
    half = sum(x -> abs.(x), v)
    z = c .- half
    o = c .+ half
    zo2μρ(z, o)
end

function merge_match(a, b, na, nb, fmatch)
    i, j = 1, 1
    while i <= na && j <= nb
        if a[i] == b[j]
            fmatch(i, j)
            i += 1
            j += 1
        elseif a[i] < b[j]
            i += 1
        else
            j += 1
        end
    end
end
function coordinates(ix, iy, iz, ϵϵ, d, z, za, dx, dy, N)
    ṫ = GL_NODES(eltype(d))[iz]
    t̃ = one(eltype(d)) - ṫ
    Ndϵ = length(ϵϵ.ϵ[1].d) # todo correct? why [1]
    x = fill(○(z), Ndϵ)
    fmatch(idϵ, id) = begin
        z̃ = z[id] + ṫ * za[id]
        d̃x = t̃ * dx[id]
        d̃y = t̃ * dy[id]
        x[idϵ] = z̃ + 2 * (ix * d̃x + iy * d̃y)
    end
    merge_match(ϵϵ.ϵ[1].d, d, Ndϵ, N, fmatch)
    # SVector{Ndϵ,typeof(d)}(x)
    x
end
function ι(ϕ)
    (iszero(ϕ) || isone(ϕ)) && return one(ϕ)
    ϕ̄ = one(ϕ) - ϕ
    one(ϕ) + (ϕ * log(ϕ) + ϕ̄ * log(ϕ̄)) / log(2)
end
# function ΦΨ(i, ΦΦ, d, z, za, dx, dy, nx, ny, nz, N)
#     Φ = fill(○(z), nx, ny, nz)
#     Α = ones(eltype(d), nx, ny, nz)
#     # Threads.@threads
#     for (ix, iy) in collect(Iterators.product(1:nx, 1:ny))
#         ĩx = eltype(d)(ix - 1) / eltype(d)(nx - 1)
#         ĩy = eltype(d)(iy - 1) / eltype(d)(ny - 1)
#         for iz = 1:nz
#             iϕ = i[ix, iy, iz]
#             iszero(iϕ) && continue
#             ẋ = coordinates(ĩx, ĩy, iz, ΦΦ, d, z, za, dx, dy, N)
#             @show ẋ
#             ϕ = Φ̇(ΦΦ, iϕ, ẋ)
#             Φ[ix, iy, iz] = ϕ
#             Α[ix, iy, iz] = ι(ϕ)
#         end
#     end
#     Φ, Α
# end
function ΦΨ(i, ϵϵ, d, z, za, dx, dy, N, nx, ny, nz)
    T = eltype(d)
    # @info nx, ny, nz
    Φ = KernelAbstractions.allocate(GPU_BACKEND, T, nx, ny, nz)
    Α = KernelAbstractions.allocate(GPU_BACKEND, T, nx, ny, nz)
    ΦΨkernel(GPU_BACKEND, GPU_BACKEND_WORKGROUPSIZE)(Φ, Α, i, ϵϵ, d, z, za, dx, dy, N, nx, ny, ndrange=(nx, ny, nz))
    KernelAbstractions.synchronize(GPU_BACKEND)
    Φ, Α
end
@kernel function ΦΨkernel(Φ, Α, i, ϵϵ, d, z, za, dx, dy, N, nx, ny)
    T = eltype(d)
    (ix, iy, iz) = @index(Global, NTuple)
    iϵ = i[ix, iy, iz]
    # @info "ΦΨkernel", iϵ, ix, iy, iz
    if iszero(iϵ)
        Φ[ix, iy, iz] = ○(T)
        Α[ix, iy, iz] = one(T)
    else
        ĩx = T(ix - 1) / T(nx - 1)
        ĩy = T(iy - 1) / T(ny - 1)
        ẋ = coordinates(ĩx, ĩy, iz, ϵϵ, d, z, za, dx, dy, N)
        ϕ = ϕ̇(ϵϵ, iϵ, ẋ)
        # @info "ΦΨkernel", iϵ, ix, iy, ẋ, ϕ, ι(ϕ)
        Φ[ix, iy, iz] = ϕ
        Α[ix, iy, iz] = ι(ϕ)
    end
end
struct ϵTuple{ϵT}
    ϵ::ϵT
end
@generated function ϕ̇(ϵϵ::ϵTuple{ϵT}, i, x) where ϵT
    N = length(ϵT.parameters)
    branches = []
    for ĩ = 1:N
        push!(branches, quote
            if i == $ĩ
                return ϵϵ.ϵ[$ĩ].ϕ(x)
            end
        end)
    end
    quote
        $(branches...)
        return ○(x)
    end
end

end
