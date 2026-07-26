module TOGOctahedron

export Octahedron, take!

# using StaticArrays
using KernelAbstractions, LinearAlgebra, Colors, FileIO, ImageIO
using TOGGPU: GPU_BACKEND, GPU_BACKEND_WORKGROUPSIZE
using TOGExist: ○, 𝕋, Ο, ∃, ○̂
using TOGRay: GL_NODES, GL_N
using TOGPrivacy
using LoopOS: Peripheral
using TOGXAI: imageb64totext
using TOGColor: scalar2rgba
import TOGRay: ∃̇
import Base: ∩, take!

"""
Observe an ∞-dim Ω with a 3d manifold Octahedron.
The octahedron spans between the observer and the focus.
All existence is smooth, yet all observation is discrete.
The pyramid from the focus times ρ[end] to the base is discretized and observed via Beer-Lambert accumulation along the pyramid height.
# Arguments
- `t::Number`: Max world time to observe.
- `∂t::Bool`: Fix observation time to the current world time.
- `vt::Int`: Speed of time change per step.
- `v::Number`: Speed for `observer` to move towards `focus`.
- `d::AbstractVector{<:Number}`: Non-zero support dims of `observer` and `focus`.
- `observer::AbstractVector{<:Number}`: Origin peak of the octahedron.
- `focus::AbstractVector{<:Number}`: End peak of the octahedron.
- `ρ::AbstractVector{<:Number}`: 3-dim size vector of the octahedron base.
- `θ::Number`: Rotation around the octahedron height line.
- `♯::NTuple`: 2d discretization of the base.
- `norm::Function`: Measure the length of a vector.
- `⚷::UInt`: Not implemented yet, to be used to encrypt/decrypt information.
"""
mutable struct Octahedron <: Peripheral
    t::Number
    ∂t::Bool
    vt::Number
    v::Number
    d::AbstractVector{<:Number}
    observer::AbstractVector{<:Number}
    focus::AbstractVector{<:Number}
    ρ::AbstractVector{<:Number}
    θ::Number
    ♯::NTuple
    norm::Function
    ⚷::UInt # todo use
end
function Octahedron(; t, d, observer, focus, ∂t=true,vt=0.0,v=0.0,ρ=○*(focus .- observer), θ=0.0, ⚷=zero(UInt), ♯=(2, 2), n̂orm=x -> sqrt(sum(x̃ -> x̃^2, x)))
    # @info "Octahedron"
    @assert length(d) == length(observer) == length(focus)
    Octahedron(t, ∂t, vt, v, d, observer, focus, ρ, θ, ♯, n̂orm, ⚷)
end

"""
Observe layers octahedra given a compositing factor for each layer. Returns a description.
example: `take!([(octahedron1,0.5),(octahedron2,0.5)])`
"""
function take!(oo::AbstractVector{Tuple{Octahedron,Number}})
    img = image(oo)
    imagetotext(img)
end
"""
Observe given an octahedron. Returns a description.
example: `take!(TOGgod.OCTAHEDRON[])`
"""
function take!(o::Octahedron)
    img = image(o)
    imagetotext(img)
end

function imagetotext(img)
    io = IOBuffer()
    save(Stream(format"PNG", io), img)
    b64 = base64encode(take!(io))
    imageb64totext(b64=b64)
end

function image(oo::AbstractVector{Tuple{Octahedron,Number}})
    ϕ, α = rawtake!(oo)
    image(ϕ, α)
end
function image(o::Octahedron)
    ϕ, α = rawtake!(o)
    image(ϕ, α)
end
function image(ϕ, α)
    img = fill(RGBA(0,0,0,0), size(ϕ)...)
    for i = CartesianIndices(ϕ)
        img[i] = scalar2rgba(ϕ[i], α[i])
    end
    img
end

composite(ψ₁, α₁, ψ₂, α₂, β) = (one(β) - β) * ψ₁ .+ β * ψ₂, ○ * (α₁ .+ α₂)
function rawtake!(oo::AbstractVector{Tuple{Octahedron,Number}})
    # @info "Octahedron, take! oo"
    ϕ, α = rawtake!(oo[1][1])
    for i = 2:length(oo)
        ϕ̃, α̃ = rawtake!(oo[i][1])
        ϕ, α = composite(ϕ, α, ϕ̃, α̃, oo[i][2])
    end
    ϕ, α
end
"""
Observe given an octahedron. Returns fuzze existence potential and its accumulation.
example: `rawtake!(TOGgod.OCTAHEDRON[])`
"""
function rawtake!(o::Octahedron)
    N, z, dx, dy, c, a, za, ca, _, ϵμ, ϵρ = pyramid(o)
    ∂ = fill(true, length(o.d))
    ϵ = ∃(o.d, ϵμ, ϵρ, ∂, ∂, ○̂)
    ϵϵϵ = ∩(ϵ, o.t)
    hasdepth = !iszero(o.ρ[end])
    nz = hasdepth ? GL_N - 1 : 1
    i = zeros(Int, o.♯[1], o.♯[2], nz)
    ϵϵ = owners!(i, o, c, ca, dx, dy, ϵϵϵ, nz, N)
    isempty(ϵϵ) && return fill(○, o.♯[1], o.♯[2]), ones(eltype(c), o.♯[1], o.♯[2])
    ϵϵϵϵ = ϵTuple(ntuple(i -> ϵϵ[i], length(ϵϵ)))
    ôneϕ, ôneα = if hasdepth
        ix = a .≠ ○
        ôneϵ = ∃(o.d[ix], a[ix], zeros(eltype(a), sum(ix)), ∂[ix], ∂[ix], ○̂)
        aϵϵϵ = ∩(ôneϵ, o.t)
        ϕ = isempty(aϵϵϵ) ? ○ : only(aϵϵϵ).ϕ(a)
        ϕ, ι(ϕ)
    else
        zero(eltype(a)), one(eltype(a))
    end
    Φ, Α = ΦΨ(i, ϵϵϵϵ, o.d, z, za, dx, dy, N, o.♯[1], o.♯[2], nz)
    hasdepth ? ∃̇(Φ, Α, ôneϕ, ôneα, o.♯[1], o.♯[2], nz) : Φ[:, :, 1], Α[:, :, 1]
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
    h = (ns - 1) * ○
    mid = (ns + 1) * ○
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

pyramid(o::Octahedron) = pyramid(o.observer, o.focus, o.ρ, o.θ, o.norm)
function pyramid(observer, focus, gρ, θ, n̂orm)
    N = length(observer)
    zo = focus .- observer
    D = n̂orm(zo)^2
    perm = sortperm(abs.(zo[(end-1):end]))
    i1, i2 = N - 3 + perm[1], N - 3 + perm[2]
    e1 = zeros(eltype(focus), N)
    e1[i1] = one(eltype(focus))
    d1 = e1 - dot(e1, zo) / D * zo
    d1 /= n̂orm(d1)
    e2 = zeros(eltype(focus), N)
    e2[i2] = one(eltype(focus))
    d2 = e2 - dot(e2, zo) / D * zo
    d2 = d2 - dot(d2, d1) * d1
    d2 /= n̂orm(d2)
    dx = cos(θ) * d1 + sin(θ) * d2
    dy = -sin(θ) * d1 + cos(θ) * d2
    dx *= gρ[end-2]
    dy *= gρ[end-1]
    c = (focus .+ observer) * ○
    ca = gρ[end] * zo * ○
    a = c .+ ca
    μ, ρ = pyramid_aabb(c, a, dx, dy)
    z = c .- dx .- dy
    za = a .- z
    N, z, dx, dy, c, a, za, ca, zo, μ, ρ
end

zo2μρ(z, o) = (o .+ z)*○, (o .- z)*○
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
    x = fill(○, Ndϵ)
    fmatch(idϵ, id) = begin
        z̃ = z[id] + ṫ * za[id]
        d̃x = t̃ * dx[id]
        d̃y = t̃ * dy[id]
        # @info z̃, ix, iy, d̃x, d̃y, id, idϵ, t̃, ṫ, za[id], z[id], dy[id], dx[id]
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
function ΦΨ(i, ϵϵ, d, z, za, dx, dy, N, nx, ny, nz)
    # T = eltype(d)
    T = Float32
    Φ = KernelAbstractions.allocate(GPU_BACKEND, T, nx, ny, nz)
    Α = KernelAbstractions.allocate(GPU_BACKEND, T, nx, ny, nz)
    ΦΨkernel(GPU_BACKEND, GPU_BACKEND_WORKGROUPSIZE)(Φ, Α, i, ϵϵ, d, z, za, dx, dy, N, nx, ny, ndrange=(nx, ny, nz))
    KernelAbstractions.synchronize(GPU_BACKEND)
    Φ, Α
end
@kernel function ΦΨkernel(Φ, Α, i, ϵϵ, d, z, za, dx, dy, N, nx, ny)
    # T = eltype(d)
    (ix, iy, iz) = @index(Global, NTuple)
    iϵ = i[ix, iy, iz]
    if iszero(iϵ)
        Φ[ix, iy, iz] = ○
        # Α[ix, iy, iz] = one(T)
        Α[ix, iy, iz] = 1//1
    else
        # ĩx = T(ix - 1) / T(nx - 1)
        # ĩy = T(iy - 1) / T(ny - 1)
        ĩx = (ix - 1) / (nx - 1)
        ĩy = (iy - 1) / (ny - 1)
        # @info ix, iy, nx, ny, ĩx, ĩy
        ẋ = coordinates(ĩx, ĩy, iz, ϵϵ, d, z, za, dx, dy, N)
        ϕ = ϕ̇(ϵϵ, iϵ, ẋ)
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
        return ○
    end
end

end
