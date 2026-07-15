module TOGCreateServer

export TOG

# using StaticArrays
using ZMQ
using TOGZMQAPIServer
using TOGOctahedron: Octahedron, pyramid, box_aabb, ∃̇
using TOG∃: ∃
import TOG∃.∃!

const SOCKET = Ref{Socket}()
const TASK = Ref{Task}()

# __init__() = atexit(sleep)
# function sleep()
#     schedule(TASK[], InterruptException(), error=true)
#     TOGZMQAPIServer.sleep(SOCKET[])
# end
function awaken(;socketlocation, ω)
    SOCKET[], TASK[] = TOGZMQAPIServer.awaken(socketlocation=socketlocation, functions=Dict(
        :create => create(ω),
    ))
end

create(ω) = (x...) -> create(x..., ω)
function create(o::Octahedron, μ, ρ, ϕ, ∂₀, ∂₁, n, ω)
    # @info "create", length(μ)
    length(μ) == 2 && return ∃!2d(o, μ, ρ, ϕ, ∂₀, ∂₁, n, ω)
    length(μ) == 3 && return ∃!3d(o, μ, ρ, ϕ, ∂₀, ∂₁, n, ω)
    length(μ) == 4 && return ∃!4d(o, μ, ρ, ϕ, ∂₀, ∂₁, n, ω)
    ∃!(o, ϕ, ∂₀, ∂₁, n, ω)
end

function ∃!(o::Octahedron, ϕ, ∂₀, ∂₁, n, ω=o.Ω)
    _, _, _, _, _, _, _, _, _, μ̃, ρ̃ = pyramid(o)
    # N, z, dx, dy, c, a, za, ca, zo, μ, ρ
    ∃!(∃(o.d, μ̃, ρ̃, ∂₀, ∂₁, ϕ), n, ω)
end
function ∃!2d(o::Octahedron, μ, ρ, ϕ, ∂₀, ∂₁, n, ω=o.Ω)
    _, z, dx, dy, _, _, _, _, zo, _, _ = pyramid(o)
    # N, z, dx, dy, c, a, za, ca, zo, μ, ρ
    # @info "∃!2d, z, dx, dy, zo", z, dx, dy, zo
    μ̃ = z .+ 2 * (μ[1] * dx .+ μ[2] * dy)
    dx̃ = 2 * dx * ρ[1]
    dỹ = 2 * dy * ρ[2]
    # dz̃ = eps(eltype(μ)) / o.norm(zo) * zo # todo ? 2d possible now
    dz̃ = zeros(eltype(dx̃),length(dx̃))
    # @show "∃!2d", eps(eltype(μ)), o.norm(zo)
    # @info "∃!2d, μ̃, dx̃, dỹ, dz̃", μ̃, dx̃, dỹ, dz̃
    # @show "∃!2d",  typeof(o.d), typeof(μ), typeof(μ̃)
    # μ̃, ρ̃ = box_aabb(μ̃, SA[dx̃, dỹ, dz̃])
    μ̃, ρ̃ = box_aabb(μ̃, [dx̃, dỹ, dz̃])
    # @info "∃!2d, μ̃, ρ̃", μ̃, ρ̃
    # @show "∃!2d",  typeof(o.d), typeof(μ), typeof(μ̃)
    ∂̃₀ = [fill(true, length(o.d) - 2)..., ∂₀...]
    ∂̃₁ = [fill(true, length(o.d) - 2)..., ∂₁...]
    # @info "∃!2d, ∂̃₀, ∂̃₁", ∂̃₀, ∂̃₁
    ϵ = ∃(o.d, μ̃, ρ̃, ∂̃₀, ∂̃₁, ϕ)
    # @info "∃!2d, ϵ, n", ϵ, n
    ∃!(ϵ, n, ω)
end
function ∃!3d(o::Octahedron, μ, ρ, ϕ, ∂₀, ∂₁, n, ω=o.Ω)
    _, z, dx, dy, _, _, za, _, _, _, _ = pyramid(o)
    # N, z, dx, dy, c, a, za, ca, zo, μ, ρ
    t̃ = one(eltype(μ)) - μ[3]
    μ̃ = z .+ μ[3] * za .+ 2 * (μ[1] * t̃ * dx .+ μ[2] * t̃ * dy)
    ρ̃ = zeros(μ̃)
    ρ̃[2] = 2 * o.norm(dx) * ρ[1] * t̃ * min(μ[1], one(eltype(μ)) - μ[1])
    ρ̃[3] = 2 * o.norm(dy) * ρ[2] * t̃ * min(μ[2], one(eltype(μ)) - μ[2])
    ρ̃[4] = o.norm(za) * ρ[3] * min(μ[3], (one(eltype(μ)) - max(μ[1], μ[2])) * t̃)
    ∃!(∃(o.d, μ̃, ρ̃, [true, ∂₀...], [true, ∂₁...], ϕ), n, ω)
end
function ∃!4d(o::Octahedron, μ, ρ, ϕ, ∂₀, ∂₁, n, ω=o.Ω)
    _, _, _, _, _, _, _, _, _, μ̃, ρ̃ = pyramid(o)
    # N, z, dx, dy, c, a, za, ca, zo, μ, ρ
    ∃!(∃(o.d, μ̃, ρ̃, ∂₀, ∂₁, ϕ), n, ω)
end

# """
# Will show 2d Typst in TOG with center μ and radius ρ. Only needs actual content, page is already setup. 
# """
# put!(::Type{TOGCreate}, typst_code::String, μ::SVector{2,T}=SA[○, ○], ρ::SVector{2,T}=SA[○, ○]) = ∃!2d(typst(typst_code), μ, ρ)
# """
# Will show 2d RGBA matrix in TOG with center μ and radius ρ.
# """
# put!(::Type{TOG}, mat::AbstractMatrix{PNGFiles.ColorTypes.RGBA}, μ::SVector{2,T}=SA[○, ○], ρ::SVector{2,T}=SA[○, ○]) = ∃!2d(rgbamatrix(mat), μ, ρ)
# """
# Will show a 2d ϕ:(t,x,y)->[0,1] in TOG with center μ and radius ρ.
# """
# put!(::Type{TOG}, ϕ::Function, μ::SVector{2,T}=SA[○, ○], ρ::SVector{2,T}=SA[○, ○]) = ∃!2d(ϕ, μ, ρ)
# # put!(::Type{TOG}, ϕ::Function, μ::SVector{2,T}=SA[○, ○], ρ::SVector{2,T}=SA[○, ○]) = ∃!2d(x -> ϕ((x[1], x[2], x[3])), μ, ρ)
# """Will travel to a random location in TOG essentially "clearing" the view."""
# function put!(::Type{TOG})
#     g = godBROWSER[].g
#     δg = g.ône.μ[end] - g.ẑero.μ[end]
#     μ̇1 = SA[g.ẑero.μ[1], rand(T, 3)...]
#     μ̇2 = SA[μ̇1[1], μ̇1[2], μ̇1[3], μ̇1[4]+δg]
#     while !valid(μ̇1, μ̇2, g.ρ, g.θ, g.norm)
#         yield()
#         μ̇1 = SA[g.ẑero.μ[1], rand(T, 3)...]
#         μ̇2 = SA[μ̇1[1], μ̇1[2], μ̇1[3], μ̇1[4]+δg]
#     end
#     move!(g, μ̇1)
#     focus!(g, μ̇2)
# end

end
