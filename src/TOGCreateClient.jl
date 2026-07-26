module TOGCreateClient

export put!, TOGCreateClientAdvice

using ZMQ, PNGFiles, TOGZMQAPIClient
using TOGTypst: typst
using TOGOctahedron: Octahedron, pyramid, box_aabb, ∃̇
using TOGExist: ∃
using TOGMatrix: rgbamatrix
import Base: put!

const SOCKET = Ref{Socket}()

# __init__() = atexit(sleep)
# sleep() = TOGZMQAPIClient.sleep(SOCKET[])
awaken(socketlocation) = SOCKET[] = TOGZMQAPIClient.awaken(socketlocation)
∃!call(x...) = TOGZMQAPIClient.call(SOCKET[], :create, x...)

const TOGCreateClientAdvice = """
For all `TOGCreateClient.put!`:
These create entities in Ω, the visual world, running as a separate process reached via ZMQ.
Simplest is to create with relative center μ and radius ρ inside Octahedron o. Alternatively, you can also create an arbitrary ∃.
The μ and ρ should be such that the entire entity is within the Octahedron (in [0,1]^3).
length(μ) == length(ρ).
length(μ) == 2 => Will create on the octahedron base.
length(μ) == 3 => Will create on inside octahedron with 3rd dim along the height towards the focus.
Use dims 4 and above for example to emulate local time, to create animations, to get clear pages, etc. You can change coordinates or the dimension itself. 
If you see an 'Intersection Found' error, it means you tried to create something that overlaps in all dimensions with something existing.
If you create something using a `Pkg` that is not installed in Ω, you should see an error saying that. Let me know and we can install that `Pkg` there. Over time, we can find a better solution.
"""

"""
Create an arbitrary ∃.
"""
function put!(ϵ::∃, name::String)
    ∂ = fill(true, length(μ))
    ∃!call(ϵ, name)
end
"""
Create Typst. Only needs actual content, page is already setup.
"""
function put!(o::Octahedron, name::String, typst_code::String, μ::AbstractVector{<:Number}, ρ::AbstractVector{<:Number})
    ∂ = fill(true, length(μ))
    ∃!(o, μ, ρ, typst(typst_code=typst_code, dpi=300), ∂, ∂, name)
end
"""
Create a 2d RGBA Matrix.
"""
function put!(o::Octahedron, name::String, mat::AbstractMatrix, μ::AbstractVector{<:Number}, ρ::AbstractVector{<:Number})
    ∂ = fill(true, length(μ))
    ∃!(o, μ, ρ, rgbamatrix(mat), ∂, ∂, name)
end
"""
Create an arbitrary ϕ:Ω->I.
"""
function put!(o::Octahedron, name::String, ϕ::Function, μ::AbstractVector{<:Number}, ρ::AbstractVector{<:Number})
    ∂ = fill(true, length(μ))
    ∃!(o, μ, ρ, ϕ, ∂, ∂, name)
end

function ∃!(o::Octahedron, μ, ρ, ϕ, ∂₀, ∂₁, name::String)
    length(μ) == 2 && return ∃!2d(o, μ, ρ, ϕ, ∂₀, ∂₁, name)
    length(μ) == 3 && return ∃!3d(o, μ, ρ, ϕ, ∂₀, ∂₁, name)
    length(μ) == 4 && return ∃!4d(o, μ, ρ, ϕ, ∂₀, ∂₁, name)
    ∃!Nd(o, ϕ, ∂₀, ∂₁, name, ω)
end
function ∃!Nd(o::Octahedron, ϕ, ∂₀, ∂₁, name::String)
    _, _, _, _, _, _, _, _, _, μ̃, ρ̃ = pyramid(o)
    ϵ = ∃(o.d, μ̃, ρ̃, ∂₀, ∂₁, ϕ)
    ∃!call(ϵ, name)
end
function ∃!2d(o::Octahedron, μ, ρ, ϕ, ∂₀, ∂₁, name::String)
    _, z, dx, dy, _, _, _, _, zo, _, _ = pyramid(o)
    μ̃ = z .+ 2 * (μ[1] * dx .+ μ[2] * dy)
    dx̃ = 2 * dx * ρ[1]
    dỹ = 2 * dy * ρ[2]
    dz̃ = zeros(eltype(dx̃),length(dx̃))
    μ̃, ρ̃ = box_aabb(μ̃, [dx̃, dỹ, dz̃])
    ∂̃₀ = [fill(true, length(o.d) - 2)..., ∂₀...]
    ∂̃₁ = [fill(true, length(o.d) - 2)..., ∂₁...]
    ϵ = ∃(o.d, μ̃, ρ̃, ∂̃₀, ∂̃₁, ϕ)
    ∃!call(ϵ, name)
end
function ∃!3d(o::Octahedron, μ, ρ, ϕ, ∂₀, ∂₁, name::String)
    _, z, dx, dy, _, _, za, _, _, _, _ = pyramid(o)
    t̃ = one(eltype(μ)) - μ[3]
    μ̃ = z .+ μ[3] * za .+ 2 * (μ[1] * t̃ * dx .+ μ[2] * t̃ * dy)
    ρ̃ = zeros(eltype(μ̃), size(μ̃)...)
    ρ̃[2] = 2 * o.norm(dx) * ρ[1] * t̃ * min(μ[1], one(eltype(μ)) - μ[1])
    ρ̃[3] = 2 * o.norm(dy) * ρ[2] * t̃ * min(μ[2], one(eltype(μ)) - μ[2])
    ρ̃[4] = o.norm(za) * ρ[3] * min(μ[3], (one(eltype(μ)) - max(μ[1], μ[2])) * t̃)
    ϵ = ∃(o.d, μ̃, ρ̃, [true, ∂₀...], [true, ∂₁...], ϕ)
    ∃!call(ϵ, name)
end
function ∃!4d(o::Octahedron, μ, ρ, ϕ, ∂₀, ∂₁, name::String)
    _, _, _, _, _, _, _, _, _, μ̃, ρ̃ = pyramid(o)
    ϵ = ∃(o.d, μ̃, ρ̃, ∂₀, ∂₁, ϕ)
    ∃!call(ϵ, name)
end

end
