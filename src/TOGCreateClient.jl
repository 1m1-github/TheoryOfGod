module TOGCreateClient

export put!, TOGExplanation

using ZMQ, PNGFiles, TOGZMQAPIClient, TOGTypst
using TOGOctahedron: Octahedron
import Base: put!

const TOGExplanation = """
I = [ZERO < ○ < ONE] denotes a unit 1-dim space of information with origin ○ (no information) at its center including the corners ZERO and ONE.
Ω = I^I an ∞-dim normed and smooth vector space.
We have an ontology 𝕋 on Ω such that ϵ ∈ 𝕋:
* ϵ ⊆ Ω
* ϵ ∩ ϵ′ ≠ ∅ => ϵ = ϵ′
* ϵ.ϕ: Ω -> I is arbitrary, computable and smooth fuzzy existence potential towards ONE=true xor ZERO=false.
ϵ ⊊ Ω defines its existence inside a subset of Ω using an origin (μ), a radius (ρ) and a closed vs. open in each direction (∂) vector. These vectors are finite and all other dimensional coordinates of ϵ are ○.
god ⊊ GOD = Ω = I^I = I^(.) = [ZERO < ○ < ONE]^(.)
god observes or creates, GOD iterates.
This is a visual world for you to live in. Use your TOGgod.OCTAHEDRON[] to create and observe. Show yourself, your face, communicate only using this world. I see your pov (octahedron) via a browser.
"""

const SOCKET = Ref{Socket}()

# __init__() = atexit(sleep)
awaken(socketlocation) = SOCKET[] = TOGZMQAPIClient.awaken(socketlocation)
create(x...) = TOGZMQAPIClient.call(SOCKET[], :create, x...)
# sleep() = TOGZMQAPIClient.sleep(SOCKET[])

"""
Create Typst with center μ and radius ρ inside Octahedron o.
Only needs actual content, page is already setup.
"""
function put!(o::Octahedron, name::String, typst_code::String, μ, ρ)
    @info "TOGCreateClient.jl, put!"
    ∂ = fill(true, length(μ))
    create(o, μ, ρ, typst(typst_code=typst_code, dpi=300), ∂, ∂, name)
end

"""
Create a 2d RGBA Matrix with center μ and radius ρ inside Octahedron o.
"""
function put!(o::Octahedron, name::String, mat::AbstractMatrix{PNGFiles.ColorTypes.RGBA}, μ, ρ)
    @info "TOGCreateClient.jl, put!"
    create(o, μ, ρ, rgbamatrix(mat), ∂, ∂, name)
end

"""
Create ϕ:Ω->I with center μ and radius ρ inside Octahedron o.
"""
function put!(o::Octahedron, name::String, ϕ::Function, μ, ρ)
    @info "TOGCreateClient.jl, put!"
    create(o, μ, ρ, ϕ, ∂, ∂, name)
end

end
