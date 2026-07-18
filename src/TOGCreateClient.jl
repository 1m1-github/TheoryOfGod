module TOGCreateClient

export put!

using ZMQ, PNGFiles, TOGZMQAPIClient, TOGTypst
using TOGOctahedron: Octahedron
import Base: put!

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
    ∂ = fill(true, length(μ))
    create(o, μ, ρ, typst(typst_code=typst_code, dpi=300), ∂, ∂, name)
end

"""
Create a 2d RGBA Matrix with center μ and radius ρ inside Octahedron o.
"""
function put!(o::Octahedron, name::String, mat::AbstractMatrix{PNGFiles.ColorTypes.RGBA}, μ, ρ)
    create(o, μ, ρ, rgbamatrix(mat), ∂, ∂, name)
end

"""
Create ϕ:Ω->I with center μ and radius ρ inside Octahedron o.
"""
function put!(o::Octahedron, name::String, ϕ::Function, μ, ρ)
    create(o, μ, ρ, ϕ, ∂, ∂, name)
end

end
