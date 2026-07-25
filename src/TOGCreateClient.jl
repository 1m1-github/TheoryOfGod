module TOGCreateClient

export put!, TOGCreateClientAdvice

using ZMQ, PNGFiles, TOGZMQAPIClient
using TOGTypst: typst
using TOGOctahedron: Octahedron
using TOGMatrix: rgbamatrix
import Base: put!

const SOCKET = Ref{Socket}()

# __init__() = atexit(sleep)
# sleep() = TOGZMQAPIClient.sleep(SOCKET[])
awaken(socketlocation) = SOCKET[] = TOGZMQAPIClient.awaken(socketlocation)
create(x...) = TOGZMQAPIClient.call(SOCKET[], :create, x...)

const TOGCreateClientAdvice = """
For all `TOGCreateClient.put!`:
Create with relative center μ and radius ρ inside Octahedron o.
The μ and ρ should be such that the entire entity is within the Octahedron (in [0,1]^3).
length(μ) == length(ρ).
length(μ) == 2 => Will create on the octahedron base.
length(μ) == 3 => Will create on inside octahedron with 3rd dim along the height towards the focus.
Use dims 4 and above for example to emulate local time.
"""

"""
Create Typst. Only needs actual content, page is already setup.
"""
function put!(o::Octahedron, name::String, typst_code::String, μ::AbstractVector{<:Number}, ρ::AbstractVector{<:Number})
    # @info "TOGCreateClient.jl, put!"
    ∂ = fill(true, length(μ))
    create(o, μ, ρ, typst(typst_code=typst_code, dpi=300), ∂, ∂, name)
end

"""
Create a 2d RGBA Matrix.
"""
function put!(o::Octahedron, name::String, mat::AbstractMatrix, μ::AbstractVector{<:Number}, ρ::AbstractVector{<:Number})
    # @info "TOGCreateClient.jl, put!"
    ∂ = fill(true, length(μ))
    create(o, μ, ρ, rgbamatrix(mat), ∂, ∂, name)
end

"""
Create an arbitrary ϕ:Ω->I.
"""
function put!(o::Octahedron, name::String, ϕ::Function, μ::AbstractVector{<:Number}, ρ::AbstractVector{<:Number})
    # @info "TOGCreateClient.jl, put!"
    ∂ = fill(true, length(μ))
    create(o, μ, ρ, ϕ, ∂, ∂, name)
end

end
