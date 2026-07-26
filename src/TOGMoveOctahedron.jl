module TOGMoveOctahedron

# using LoopOS: @whiletrue
# using TOGZMQAPIServer: push!
using TOGOctahedron: Octahedron, pyramid
using TOGExist
using TOGExist: 𝕋

δN = 1
δ(o::Octahedron) = δN * 0.001

# function awaken(;socket, o::Octahedron, ω::𝕋)
# function awaken(;socket, o::Octahedron)
    # TOGZMQAPIServer.start(socket)
    # push!(μ -> move!(o, μ))
    # push!(μ -> focus!(o, μ))
    # push!(ρ -> scale!(o, ρ))
    # push!(θ -> rotate!(o, θ))
    # push!(v -> speed!(o, v))
    # push!(δ -> accelerate!(o, δ))
    # push!(δ -> jerk!(o, δ))
    # push!(_ -> present!(o, ω))
    # push!(_ -> freezetime!(o))
    # push!(δ -> speeduptime!(o, δ))
    # push!(δ -> slowdowntime!(o, δ))
    # @async @whiletrue step!(o, ω)
# end

function step!(o::Octahedron, ω::𝕋)
    o.∂t ? present!(o, ω) : o.t += o.vt
    μ = o.observer .+ o.v * (o.observer .- o.focus)
    move!(g, μ)
end
valid(o::Octahedron) = valid(o.observer, o.focus, o.ρ, o.θ, o.norm)
function valid(observer, focus, ρ, θ, norm)
    _, _, _, _, _, _, _, _, _, μ, ρ = pyramid(observer, focus, ρ, θ, norm)
    all(zero(eltype(μ)) .≤ μ .- ρ .≤ μ .+ ρ .≤ one(eltype(μ)))
end
function move!(o::Octahedron, observer)
    # @info "move!", ẑeroμ
    valid(observer, o.focus, o.ρ, o.θ, o.norm) || return
    # @info "valid"
    o.observer = observer
end
function focus!(o::Octahedron, focus)
    # @info "focus!", focus
    valid(o.observer, focus, o.ρ, o.θ, o.norm) || return
    # @info "valid"
    o.focus = focus
end
function scale!(o::Octahedron, ρ)
    # @info "scale!", ρ
    valid(o) || return
    # @info "valid"
    o.ρ = ρ
end
function rotate!(o::Octahedron, θ)
    # @info "rotate!", θ
    valid(o) || return
    # @info "valid"
    g.θ = θ
end

# scale!(o::Octahedron, i, δ) = scale!(o, ntuple(ĩ -> begin
#         ĩ == i && return o.ρ[ĩ] + δ
#         o.ρ[ĩ]
#     end, length(o.ρ)))
scale!(o::Octahedron, i, δ) = scale!(o, [ĩ == i ? o.ρ[ĩ] + δ : o.ρ[ĩ] for ĩ = 1:length(o.ρ)])
# move!(o::Octahedron, i, δ) =
#     move!(o, SVector(ntuple(ĩ -> begin
#             ĩ == i && return o.ẑeroμ[ĩ] + δ
#             o.ẑeroμ[ĩ]
#         end, length(o.ẑeroμ))))
move!(o::Octahedron, i, δ) = move!(o, [ĩ == i ? o.observer[ĩ] + δ : o.observer[ĩ] for ĩ = 1:length(o.observer)])
# focus!(o::Octahedron, i, δ) = focus!(o, SVector(ntuple(ĩ -> begin
#         ĩ == i && return o.focus[ĩ] + δ
#         o.focus[ĩ]
#     end, length(o.focus))))
focus!(o::Octahedron, i, δ) = focus!(o, [ĩ == i ? o.focus[ĩ] + δ : o.focus[ĩ] for ĩ = 1:length(o.focus)])

focusup!(o::Octahedron, i) = focus!(o, i, δ(o))
focusdown!(o::Octahedron, i) = focus!(o, i, -δ(o))
moveup!(o::Octahedron, i) = move!(o, i, δ(o))
movedown!(o::Octahedron, i) = move!(o, i, -δ(o))
jerkup!(o::Octahedron) = jerk!(o, δ(o))
jerkdown!(o::Octahedron) = jerk!(o, -δ(o))
scaleup!(o::Octahedron, i) = scale!(o, i, δ(o))
scaledown!(o::Octahedron, i) = scale!(o, i, -δ(o))
rotateup!(o::Octahedron) = rotate!(o, o.θ + δ(o))
rotatedown!(o::Octahedron) = rotate!(o, o.θ - δ(o))

speed!(o::Octahedron, v) = o.v = clamp(v, zero(o.v), one(o.v))
accelerate!(o::Octahedron, δ) = speed!(o, iszero(o.v) ? δ : o.v * exp(δ))
jerk!(o::Octahedron, δ) = accelerate!(o, o.v * exp(δ))

present!(o::Octahedron, ω::𝕋) = o.t = TOGExist.t(ω.Ο[ω])
freezetime!(o::Octahedron) = o.vt = zero(o.vt)
speeduptime!(o::Octahedron, δ) = o.vt += δ
slowdowntime!(o::Octahedron, δ) = o.vt -= δ

end
