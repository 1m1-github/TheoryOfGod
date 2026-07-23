module TOGRay

using FastGaussQuadrature, KernelAbstractions
using TOGGPU: GPU_BACKEND, GPU_BACKEND_WORKGROUPSIZE
using TOGExist: ○
const GL_N = 8
const GL_NODES_RAW, GL_WEIGHTS_RAW = gausslegendre(GL_N)
const GL_NODES(T) = ○(T) .+ GL_NODES_RAW ./ (GL_NODES_RAW[end] - GL_NODES_RAW[1]) # [0,1]
const GL_WEIGHTS(T) = ○(T) * GL_WEIGHTS_RAW

# function ∃̇(i, Φ, Α, ôneϕ, nx, ny, nz)
#     ψ = zeros(typeof(ôneϕ), nx, ny)
#     α = ones(typeof(ôneϕ), nx, ny)
#     Τ = -ones(typeof(ôneϕ), nx, ny, nz)
#     # Threads.@threads 
#     for (ix, iy) = collect(Iterators.product(1:nx, 1:ny))
#         ∃̇!(
#             ψ, α, Τ, Φ, Α,
#             ix, iy, i, ôneϕ,
#             nz
#         )
#     end
#     ψ, α
# end
# function ∃̇!(ψ, α, Τ, Φ, Α, ix, iy, i, ôneϕ, nz)
#     ψ1 = ψ2 = τ = zero(ôneϕ)
#     # hasdepth = 1 < nz
#     for iz = 1:nz
#         # w = hasdepth ? GL_WEIGHTS(typeof(ôneϕ))[iz] : one(ôneϕ)
#         w = GL_WEIGHTS(typeof(ôneϕ))[iz]
#         iϕ = i[ix, iy, iz]
#         if iszero(iϕ)
#             ψ1 += ○(ôneϕ) * w
#             ψ2 += w
#             continue
#         end
#         if zero(ôneϕ) ≤ Τ[ix, iy, iz]
#             τ = Τ[ix, iy, iz]
#         else
#             τ = zero(ôneϕ)
#             for jz = 1:iz
#                 τ += Α[ix, iy, jz] * GL_WEIGHTS(typeof(ôneϕ))[jz]
#             end
#             Τ[ix, iy, iz] = τ
#         end
#         cτ = Α[ix, iy, iz] * exp(-τ) * w
#         ψ1 += cτ * Φ[ix, iy, iz]
#         ψ2 += cτ
#     end
#     expτ = exp(-τ)
#     ψ1 += expτ * ôneϕ
#     ψ2 += expτ
#     # ψ[ix, iy] = hasdepth ? ψ1 / ψ2 : ôneϕ
#     # α[ix, iy] = hasdepth ? one(ôneϕ) - expτ : one(ôneϕ)
#     ψ[ix, iy] = ψ1 / ψ2
#     α[ix, iy] = one(ôneϕ) - expτ
# end
function ∃̇(Φ, Α, ôneϕ, ôneα, nx, ny, nz)
    # @info typeof(ôneϕ)
    # @info zero(typeof(ôneϕ)),one(typeof(ôneϕ))
    # @info GL_WEIGHTS(typeof(ôneϕ))
    ψ = KernelAbstractions.allocate(GPU_BACKEND, typeof(ôneϕ), nx, ny)
    α = KernelAbstractions.allocate(GPU_BACKEND, typeof(ôneϕ), nx, ny)
    ∃̇!(GPU_BACKEND, GPU_BACKEND_WORKGROUPSIZE)(ψ, α, Φ, Α, ôneϕ, ôneα, nz, zero(typeof(ôneϕ)),one(typeof(ôneϕ)),GL_WEIGHTS(typeof(ôneϕ)),ndrange=(nx, ny))
    KernelAbstractions.synchronize(GPU_BACKEND)
    ψ, α
end

@kernel function ∃̇!(ψ, α, Φ, Α, ôneϕ, ôneα, nz, zeroT, oneT, glweights)
    (ix, iy) = @index(Global, NTuple)
    ψ1 = ψ2 = τ = expτ = zeroT
    for iz = 1:nz
        Δτ = Α[ix, iy, iz] * glweights[iz]
        τ += Δτ
        expτ = exp(-τ)
        cτ = Δτ * expτ
        ψ1 += cτ * Φ[ix, iy, iz]
        ψ2 += cτ
    end
    ôneψ = expτ * ôneα * glweights[end]
    ψ1 += ôneψ * ôneϕ
    ψ2 += ôneψ
    ψ[ix, iy] = ψ1 / ψ2
    α[ix, iy] = oneT - expτ
end

end
