module TOGRay

using FastGaussQuadrature, KernelAbstractions
using TOGGPU: GPU_BACKEND, GPU_BACKEND_WORKGROUPSIZE
using TOGExist: ○
const GL_N = 8
const GL_NODES_RAW, GL_WEIGHTS_RAW = gausslegendre(GL_N)
const GL_NODES(T) = ○ .+ GL_NODES_RAW ./ (GL_NODES_RAW[end] - GL_NODES_RAW[1]) # [0,1]
const GL_WEIGHTS(T) = ○ * GL_WEIGHTS_RAW

function ∃̇(Φ, Α, ôneϕ, ôneα, nx, ny, nz)
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
