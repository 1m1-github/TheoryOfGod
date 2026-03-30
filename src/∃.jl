abstract type ∀ end
struct ∃{N,F,P<:∀} <: ∀
    ϵ̂::P
    d::SVector{N,T}
    μ::SVector{N,T}
    ρ::SVector{N,T}
    ∂::SVector{N,Tuple{Bool,Bool}}
    Φ::F
    h::UInt
    function ∃(ϵ̂::∀, d::SVector{N,T}, μ::SVector{N,T}, ρ::SVector{N,T}, ∂::SVector{N,Tuple{Bool,Bool}}, Φ::F) where {N,F}
        p = sortperm(d)
        ḋ, μ̇, ρ̇, ∂̇ = map(x -> x[p], (d, μ, ρ, ∂))
        for i = 1:N
            @assert zero(T) ≤ ḋ[i] ≤ one(T)
            1 < i && @assert ḋ[i-1] ≠ ḋ[i]
            @assert zero(T) ≤ μ̇[i] - ρ̇[i] ≤ μ̇[i] + ρ̇[i] ≤ one(T)
        end
        μΩ, ρΩ = μρΩ(ϵ̂, μ̇, ρ̇)
        ϕ = Φ̂(Φ̇(Φ), μΩ .- ρΩ, μΩ .+ ρΩ)
        @assert gpu_safe(ϕ, Val(N))
        h = hash(ḋ, hash(μ̇, hash(ρ̇, hash(∂̇, hash(ϵ̂)))))
        new{N,typeof(ϕ),typeof(ϵ̂)}(ϵ̂, ḋ, μ̇, ρ̇, ∂̇, ϕ, h)
    end
end
Base.hash(ϵ::∃, h::UInt) = hash(ϵ.h, h)
struct 𝕋 <: ∀
    ϵ̃::Dict{∀,Vector{∃}}
    Ο::Dict{∀,UInt}
    L::ReentrantLock
    s::Ref{UInt}
    function 𝕋()
        ϵ̃ = Dict{∀,Vector{∃}}()
        Ο = Dict{∀,UInt}()
        Ω = new(ϵ̃, Ο, ReentrantLock(), Ref(UInt(1)))
        Ω.ϵ̃[Ω] = ∃[]
        Ω.Ο[Ω] = Ω.s[]
        Ω
    end
end
Base.hash(::𝕋, h::UInt) = hash(:Ω, h)
t(Ο::UInt) = one(T) - one(T) / (one(T) + T(log(Ο)))
t(ϵ::∀, ω::𝕋) = t(ω.Ο[ϵ])
t(ω::𝕋) = t(ω, ω)
const ○ = one(T) / (one(T) + one(T))
const ○̂ = x -> ○
struct Φ̂{N,F}
    Φ::F
    z::SVector{N,T}
    o::SVector{N,T}
end
function (ϕ::Φ̂{N})(x) where N
    for i = 1:N
        ϕ.o[i] == ϕ.z[i] && return ○
        x[i] ≤ ϕ.z[i] && return ○
        ϕ.o[i] ≤ x[i] && return ○
    end
    ẋ = (x .- ϕ.z) ./ (ϕ.o .- ϕ.z)
    ϕ.Φ(ẋ)
end
Φ̇(Φ::Φ̂) = Φ.Φ
Φ̇(Φ) = Φ
function gpu_safe(Φ, ::Val{N}) where N
    @kernel function gpu_test(Φ, y, ::Val{N}) where N
        x = @SVector zeros(T, N)
        y[1] = Φ(x)
    end
    try
        y = KernelAbstractions.zeros(GPU_BACKEND, T, 1)
        gpu_test(GPU_BACKEND, GPU_BACKEND_WORKGROUPSIZE)(Φ, y, Val(N), ndrange=1)
        KernelAbstractions.synchronize(GPU_BACKEND)
        true
    catch e
        showerror(stderr, e, catch_backtrace())
        false
    end
end
function Base.copy!(ϵ::∃, ḋ, μ̇, ρ̇, ω::𝕋)
    ϵ̇ = ∃!(ḋ, μ̇, ρ̇, ϵ.∂, ϵ.Φ, ω)
    for ϵ̃ = ω.ϵ̃[ϵ]
        μ̃ = μ̂̂(μ̇, ρ̇, ϵ̃.μ)
        ρ̃ = ρ̂̂(ρ̇, ϵ̃.ρ)
        copy!(ϵ̃, ḋ, μ̃, ρ̃, ω)
    end
    ϵ̇
end
ρ̂̂(ρ̂, ρ) = T(2) .* ρ̂ .* ρ
μ̂̂(μ̂, ρ̂, μ) = μ̂ .+ ρ̂ .* (T(2) .* μ .- one(T))
ρ̃̃(ρ, ρ̃) = ○ * ρ ./ ρ̃
μ̃̃(μ, ρ, μ̃) = ○ * ((μ .- μ̃) ./ ρ .+ one(T))
function μρΩ(ϵ::∃)
    ϵ.ϵ̂ isa 𝕋 && return ϵ.μ, ϵ.ρ
    μ̂, ρ̂ = μρΩ(ϵ.ϵ̂)
    μ̂̂(μ̂, ρ̂, ϵ.μ), ρ̂̂(ρ̂, ϵ.ρ)
end
μρΩ(::𝕋, μ, ρ) = μ, ρ
function μρΩ(ϵ::∃, μ, ρ)
    μ̇, ρ̇ = μρΩ(ϵ)
    μ̂̂(μ̇, ρ̇, μ), ρ̂̂(ρ̇, ρ)
end
function μρ(ϵd, ϵμ, ϵρ, ϵ∂, d)
    i = searchsortedfirst(ϵd, d)
    N = length(ϵd)
    i ≤ N && ϵd[i] == d && return ϵμ[i], ϵρ[i], ϵ∂[i]
    d₀ = d₁ = ϵd[1]
    dₙ = ϵd[N]
    μ₀ = μ₁ = ○
    if d < d₀
        d₀ = zero(T)
        μ₁ = ϵμ[1]
    elseif dₙ < d
        d₀, d₁ = dₙ, one(T)
        μ₀ = ϵμ[N]
    else
        d₀, d₁ = ϵd[i-1], ϵd[i]
        μ₀, μ₁ = ϵμ[i-1], ϵμ[i]
    end
    d = (d - d₀) / (d₁ - d₀)
    μ₀ + (μ₁ - μ₀) * d, zero(T), (true, true)
end
function ⊂(z, o, ∂, ẑ, ô, ∂̂)
    # ż = ẑ < z || (ẑ == z && (!∂[1] && ∂̂[1]))
    # ȯ = o < ô || (o == ô && (!∂[1] && ∂̂[1]))
    ż = ẑ < z || (ẑ == z && (!∂[1] || ∂̂[1]))
    ȯ = o < ô || (o == ô && (!∂[1] || ∂̂[1]))
    ż && ȯ
end
function ⫉(ϵd, ϵμ, ϵρ, ϵ∂, ϵ̂d, ϵ̂μ, ϵ̂ρ, ϵ̂∂)
    for (î, d̂) in enumerate(ϵ̂d)
        ρ̂ = ϵ̂ρ[î]
        μ̂ = ϵ̂μ[î]
        μ, ρ, ∂ = μρ(ϵd, ϵμ, ϵρ, ϵ∂, d̂)
        z, o = μ - ρ, μ + ρ
        ẑ, ô = μ̂ - ρ̂, μ̂ + ρ̂
        !⊂(z, o, ∂, ẑ, ô, ϵ̂∂[î]) && return false
    end
    true
end
β(d, μ, ρ, ∂, ω::𝕋) = β(d, μ, ρ, ∂, ω, ω)
function β(d, μ, ρ, ∂, ϵ̂::∀, ω::𝕋)
    ϵ̃ = filter(ϵ -> begin
        ϵμ, ϵρ = μρΩ(ϵ)
        ⫉(d, μ, ρ, ∂, ϵ.d, ϵμ, ϵρ, ϵ.∂)
    end, ω.ϵ̃[ϵ̂])
    isempty(ϵ̃) && return ϵ̂
    β(d, μ, ρ, ∂, only(ϵ̃), ω)
end
function Base.:∩(z₁, o₁, ∂₁, z₂, o₂, ∂₂)
    ż = max(z₁, z₂)
    ȯ = min(o₁, o₂)
    ż < ȯ && return true
    ż ≠ ȯ && return false
    ∂₀₀ = z₂ < z₁ ? ∂₁[1] : (z₁ < z₂ ? ∂₂[1] : ∂₁[1] && ∂₂[1])
    ∂₀₁ = o₁ < o₂ ? ∂₁[2] : (o₂ < o₁ ? ∂₂[2] : ∂₁[2] && ∂₂[2])
    ∂₀₀ && ∂₀₁
end
function Base.:∩(ϵ₁d, ϵ₁μ, ϵ₁ρ, ϵ₁∂, ϵ₂d, ϵ₂μ, ϵ₂ρ, ϵ₂∂, ϵ₂ϵ̂d)
    for d = ϵ₂ϵ̂d
        μ₁, ρ₁, ∂₁ = μρ(ϵ₁d, ϵ₁μ, ϵ₁ρ, ϵ₁∂, d)
        μ₂, ρ₂, ∂₂ = μρ(ϵ₂d, ϵ₂μ, ϵ₂ρ, ϵ₂∂, d)
        z₁, o₁ = μ₁ - ρ₁, μ₁ + ρ₁
        z₂, o₂ = μ₂ - ρ₂, μ₂ + ρ₂
        !∩(z₁, o₁, ∂₁, z₂, o₂, ∂₂) && return false
    end
    true
end
function ∃!(d, μ, ρ, ∂, ϕ, ω::𝕋)
    lock(ω.L)
    ϵ̂ = β(d, μ, ρ, ∂, ω)
    any(ϵ̃ -> begin
        ϵ̃ϵ̂d = ϵ̂ isa ∃ ? ϵ̂.d : unique(sort(d ∪ ϵ̃.d))
        ∩(d, μ, ρ, ∂, ϵ̃.d, ϵ̃.μ, ϵ̃.ρ, ϵ̃.∂, ϵ̃ϵ̂d)
    end, ω.ϵ̃[ϵ̂]) && (unlock(ω.L); return nothing)
    μ̃, ρ̃ = if ϵ̂ === ω
        μ, ρ
    else
        ϵ̂μ, ϵ̂ρ = μρΩ(ϵ̂)
        μ̃̃(ϵ̂μ, ϵ̂ρ, μ), ρ̃̃(ρ, ϵ̂ρ)
    end
    ϵ = ∃(ϵ̂, d, μ̃, ρ̃, ∂, ϕ)
    while Sys.free_memory() < ω.s[] + sizeof(ϵ)
        rm!(ω)
    end
    ω.s[] += sizeof(ϵ)
    ω.Ο[ω] += 1
    ω.Ο[ϵ] = ω.Ο[ω]
    push!(ω.ϵ̃[ϵ̂], ϵ)
    ω.ϵ̃[ϵ] = ∃[]
    unlock(ω.L)
    ϵ
end
function rm!(ω::𝕋)
    # todo
end
