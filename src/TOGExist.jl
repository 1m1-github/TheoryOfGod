module TOGExist

export ∃!, ∩, ∃, copy!

using KernelAbstractions, IntervalTrees, TOGGPU
import Base.:∩

abstract type ∀ end
struct ∃{F} <: ∀
    d::Vector{Float64}
    μ::Vector{Float64}
    ρ::Vector{Float64}
    ∂₀::Vector{Bool}
    ∂₁::Vector{Bool}
    ϕ::F
    """
    I = [ZERO < ○ < ONE] is a unit 1-dim space of information with origin ○ (no information) at its center including the corners ZERO and ONE.
    Ω = I^I is an ∞-dim normed vector space.
    We have an ontology 𝕋 on Ω such that ϵ ∈ 𝕋:
    * ϵ ⊆ Ω
    * ϵ ∩ ϵ′ ≠ ∅ => ϵ = ϵ′ <=> ϵ ≠ ϵ′ => ∃dim: ϵ|dim ∩ ϵ′|dim = ∅
    ϵ defines its existence inside a finite-dim box of non-zero support with origin, radius and closed xor open in borders each direction. All undeclared dims of ϵ have μ=○ and ρ=0.
    Entities from outside Ω that can create ϵ ⊊ Ω are called `god`s.
    god ⊊ GOD = Ω = I^I = I^(.) = [ZERO < ○ < ONE]^(.)
    god observes or creates, GOD iterates.
    The complexity Ο::BigInt counts the number of existences in Ω. Time is derived from complexity as t = log(Ο)/(1+log(Ο)).
    # Arguments
    - `d::AbstractVector{<:Number}`: Dimensions of the finite non-zero support of ϵ.
    - `μ::AbstractVector{<:Number}`: Center of ϵ.
    - `ρ::AbstractVector{<:Number}`: Radius of ϵ.
    - `∂₀::AbstractVector{Bool}`: Closed xor open left borders of ϵ.
    - `∂₁::AbstractVector{Bool}`: Closed xor open right borders of ϵ.
    - `ϕ`: Ω -> I is arbitrary, computable and smooth fuzzy existence potential.
    length(d) == length(μ) == length(ρ) == length(∂₀) == length(∂₁)
    """
    function ∃(d::AbstractVector{<:Number}, μ::AbstractVector{<:Number}, ρ::AbstractVector{<:Number}, ∂₀::AbstractVector{Bool}, ∂₁::AbstractVector{Bool}, ϕ::F) where {F}
        @assert length(d) == length(μ) == length(ρ) == length(∂₀) == length(∂₁)
        N = length(d)
        p = sortperm(d)
        ḋ = [d[p[i]] for i = eachindex(p)]
        μ̇ = [μ[p[i]] for i = eachindex(p)]
        ρ̇ = [ρ[p[i]] for i = eachindex(p)]
        ∂̇₀ = [∂₀[p[i]] for i = eachindex(p)]
        ∂̇₁ = [∂₁[p[i]] for i = eachindex(p)]
        containscenter = true
        z, o = μ̇ .- ρ̇, μ̇ .+ ρ̇
        @assert all(0 .≤ z)
        @assert all(o .≤ 1)
        @assert all(0 .≤ ḋ .≤ 1)
        for i = 1:N
            1 < i && @assert ḋ[i-1] ≠ ḋ[i]
            @assert !iszero(ρ̇[i]) || (∂̇₀[i] && ∂̇₁[i] && μ̇[i] ≠ ○)
            containscenter &= z[i] ≤ ○ ≤ o[i]
        end
        @assert !containscenter # ∃({}) ∈ Ω
        ϕ̂ = Φ(Φ̇(ϕ), ḋ, z, o)
        @assert iscomputable(ϕ, N)
        new{typeof(ϕ̂)}(ḋ, μ̇, ρ̇, ∂̇₀, ∂̇₁, ϕ̂)
    end
end
struct Φ{F}
    ϕ::F
    d::Vector{Float64}
    z::Vector{Float64}
    o::Vector{Float64}
    Φ(_ϕ, _d, _z, _o) = new{typeof(_ϕ)}(_ϕ, _d, _z, _o)
end
(ϕ::∀)(x) = ○
function (ϕ::Φ)(x)
    ẋ = similar(x)
    for i = eachindex(x)
        ρ = ϕ.o[i] - ϕ.z[i]
        ẋ[i] = iszero(ρ) ? ○ : (x[i] - ϕ.z[i]) / ρ
    end
    ϕ.ϕ(ẋ) # todo Base.invokelatest?
end
Φ̇(ϕ::Φ) = ϕ.ϕ
Φ̇(ϕ) = ϕ
function iscomputable(ϕ, N)
    try
        @kernel gpu(ϕ̃, x) = ϕ̃(x)
        x = KernelAbstractions.zeros(TOGGPU.GPU_BACKEND, Float64, N)
        gpu(TOGGPU.GPU_BACKEND, TOGGPU.GPU_BACKEND_WORKGROUPSIZE)(ϕ, x, ndrange=1)
        true
    catch e
        bt = catch_backtrace()
        showerror(stderr, e, bt)
        false
    end
end
struct 𝕋 <: ∀
    ϵ::Dict{Float64,IntervalMap{Float64,Set{∃}}}
    Ο::Dict{∀,BigInt}
    n::Dict{String,∀}
    L::ReentrantLock
    function 𝕋()
        ϵ = Dict{Float64,IntervalMap{Float64,Set{∃}}}()
        Ο = Dict{∀,BigInt}()
        n = Dict{String,∀}()
        L = ReentrantLock()
        Ω = new(ϵ, Ο, n, L)
        Ω.Ο[Ω] = BigInt(1)
        Ω.n["Ω"] = Ω
        Ω
    end
end
t(Ο::Integer) = 1 - 1 / (1 + log(Ο))
t(ϵ::∀, ω::𝕋) = t(ω.Ο[ϵ])
t(ω::𝕋) = t(ω, ω)
Ο(t::Number) = round(BigInt, exp(t / (1 - t)))
○ = 0.5
○̂ = x -> ○
Base.copy!(ϵ::∃, ḋ, μ̇, ρ̇, n, ω::𝕋) = ∃!(∃(ḋ, μ̇, ρ̇, ϵ.∂₀, ϵ.∂₁, ϵ.ϕ), n, ω)
function ∩ᵢ(z₁, o₁, ∂₁₀, ∂₁₁, z₂, o₂, ∂₂₀, ∂₂₁)
    ż = max(z₁, z₂)
    ȯ = min(o₁, o₂)
    ż < ȯ && return true
    ż ≠ ȯ && return false
    ∂₀₀ = z₂ < z₁ ? ∂₁₀ : (z₁ < z₂ ? ∂₂₀ : ∂₁₀ && ∂₂₀)
    ∂₀₁ = o₁ < o₂ ? ∂₁₁ : (o₂ < o₁ ? ∂₂₁ : ∂₁₁ && ∂₂₁)
    ∂₀₀ && ∂₀₁
end
function ∩ᵢ(ϵ::∃, ω::𝕋)
    β = Dict{Float64,Tuple{Union{Set{∃},Nothing},Vector{Function}}}()
    addonlyϵ(z, o, d) = ω.ϵ[d][(z, o)] = Set{∃}((ϵ,))
    addonlyϵ(x, d) = addonlyϵ(x, x, d)
    addintervalϵ(z, o, d) = push!(ω.ϵ[d][Interval(z, o)].value, ϵ)
    addintervalϵ(x, d) = addintervalϵ(x, x, d)
    for (i, d) = enumerate(ϵ.d)
        command = Function[]
        z, o = ϵ.μ[i] - ϵ.ρ[i], ϵ.μ[i] + ϵ.ρ[i]
        if !haskey(ω.ϵ, d)
            β̇ = z ≤ ○ ≤ o ? nothing : Set{∃}()
            push!(command, () -> begin
                ω.ϵ[d] = IntervalMap{Float64,Set{∃}}()
                addonlyϵ(z, o, d)
            end)
            if z ≠ o
                if ϵ.∂₀[i]
                    push!(command, () -> addonlyϵ(z, d))
                end
                if ϵ.∂₁[i]
                    push!(command, () -> addonlyϵ(o, d))
                end
            end
            β[d] = β̇, command
            continue
        end
        if haskey(ω.ϵ[d], (z, o))
            push!(command, () -> addintervalϵ(z, o, d))
        else
            push!(command, () -> addonlyϵ(z, o, d))
        end
        if z ≠ o
            if ϵ.∂₀[i]
                if haskey(ω.ϵ[d], (z, z))
                    push!(command, () -> addintervalϵ(z, d))
                else
                    push!(command, () -> addonlyϵ(z, d))
                end
            end
            if ϵ.∂₁[i]
                if haskey(ω.ϵ[d], (o, o))
                    push!(command, () -> addintervalϵ(o, d))
                else
                    push!(command, () -> addonlyϵ(o, d))
                end
            end
        end
        β̇s = map(p -> p.value, collect(intersect(ω.ϵ[d], (z, o))))
        β̇ = isempty(β̇s) ? (z ≤ ○ ≤ o ? nothing : Set{∃}()) : ∪(β̇s...)
        β[d] = β̇, command
    end
    β
end
function ∩ᵢ(ϵ::∃, β, ω::𝕋, t::Number)
    β̃ = filter(β̇ -> !isnothing(β̇[2][1]), β)
    β̇ = Set{∃}()
    isempty(β̃) && return β̇
    # Ο̇ = Ο(t)
    for ϵ̃ = ∩(map(t -> t[1], collect(values(β̃)))...)
        # Οϵ̃ = ω.Ο[ϵ̃]
        for i = eachindex(ϵ̃.d)
            # if Οϵ̃ ≤ Ο̇ && ∩ᵢ(ϵ̃.μ[i] - ϵ̃.ρ[i], ϵ̃.μ[i] + ϵ̃.ρ[i], ϵ̃.∂₀[i], ϵ̃.∂₁[i], ϵ.μ[i] - ϵ.ρ[i], ϵ.μ[i] + ϵ.ρ[i], ϵ.∂₀[i], ϵ.∂₁[i])
            if ∩ᵢ(ϵ̃.μ[i] - ϵ̃.ρ[i], ϵ̃.μ[i] + ϵ̃.ρ[i], ϵ̃.∂₀[i], ϵ̃.∂₁[i], ϵ.μ[i] - ϵ.ρ[i], ϵ.μ[i] + ϵ.ρ[i], ϵ.∂₀[i], ϵ.∂₁[i])
                push!(β̇, ϵ̃)
                break
            end
        end
    end
    β̇
end
∩(ϵ::∃, ω::𝕋, t::Number) = ∩(∩ᵢ(ϵ, ∩ᵢ(ϵ, ω), ω, t))
function ∃!(ϵ::∃, n::AbstractString, ω::𝕋)
    lock(ω.L) do
        β = ∩ᵢ(ϵ, ω)
        β̇ = ∩ᵢ(ϵ, β, ω, 0)
        if !isempty(β̇)
            @error "Intersection found."
            return
        end
        for (_, f) = values(β), ḟ = f ḟ() end
        # while Sys.free_memory() < Base.summarysize(ω) + Base.summarysize(ϵ) # todo rm oldest ϵ
        #     rm!(ω)
        # end
        ω.Ο[ω] += 1
        ω.Ο[ϵ] = ω.Ο[ω]
        ω.n[n] = ϵ
        ϵ
    end
end

end
