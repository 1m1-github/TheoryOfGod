"""
I = [ZERO < ○ < ONE] denotes a unit 1-dim space of information with origin ○ (no information) at its center including the corners ZERO and ONE.
Ω = I^I an ∞-dim normed and smooth vector space.
We have an ontology 𝕋 on Ω such that ϵ ∈ 𝕋:
* ϵ ⊆ Ω
* ϵ ∩ ϵ′ ≠ ∅ => ϵ = ϵ′
* ϵ.ϕ: Ω -> I is arbitrary, computable and smooth fuzzy existence potential towards ONE=true xor ZERO=false.
ϵ ⊊ Ω defines its existence inside a subset of Ω using an origin (μ), a radius (ρ) and a closed vs. open in each direction (∂) vector. These vectors are finite and all other dimensional coordinates of ϵ are ○.
god ⊊ GOD = Ω = I^I = I^(.) = [ZERO < ○ < ONE]^(.)
god observes or creates, GOD iterates.
"""
module TOG
using KernelAbstractions, IntervalTrees, TOGGPU
import Base.:∩
abstract type ∀ end
struct ∃{V<:AbstractVector,P<:AbstractVector{Bool},F} <: ∀
    d::V
    μ::V
    ρ::V
    ∂₀::P
    ∂₁::P
    ϕ::F
    function ∃(d::V, μ::V, ρ::V, ∂₀::P, ∂₁::P, ϕ::F) where {V<:AbstractVector,P<:AbstractVector{Bool},F}
        # @show length(d), length(μ), length(ρ), length(∂₀), length(∂₁)
        # @info typeof(d), typeof(μ), typeof(ρ), typeof(∂₀), typeof(∂₁)
        # @info length(d) ,length(μ) ,length(ρ) ,length(∂₀) ,length(∂₁)
        @assert length(d) == length(μ) == length(ρ) == length(∂₀) == length(∂₁)
        @assert eltype(μ) == eltype(ρ) == eltype(d) # todo needed?
        # @info "A"
        T = eltype(μ)
        # @info T
        N = length(d)
        # @info N
        p = sortperm(d)
        # @info p
        ḋ = [d[p[i]] for i = eachindex(p)]
        μ̇ = [μ[p[i]] for i = eachindex(p)]
        ρ̇ = [ρ[p[i]] for i = eachindex(p)]
        ∂̇₀ = [∂₀[p[i]] for i = eachindex(p)]
        ∂̇₁ = [∂₁[p[i]] for i = eachindex(p)]
        containscenter = true
        z, o = μ̇ .- ρ̇, μ̇ .+ ρ̇
        @assert all(zero(T) .≤ z)
        @assert all(o .≤ one(T))
        @assert all(zero(T) .≤ ḋ .≤ one(T))
        for i = 1:N
            # @info i
            # @show "TOG ∃", i, ḋ[i], μ̇[i], ρ̇[i], ∂̇₀[i], ∂̇₁[i], z, o
            1 < i && @assert ḋ[i-1] ≠ ḋ[i]
            @assert !iszero(ρ̇[i]) || (∂̇₀[i] && ∂̇₁[i] && μ̇[i] ≠ ○(T))
            containscenter &= z[i] ≤ ○(T) ≤ o[i]
        end
        @assert !containscenter # ∃({}) ∈ Ω
        # @info "!containscenter"
        ϕ̂ = Φ(Φ̇(ϕ), ḋ, z, o)
        # assert Φ̇(Φ)(zeros(T, N)) isa T # 
        @assert iscomputable(ϕ, T, N)
        # @info "iscomputable"
        new{V,P,typeof(ϕ̂)}(ḋ, μ̇, ρ̇, ∂̇₀, ∂̇₁, ϕ̂)
    end
end
function iscomputable(ϕ, T, N)
    try
        @kernel gpu(ϕ̃, x) = ϕ̃(x)
        x = KernelAbstractions.zeros(TOGGPU.GPU_BACKEND, T, N)
        gpu(TOGGPU.GPU_BACKEND, TOGGPU.GPU_BACKEND_WORKGROUPSIZE)(ϕ, x, ndrange=1)
        true
    catch e
        bt = catch_backtrace()
        showerror(stderr, e, bt)
        false
    end
end
struct 𝕋{T} <: ∀
    ϵ::Dict{T,IntervalMap{T,Set{∃}}}
    Ο::Dict{∀,UInt}
    s::Ref{UInt}
    n::Dict{String,∀}
    L::ReentrantLock
    function 𝕋(T)
        ϵ = Dict{T,IntervalMap{T,Set{∃}}}()
        Ο = Dict{∀,UInt}()
        s = Ref(UInt(1))
        n = Dict{String,∀}()
        L = ReentrantLock()
        Ω = new{T}(ϵ, Ο, s, n, L)
        Ω.Ο[Ω] = Ω.s[]
        Ω.n["Ω"] = Ω
        Ω
    end
end
t(Ο::UInt, T) = one(T) - one(T) / (one(T) + T(log(Ο)))
t(ϵ::∀, ω::𝕋) = t(ω.Ο[ϵ], T(ω))
t(ω::𝕋) = t(ω, ω)
Ο(t) = round(UInt, exp(t / (one(t) - t)))
○(T::DataType) = one(T) / (one(T) + one(T))
○(T) = ○(eltype(T))
○̂(T) = x -> ○(T)
T(ω::𝕋) = first(typeof(ω).parameters)
struct Φ{D<:AbstractVector,F}
    ϕ::F
    d::D
    z::D
    o::D
end
(ϕ::∀)(x) = ○(x)
function (ϕ::Φ)(x)
    # ○̇ = ○(x)
    # @show "ϕ", x
    # for i = 1:length(ϕ.d)
    #     ϕ.o[i] == ϕ.z[i] && return ○̇
    #     x[i] ≤ ϕ.z[i] && return ○̇
    #     ϕ.o[i] ≤ x[i] && return ○̇
    # end
    # Base.invokelatest(ϕ.Φ, x)
    # @show "ϕ", ϕ.Φ(x)
    ẋ = (x .- ϕ.z) ./ (ϕ.o .- ϕ.z)
    ϕ.ϕ(ẋ)
    # ϕ.Φ(x)
end
Φ̇(ϕ::Φ) = ϕ.ϕ
Φ̇(ϕ) = ϕ
Base.copy!(ϵ::∃, ḋ, μ̇, ρ̇, n, ω::𝕋) = ∃!(∃(ḋ, μ̇, ρ̇, ϵ.∂₀, ϵ.∂₁, ϵ.ϕ), n, ω)
function ∩(z₁, o₁, ∂₁₀, ∂₁₁, z₂, o₂, ∂₂₀, ∂₂₁)
    ż = max(z₁, z₂)
    ȯ = min(o₁, o₂)
    ż < ȯ && return true
    ż ≠ ȯ && return false
    ∂₀₀ = z₂ < z₁ ? ∂₁₀ : (z₁ < z₂ ? ∂₂₀ : ∂₁₀ && ∂₂₀)
    ∂₀₁ = o₁ < o₂ ? ∂₁₁ : (o₂ < o₁ ? ∂₂₁ : ∂₁₁ && ∂₂₁)
    ∂₀₀ && ∂₀₁
end
function ∩ᵢ(ϵ::∃, ω::𝕋)
    β = Dict{T(ω),Tuple{Union{Set{∃},Nothing},Vector{Function}}}()
    addonlyϵ(z, o, d) = ω.ϵ[d][(z, o)] = Set{∃}((ϵ,))
    addonlyϵ(x, d) = addonlyϵ(x, x, d)
    addintervalϵ(z, o, d) = push!(ω.ϵ[d][Interval(z, o)].value, ϵ)
    addintervalϵ(x, d) = addintervalϵ(x, x, d)
    # Threads.@threads
    for (i, d) = enumerate(ϵ.d)
        command = Function[]
        z, o = ϵ.μ[i] - ϵ.ρ[i], ϵ.μ[i] + ϵ.ρ[i]
        if !haskey(ω.ϵ, d)
            β̇ = z ≤ ○(T(ω)) ≤ o ? nothing : Set{∃}()
            push!(command, () -> begin
                ω.ϵ[d] = IntervalMap{T(ω),Set{∃}}()
                addonlyϵ(z, o, d)
                # ω.ϵ[d][(z, o)] = Set{∃}((ϵ,))
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
        β̇ = isempty(β̇s) ? (z ≤ ○(T(ω)) ≤ o ? nothing : Set{∃}()) : ∪(β̇s...)
        β[d] = β̇, command
    end
    β
end
function ∩(ϵ::∃, β, Ο, ω::𝕋)
    β̃ = filter(β̇ -> !isnothing(β̇[2][1]), β)
    β̇ = Set{∃}()
    isempty(β̃) && return β̇
    for ϵ̃ = ∩(map(t -> t[1], collect(values(β̃)))...)
        for i = eachindex(ϵ̃.d)
            if Ο ≤ ω.Ο[ϵ̃] && ∩(ϵ̃.μ[i] - ϵ̃.ρ[i], ϵ̃.μ[i] + ϵ̃.ρ[i], ϵ̃.∂₀[i], ϵ̃.∂₁[i], ϵ.μ[i] - ϵ.ρ[i], ϵ.μ[i] + ϵ.ρ[i], ϵ.∂₀[i], ϵ.∂₁[i])
                push!(β̇, ϵ̃)
            end
        end
    end
    β̇
end
∩(ϵ::∃, ω::𝕋, Ο) = ∩(∩(ϵ, ∩ᵢ(ϵ, ω), Ο, ω))
function ∃!(ϵ::∃, n, ω::𝕋)
    lock(ω.L)
    β = ∩ᵢ(ϵ, ω)
    # @info typeof(ϵ), typeof(β)
    β̇ = ∩(ϵ, β, one(UInt), ω)
    isempty(β̇) || (unlock(ω.L); error("Intersection found."))
    for (_, f) = values(β), ḟ = f
        ḟ()
    end
    # @info "∃!", ϵ.μ, ϵ.ρ # DEBUG
    # while Sys.free_memory() < ω.s[] + sizeof(ϵ)
    #     rm!(ω)
    # end
    ω.s[] += sizeof(ϵ)
    ω.Ο[ω] += 1
    ω.Ο[ϵ] = ω.Ο[ω]
    ω.n[n] = ϵ
    unlock(ω.L)
    ϵ
end
# function ∃̇(ϵ::∃, ω::𝕋) # todo parallel global running β̇ for speed
#     β = Dict{T(ω),Set{∃}}()
#     # Threads.@threads 
#     for (i, d) = enumerate(ϵ.d) # todo sorted by size of dim for speed
#         haskey(ω.ϵ, d) || return ○(T(ω))
#         possible𝕋1 = collect(intersect(ω.ϵ[d], (ϵ.μ[i], ϵ.μ[i])))
#         filter!(p -> p.first == p.last || p.first < ϵ.μ[i] < p.last, possible𝕋1)
#         β̇ = ∪(map(p -> p.value, possible𝕋1)...)
#         # @show "∃̇", 1, length(β̇) # DEBUG
#         ϵ.μ[i] ≠ ○(T(ω)) && isempty(β̇) && return ○(T(ω))
#         β[d] = β̇
#     end
#     β̇ = Set{∃}()
#     for ϵ̃ = ∩(values(β)...)
#         isgood = true
#         for (i, d) = enumerate(ϵ̃.d)
#             d ∈ ϵ.d && continue
#             if !∩(ϵ̃.μ[i] - ϵ̃.ρ[i], ϵ̃.μ[i] + ϵ̃.ρ[i], ϵ̃.∂₀[i], ϵ̃.∂₁[i], ○(T(ω)), ○(T(ω)), true, true)
#                 # push!(β̇, ϵ̃)
#                 isgood = false
#             end
#         end
#         isgood && push!(β̇, ϵ̃)
#     end
#     # @show "∃̇", 2, length(β̇) # DEBUG
#     isempty(β̇) && return ○(T(ω))
#     only(β̇).ϕ(ϵ.μ)
# end
end
