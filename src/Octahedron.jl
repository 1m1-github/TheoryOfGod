# z, o=SVector{3}(g.ẑero.μ[end-2:end]), SVector{3}(g.ône.μ[end-2:end])
function calc_ew(z, o)
    d = o - z
    u = SVector(iszero(d[1]) && iszero(d[2]) ? (zero(T), d[3], -d[2]) : (d[2], -d[1], zero(T)))
    u = u / norm(u)
    v = cross(d, u)
    v = v / norm(v)
    a = abs.(u)
    b = abs.(v)
    L = abs.(d)
    A = []
    # i = 1
    for i = 1:3
        zero(T) < a[i] || continue
        zero(T) < b[i] || continue
        wx = L[i] / a[i] / 4
        wy = L[i] / b[i] / 4
        a[i] * wx + b[i] * wy ≤ L[i] / 2 || continue
        j = (i + 1) % 3 + 1
        a[j] * wx + b[j] * wy ≤ L[j] / 2 || continue
        j = (j + 1) % 3 + 1
        a[j] * wx + b[j] * wy ≤ L[j] / 2 || continue
        push!(A, (wx, wy, wx * wy))
    end
    # i = 1
    # j = 2
    for i = 1:2, j = i+1:3
        D = a[i] * b[j] - a[j] * b[i]
        wx = (L[i] * b[j] - L[j] * b[i]) / D / 2
        wy = (L[j] * a[i] - L[i] * a[j]) / D / 2
        zero(T) < wx || continue
        zero(T) < wy || continue
        k = 6 - i - j
        a[k] * wx + b[k] * wy ≤ L[k] / 2 || continue
        push!(A, (wx, wy, wx * wy))
    end
    _, maxi = findmax(a -> a[3], A)
    wx = A[maxi][1]
    wy = A[maxi][2]
    ex = u * wx
    ey = v * wy
    @show wx,wy
    ex, ey, wx, wy
end

"""
Tighten an interval [lo, hi] given constraint L ≤ s*coeff + s_other*coeff_other ≤ U,
relaxing over s_other ∈ [-1,1]. Returns (lo, hi, empty).
"""
@inline function tighten(lo, hi, coeff, abs_other, L, U)
    if zero(T) < coeff
        lo = max(lo, (L - abs_other) / coeff)
        hi = min(hi, (U + abs_other) / coeff)
    elseif coeff < zero(T)
        lo = max(lo, (U + abs_other) / coeff)
        hi = min(hi, (L - abs_other) / coeff)
    else
        (L > abs_other || U < -abs_other) && return (lo, hi, true)
    end
    return (lo, hi, lo > hi)
end

"""
Convert continuous s ∈ [-1,1] bounds to 1-based grid indices on an axis of size ns.
s_i = (2i - 1 - ns) / (ns - 1)  ⟹  i = (s*(ns-1) + ns + 1) / 2
Returns (ilo, ihi) or (0, 0) if empty.
"""
@inline function s_to_indices(slo, shi, ns)
    h = (ns - 1) * ○
    mid = (ns + 1) * ○
    ilo = clamp(ceil(Int, slo * h + mid), 1, ns)
    ihi = clamp(floor(Int, shi * h + mid), 1, ns)
    ilo > ihi && return (0, 0)
    return (ilo, ihi)
end

"""
    pyramid_box_intersection(z, o, ex, ey, wx, wy, a, b, nx, ny, nz)

Per pyramid slice k=1:nz, find (i_min, i_max, j_min, j_max) index ranges of grid
points inside axis-aligned box [a,b]. Empty slices → (0,0,0,0). Cost: O(nz·n).
"""
# i=i
# ĩ=length(ΦΦ)
# z=g.ẑero.μ
# o=g.ône.μ
# ex=
# ey=
# wx=
# wy=
# a=ϵ.μ .- ϵ.ρ
# b=ϵ.μ .+ ϵ.ρ
# nx=g.♯[1]
# ny=g.♯[2]
# nz=GL_N
function pyramid_box_intersection(
    i, ĩ,
    z::SVector{N,T}, o::SVector{N,T},
    ex, ey, wx, wy, a, b, nx, ny, nz
) where {N}
    c = (z + o) * ○
    d = o - c

    intersects = false

    @inbounds for k = 1:nz
        t = k / (nz + 1)
        t = GL_NODES[k]
        omt = 1 - t
        wxk = wx * omt
        wyk = wy * omt

        si_lo, si_hi = -one(T), one(T)
        sj_lo, sj_hi = -one(T), one(T)
        empty = false

        for m = 1:N
            ck = c[m] + t * d[m]
            α = wxk * ex[m]
            β = wyk * ey[m]
            L = a[m] - ck
            U = b[m] - ck

            si_lo, si_hi, empty = tighten(si_lo, si_hi, α, abs(β), L, U)
            empty && break
            sj_lo, sj_hi, empty = tighten(sj_lo, sj_hi, β, abs(α), L, U)
            empty && break
        end
        empty && continue

        ilo, ihi = s_to_indices(si_lo, si_hi, nx)
        jlo, jhi = s_to_indices(sj_lo, sj_hi, ny)
        (iszero(ilo) || iszero(jlo)) && continue
        i[ilo:ihi, jlo:jhi, k] .= ĩ
        intersects = true
    end

    intersects
end

# wx = W
# wy = H
# nx, ny, nz = 2000, 2000, 8
# z = SVector(0.0, 0.0, 0.0)
# o = SVector(1.0, 1.0, 1.0)
# a = SVector(0.75, 0.75, 0.75)
# b = SVector(1.0, 1.0, 1.0)
# @time out = pyramid_box_intersection(z, o, ex, ey, wx, wy, a, b, nx, ny, nz)

