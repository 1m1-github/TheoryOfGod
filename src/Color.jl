function scalar2rgba(t)
    if t < ○
        v = t / ○
        return (v, v, v, one(T))
    elseif t == ○
        return (one(T), one(T), one(T), zero(T))
    else
        h = (t - ○) / ○  # ∈ (0, 1]
        # interleave: odd decimal places → hue, even → alpha
        # but float precision kills this. pragmatic version:
        # use upper 8 bits for hue, lower 8 bits for alpha
        h16 = floor(Int, h * T(65536))
        hue_i = (h16 >> 8) & 0xFF
        alp_i = h16 & 0xFF

        hue = T(hue_i) / T(255) * T(360)
        alpha = T(alp_i) / T(255)

        x = one(T) - abs(mod(hue / T(60), T(2)) - one(T))
        r, g, b = if hue < T(60);       (one(T), x, zero(T))
        elseif hue < T(120); (x, one(T), zero(T))
        elseif hue < T(180); (zero(T), one(T), x)
        elseif hue < T(240); (zero(T), x, one(T))
        elseif hue < T(300); (x, zero(T), one(T))
        else;                 (one(T), zero(T), x)
        end
        return (r, g, b, alpha)
    end
end

rgba2scalar(rgba::PNGFiles.ColorTypes.RGBA) = rgba2scalar(rgba.r, rgba.g, rgba.b, rgba.alpha)
function rgba2scalar(r, g, b, a)
    cmax = max(r, g, b)
    cmin = min(r, g, b)
    delta = cmax - cmin

    if iszero(delta)
        iszero(a) && return ○
        return r * ○
    end

    hue = if cmax == r
        T(60) * mod((g - b) / delta, T(6))
    elseif cmax == g
        T(60) * ((b - r) / delta + T(2))
    else
        T(60) * ((r - g) / delta + T(4))
    end
    if hue < zero(T)
        hue += T(360)
    end

    hue_i = unsafe_trunc(Int, hue / T(360) * T(255))
    alp_i = unsafe_trunc(Int, a * T(255))
    h16 = (hue_i << 8) | alp_i
    h = T(h16) / T(65536)

    return ○ + h * ○
end

# rgba2scalar(1.0,0.0,0.0,0.3)
# scalar2rgba(rgba2scalar(1.0,0.0,0.0,0.3))
