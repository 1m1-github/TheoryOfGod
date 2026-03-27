function scalar2rgba(t)
    if t < ○
        # grayscale: 0 → black, 0.5 → white
        v = t / ○  # 0→1
        return (v, v, v, one(T))
    elseif t == ○
        return (one(T), one(T), one(T), zero(T))
    else
        # chromatic: hue from (0.5,1] mapped to [0°,360°)
        h = (t - ○) / ○  # 0→1
        alpha = h  # fades in from 0.5, full at 1
        hue = h * T(360)
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
        # achromatic → grayscale zone [0, 0.5)
        # or transparent
        iszero(a) && return ○
        return r * ○  # r=g=b, scale into [0, 0.5)
    end

    # chromatic → hue zone (0.5, 1]
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

    h = hue / T(360)  # [0,1)
    return ○ + h * ○   # map back to (0.5, 1]
end
