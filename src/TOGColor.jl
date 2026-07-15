module TOGColor

export rgba2scalar

using Colors
using TOG: ○

const HUE = 360

function rgba2scalar(rgba)
    iszero(rgba.r) && iszero(rgba.g) && iszero(rgba.b) && return rgba.r
    isone(rgba.r) && isone(rgba.g) && isone(rgba.b) && return rgba.r
    HSVA(rgba).h / HUE
end
function scalar2rgba(ϕ, α)
    # @info ϕ, typeof(ϕ)
    # @info ○(typeof(ϕ))
    # @info ϕ == ○(typeof(ϕ))
    ϕ == ○(typeof(ϕ)) && return RGBA(1, 1, 1, α)
    iszero(ϕ) && return RGBA(0, 0, 0, α)
    isone(ϕ) && return RGBA(1, 1, 1, α)
    RGBA(HSVA(ϕ * HUE, 1, 1, α))
end

end
