module TOGColor

export rgba2scalar, scalar2rgba

using Colors
using TOGExist: ○

const HUE = 360

"""
Convert Colors.RGBA to a scalar value in [0,1].
example: `rgba2scalar(RGBA(0.4,0.3,0.5,1.0))`
"""
function rgba2scalar(rgba)
    iszero(rgba.r) && iszero(rgba.g) && iszero(rgba.b) && return rgba.r
    isone(rgba.r) && isone(rgba.g) && isone(rgba.b) && return rgba.r
    HSVA(rgba).h / HUE
end

"""
Convert a scalar value in [0,1] to Colors.RGBA.
example: `scalar2rgba(scalarvalue, alphavalue)`
"""
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
