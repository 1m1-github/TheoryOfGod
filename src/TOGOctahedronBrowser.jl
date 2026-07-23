module TOGOctahedronBrowser

using TOGBroadcastBrowser: BroadcastBrowser
using TOGOctahedron: Octahedron, take!
using TOGMoveOctahedron: scaleup!, scaledown!, moveup!, movedown!, focusup!, focusdown!, jerkup!, jerkdown!, rotateup!, rotatedown!, step!, δN
using TOGExist: ○
using TOGColor: scalar2rgba
using LoopOS: @whiletrue

const JSKEYPRESS = """document.addEventListener('keydown', (e) => {fetch('/keypress', {method: 'POST',body: e.key})})"""

function awaken(; octahedron, browser)
    @info "TOGOctahedronBrowser, awaken"
    octahedron.♯ = (Int(browser.width), Int(browser.height))
    BROWSER[] =
        Browser(
            octahedron,
            browserlooptask(octahedron),
            browser)
end
mutable struct Browser
    o::Octahedron
    loop::Union{Task,Nothing}
    browser::Union{BroadcastBrowser,Nothing}
end
const BROWSER = Ref{Browser}()
const OBSERVE = Ref(true)
browserlooptask(octahedron) = errormonitor(Threads.@spawn begin
    @info "TOGOctahedronBrowser, browserlooptask"
    # t = time()
    put!(BroadcastBrowser, JS(octahedron.♯[1], octahedron.♯[2]))
    put!(BroadcastBrowser, JSKEYPRESS)
    # # put!(browser.processor, JS(o.♯[1], o.♯[2]))
    T = first(typeof(octahedron).parameters)
    ϕ = fill(○(T), octahedron.♯[1], octahedron.♯[2])
    α = zeros(T, octahedron.♯[1], octahedron.♯[2])
    @whiletrue begin
                try
        #     #         # t̃ = time()
        #     #         # dt = t̃ - t
        #     #         # t = t̃
        #     #         # step!(o)
        # sleep(1) # DEBUG
        OBSERVE[] || continue
        ϕ̇, α̇ = Base.invokelatest() do
            take!(octahedron)
        end
        @info "unique(ϕ̇)", unique(ϕ̇) # DEBUG
        δ = Δ!(ϕ, α, ϕ̇, α̇)
        @info "length(δ)", length(δ)
        isempty(δ) && continue
        js = "pixel=" * writeδ(δ, octahedron.♯[2]) * "\n" * SET_PIXELS_JS
        # @info "length(js)", length(js)
        put!(BroadcastBrowser, js)
                catch e
                    bt = catch_backtrace()
                    showerror(stderr, e, bt)
                    sleep(1)
                end
                # break # DEBUG
    end
end)

function Δ!(ϕ, α, ϕ̇, α̇)
    T = eltype(ϕ)
    δ = Tuple{CartesianIndex{2},Tuple{T,T,T,T}}[]
    for i = CartesianIndices(ϕ̇)
        ϕ[i] == ϕ̇[i] && α[i] == α̇[i] && continue
        ϕ[i] = ϕ̇[i]
        α[i] = α̇[i]
        rgba = scalar2rgba(ϕ̇[i], α̇[i])
        if rgba.r<1.0 || rgba.g<1.0 || rgba.b<1.0
            @info "Δ!, non-white", i, rgba
        end
        push!(δ, (i, (T(rgba.r), T(rgba.g), T(rgba.b), T(rgba.alpha))))
    end
    δ
end
function writeδ(δ, height)
    result = []
    for (i, color) = δ
        push!(result, (i[1] - 1, height - 1 - (i[2] - 1), round.(UInt8, typemax(UInt8) .* color)...))
    end
    bracket(x) = "[" * x * "]"
    bracket(join(map(r -> bracket(join(r, ',')), result), ','))
end
const JS(width, height) = """
document.body.style.margin = '0'
document.body.style.display = 'flex'
document.body.style.justifyContent = 'center'
document.body.style.alignItems = 'center'
document.body.style.minHeight = '100vh'
canvas = document.createElement('canvas')
canvas.width = $(width)
canvas.height = $(height)
document.body.appendChild(canvas)
ctx = canvas.getContext('2d')
imageData = ctx.createImageData(canvas.width, canvas.height)
setPixel = (x, y, r, g, b, a) => {
    let i = (y * canvas.width + x) * 4
    imageData.data[i] = r
    imageData.data[i+1] = g
    imageData.data[i+2] = b
    imageData.data[i+3] = a
}
"""
const SET_PIXELS_JS = """
for (let [x,y,r,g,b,a] of pixel) setPixel(x,y,r,g,b,a)
ctx.putImageData(imageData, 0, 0)
"""

const CHANGE_MODE = Ref(2) # 0=ρ, 1=zero, 2=focus, 3=zero+focus
const CHANGE_DIM_INDEX = Ref(2)
function keypress(key)
    # try # DEBUG
    o = BROWSER[].o
    if key == "ArrowUp"
        if CHANGE_MODE[] == 0
            scaleup!(o, CHANGE_DIM_INDEX[])
        elseif CHANGE_MODE[] == 1
            moveup!(o, CHANGE_DIM_INDEX[])
        elseif CHANGE_MODE[] == 2
            focusup!(o, CHANGE_DIM_INDEX[])
        elseif CHANGE_MODE[] == 3
            moveup!(o, CHANGE_DIM_INDEX[])
            focusup!(o, CHANGE_DIM_INDEX[])
        end
    elseif key == "ArrowDown"
        if CHANGE_MODE[] == 0
            scaledown!(o, CHANGE_DIM_INDEX[])
        elseif CHANGE_MODE[] == 1
            movedown!(o, CHANGE_DIM_INDEX[])
        elseif CHANGE_MODE[] == 2
            focusdown!(o, CHANGE_DIM_INDEX[])
        elseif CHANGE_MODE[] == 3
            movedown!(o, CHANGE_DIM_INDEX[])
            focusdown!(o, CHANGE_DIM_INDEX[])
        end
    elseif key == "0"
        global CHANGE_MODE[] = (CHANGE_MODE[] + 1) % 4
    elseif key == "["
        δN += 1
    elseif key == "]"
        δN -= 1
    elseif key == "q"
        jerkup!(o)
    elseif key == "w"
        jerkdown!(o)
    elseif key == "d"
        rotateup!(o)
    elseif key == "f"
        rotatedown!(o)
    # elseif key == " "
        # put!(TOG)
    # elseif key == "t"
        # step!(o)
    elseif key == "Backspace"
        o.∂Ο = !o.∂Ο
    elseif key ∈ ["1", "2", "3", "4", "5", "6", "7", "8", "9"]
        global CHANGE_DIM_INDEX[] = parse(UInt, key)
    else
        return
    end
    println("key=$key")
    println("CHANGE_MODE=$CHANGE_MODE[]")
    println("CHANGE_DIM_INDEX=$CHANGE_DIM_INDEX[]")
    println("o.Ο=$(o.Ο)")
    println("o.∂Ο=$(o.∂Ο)")
    println("o.vΟ=$(o.vΟ)")
    println("o.v=$(o.v)")
    println("ẑeroμ=$(o.ẑeroμ)")
    println("ôneμ=$(o.ôneμ)")
    println("o.ρ=$(o.ρ)")
    println("o.θ=$(o.θ)")
    # println("δN=$δN")
    println("o.norm(o.ôneμ.-o.ẑeroμ)=$(o.norm(o.ôneμ.-o.ẑeroμ))")
# catch e 
#     # @info e 
#     bt = catch_backtrace()
#     showerror(stdout, e, bt)
# end
end

end
