struct godBrowser
    g::god
    loop::Task
    browser::BroadcastBrowser
end
godbrowser(g, browser) =
    begin
        looptask = Threads.@spawn begin
            t = time()
            put!(browser.processor, JS(g.♯[1], g.♯[2]))
            ϕ = zeros(T, g.♯[1], g.♯[2])
            while true
                try
                    yield()
                    t̃ = time()
                    dt = t̃ - t
                    t = t̃
                    step!(g, dt)
                    # sleep(1) # todo rm
                    ϕ̇ = Base.invokelatest() do
                        ∃̇(g, Ω[])
                    end
                    δ = Δ!(ϕ, ϕ̇)
                    isempty(δ) && continue
                    js = "pixel=" * writeδ(δ, g.♯[2]) * "\n" * SET_PIXELS_JS
                    put!(browser.processor, js)
                catch e
                    bt = catch_backtrace()
                    showerror(stderr, e, bt)
                    sleep(1)
                end
            end
        end
        godBrowser(g, looptask, browser)
    end
function godbrowserstart(browser)
    g = god(
        t=zero(T),
        d=sort(SA[invϕ, invϕ^2, one(T)]), # t, x, y, z
        ẑeroμ=SA[○-T(0.0), ○-T(0.0), ○],
        ôneμ=SA[○+T(0.0), ○+T(0.0), ○+T(0.1)],
        ρ=(T(0.1), T(0.1), zero(T)),
        # ♯=(10, 10))
        ♯=(Int(browser.width), Int(browser.height)))
    global godBROWSER = Ref(godbrowser(g, browser))
    @show "got godBROWSER"
end
const CHANGE_MODE = Ref(2) # 0=zero, 1=focus, 2=ρ
const CHANGE_DIM_INDEX = Ref(2)
function godbrowserkeypress(key)
    if key == "ArrowUp"
        if CHANGE_MODE[] == 0
            moveup!(godBROWSER[].g, CHANGE_DIM_INDEX[])
        elseif CHANGE_MODE[] == 1
            focusup!(godBROWSER[].g, CHANGE_DIM_INDEX[])
        else
            scaleup!(godBROWSER[].g, CHANGE_DIM_INDEX[])
        end
    elseif key == "ArrowDown"
        if CHANGE_MODE[] == 0
            movedown!(godBROWSER[].g, CHANGE_DIM_INDEX[])
        elseif CHANGE_MODE[] == 1
            focusdown!(godBROWSER[].g, CHANGE_DIM_INDEX[])
        else
            scaledown!(godBROWSER[].g, CHANGE_DIM_INDEX[])
        end
    elseif key == "0"
        global CHANGE_MODE[] = (CHANGE_MODE[] + 1) % 3
    elseif key == "["
        global CHANGEΔ *= ○
    elseif key == "]"
        global CHANGEΔ *= T(2)
    elseif key == "q"
        jerkup!(godBROWSER[].g)
    elseif key == "w"
        jerkdown!(godBROWSER[].g)
    else
        try
            global CHANGE_DIM_INDEX[] = parse(UInt, key)
        catch end
    end
    println("CHANGE_MODE=$CHANGE_MODE[]")
    println("CHANGE_DIM_INDEX=$CHANGE_DIM_INDEX[]")
    println("ẑero.μ=$(godBROWSER[].g.ẑero.μ)")
    println("ône.μ=$(godBROWSER[].g.ône.μ)")
    println("g.ρ=$(godBROWSER[].g.ρ)")
    println("g.θ=$(godBROWSER[].g.θ)")
end
# put!(::godBrowser) = nothing # todo ?
godBROWSER = nothing

function Δ!(ϕ, ϕ̇)
    δ = Tuple{CartesianIndex{2},Tuple{T,T,T,T}}[]
    for i = CartesianIndices(ϕ̇)
        ϕ[i] == ϕ̇[i] && continue
        ϕ[i] = ϕ̇[i]
        push!(δ, (i, scalar2rgba(ϕ̇[i])))
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
