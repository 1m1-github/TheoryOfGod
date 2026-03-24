struct godBrowser
    g::god
    loop::Task
    browser::BroadcastBrowser
end
godBrowser2(g, browser) =
    begin
        @show "godBrowser2", typeof(browser)
        looptask = Threads.@spawn begin
            try
                put!(browser.processor, JS(g.♯[1], g.♯[2]))
                ϕ = zeros(T, g.♯[1], g.♯[2])
                while true
                    yield()
                    sleep(1) # todo rm
                    ϕ̇ = ∃̇(g, Ω[])
                    # @show p̂ixel
                    δ = Δ!(ϕ, ϕ̇)
                    # @show δ
                    # @show isempty(δ)
                    isempty(δ) && continue
                    js = "pixel=" * writeδ(δ, g.♯[1]) * "\n" * SET_PIXELS_JS
                    put!(browser.processor, js)
                end
            catch e
                bt = catch_backtrace()
                showerror(stderr, e, bt)
            end
        end
        godBrowser(g, looptask, browser)
    end
function godBrowser1(browser)
    @show "godBrowser1", typeof(browser)
    g = god(
        d=sort(SA[zero(T), invϕ, invϕ^2, one(T)]), # t, x, y, z
        ẑeroμ=SA[t(), ○, ○, ○],
        f̂ocusμ=SA[t(), ○*exp(T(0.1)), ○*exp(T(0.1)), ○*exp(T(0.1))],
        ρ=(T(0.1), T(0.1), one(T)),
        # ♯=(10, 10))
        ♯=(Int(browser.width), Int(browser.height)))
    @show "godBrowser1", typeof(g)
    gb = godBrowser2(g, browser)
    push!(godBROWSER[], gb)
    gb
end
# put!(::godBrowser) = nothing # todo ?
const godBROWSER = Ref(Set{godBrowser}())

# all(==(ntuple(_->one(T),4)),p̂ixel)
# g=gb.g
# gb=only(values(godBROWSER[]))
# browser=gb.browser

function Δ!(ϕ, ϕ̇)
    δ = Tuple{CartesianIndex{2},Tuple{T,T,T,T}}[]
    # i = collect(CartesianIndices(ϕ̇))[1]
    for i = CartesianIndices(ϕ̇)
        ϕ[i] == ϕ̇[i] && continue
        ϕ[i] = ϕ̇[i]
        push!(δ, (i, (c2rgb(ϕ̇[i])..., one(T))))
    end
    δ
end
# height=g.♯[2]
# color=δ[1][2]
# typemax(UInt8) .* color
function writeδ(δ, height)
    result = []
    for (i, color) = δ
        # @show i, color
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
