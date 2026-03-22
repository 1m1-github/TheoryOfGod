const MAX_RGB = T(mfb_rgb(255, 255, 255))
rgb2c(r, g, b) = T(mfb_rgb(r * 255, g * 255, b * 255)) / MAX_RGB
c2rgb(c2) = begin
    c = floor(UInt32, c2 * MAX_RGB)
    ((c >> 16) & 0xFF, (c >> 8) & 0xFF, c & 0xFF) ./ 255
end

window = mfb_open_ex("tog", g.♯..., MiniFB.WF_RESIZABLE)
buffer = zeros(UInt32, prod(g.♯))
MINIFBTASK = @async while true
    yield()
    state = mfb_update(window, buffer)
    state != MiniFB.STATE_OK && break
end
viewer_zero_vs_one_mode = 0
viewer_dim = g.ẑero.d[2]
function keyboard_cb(window::Ptr{Cvoid}, key::Int32, mod::Int32, is_pressed::Bool)::Cvoid
    try
        global g
        global viewer_dim
        global viewer_zero_vs_one_mode
        # global buffer
        println("Key pressed: $key, $mod, $(is_pressed)")
        println("viewer_zero_vs_one_mode=$viewer_zero_vs_one_mode")
        println("viewer_dim=$viewer_dim")
        println("ẑero.μ=$(g.ẑero.μ)")
        println("ône.μ=$(g.ône.μ)")
        if is_pressed
            update = false
            if key == 48 # 0
                viewer_zero_vs_one_mode = 1 - viewer_zero_vs_one_mode
            elseif 49 ≤ key ≤ 49 + 8 # numbers
                viewer_dim = g.ẑero.d[key+1-49]
            elseif key == 265 # up
                if iszero(viewer_zero_vs_one_mode)
                    g = moveup(g, viewer_dim)
                else
                    g = focusup(g, viewer_dim)
                end
                update = true
            elseif key == 264 # down
                if iszero(viewer_zero_vs_one_mode)
                    g = movedown(g, viewer_dim)
                else
                    g = focusdown(g, viewer_dim)
                end
                update = true
            elseif key == 81 # q
                g = jerkdown(g)
                update = true
            elseif key == 87 # w
                g = jerkup(g)
                update = true
            elseif key == 69 # e
                g = scaledown(g)
                update = true
            elseif key == 82 # r
                g = scaleup(g)
                update = true
            end
            update && ( global buffer = floor.(UInt32, reshape(∃̇(g, 10), prod(g.♯)) .* MAX_RGB) )
        end
    catch e
        bt = catch_backtrace()
        showerror(stderr, e, bt)
        println("after:")
        println("viewer_zero_vs_one_mode=$viewer_zero_vs_one_mode")
        println("viewer_dim=$viewer_dim")
        println("ẑero.μ=$(g.ẑero.μ)")
        println("ône.μ=$(g.ône.μ)")
    end
    return nothing
end
kb_cfunc = @cfunction(keyboard_cb, Cvoid, (Ptr{Cvoid}, Int32, Int32, Bool))
ccall((:mfb_set_keyboard_callback, MiniFB.libminifb), Cvoid,
    (Ptr{Cvoid}, Ptr{Cvoid}),
    window, kb_cfunc)
# schedule(MINIFBTASK, InterruptException(), error=true)
# mfb_close(window)
