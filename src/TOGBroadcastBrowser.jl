module TOGBroadcastBrowser

export put!, anybrowserconnected

using HTTP, URIs, Sockets
using LoopOS: BatchProcessor, start!, Peripheral
using TOGAwaken, TOGPort
import Base.put!

"Serve and execute JavaScript on an HTTP client using SSE"
struct BroadcastBrowser <: Peripheral
    stream::HTTP.Stream
    width::Int
    height::Int
    processor::BatchProcessor{String}
    BroadcastBrowser(stream, width, height) = new(stream, width, height, BatchProcessor{String}())
end
# const BROADCASTBROWSER = Ref{BroadcastBrowser}()
const CLIENTS = Ref(Set{BroadcastBrowser}())
"""
Run js on all connected browsers.
example: `put!(BroadcastBrowser, "console.log('hi')")`.
"""
put!(::Type{BroadcastBrowser}, js::String) = [put!(client.processor, js) for client = CLIENTS[]]
# put!(::Type{BroadcastBrowser}, js) = begin
# @info "TOGBroadcastBrowser.jl, put!"    
#     put!(BROADCASTBROWSER[].processor, js)
# end

# window.WS = new WebSocket('wss://studio.tail16337b.ts.net')
const HTMLINIT(port) = """
<!DOCTYPE html>
<html>
<body>
<script>
window.SSE = new EventSource(`/events?width=\${document.documentElement.clientWidth}&height=\${document.documentElement.clientHeight}`)
window.SSE.onmessage = (e) => {console.log(e.data);eval(e.data);}
</script>
</body>
</html>
"""

function safe_write(stream, js)
    try
        write(stream, js)
        flush(stream)
        true
    catch e
        e isa Base.IOError || rethrow()
        false
    end
end

function handle_sse(bb)
    HTTP.setstatus(bb.stream, 200)
    HTTP.setheader(bb.stream, "Content-Type" => "text/event-stream")
    HTTP.setheader(bb.stream, "Cache-Control" => "no-cache")
    HTTP.startwrite(bb.stream)
    start!(bb.processor) do input
        for js = input
            for line = split(js, '\n')
                safe_write(bb.stream, "data: $line\n") || return
            end
            safe_write(bb.stream, "\n") || return
        end
    end
end

function handle_ws(stream, f)
    HTTP.WebSockets.upgrade(stream) do ws
        for msg in ws HTTP.WebSockets.send(ws, f(msg)) end
    end
end

anybrowserconnected = false
function awaken(; root::Function, port=TOGPort.openport(), functions=Dict("/websocket"=>identity))
    # @info "TOGBroadcastBrowser.jl, awaken"    
    TOGAwaken.writebroadcastbrowserport(port=port)
        @async HTTP.listen!("127.0.0.1", port) do stream
        # if HTTP.WebSockets.isupgrade(stream.message)
        #     @async handle_ws(stream, functions["/websocket"])
        #     return
        # end
        target = stream.message.target
        uri = URI(target)
        # @info "TOGBroadcastBrowser, target, uri", target, uri, uri.path, haskey(functions, uri.path)
        if target == "/"
            HTTP.setstatus(stream, 200)
            HTTP.setheader(stream, "Content-Type" => "text/html")
            HTTP.startwrite(stream)
            write(stream, HTMLINIT(port))
        elseif uri.path == "/events"
            params = queryparams(uri)
            width = parse(Int, params["width"])
            height = parse(Int, params["height"])
            # BROADCASTBROWSER[] = BroadcastBrowser(stream, width, height)
            bb = BroadcastBrowser(stream, width, height)
            push!(CLIENTS[], bb)
            anybrowserconnected = true
            # root(port, BROADCASTBROWSER[])
            # handle_sse(BROADCASTBROWSER[])
            root(port, bb)
            handle_sse(bb)
            delete!(CLIENTS[], bb)
            anybrowserconnected = !isempty(CLIENTS[])
        elseif haskey(functions, uri.path)
            functions[uri.path](read(stream))
            HTTP.setstatus(stream, 204)
            HTTP.startwrite(stream)
        else
            HTTP.setstatus(stream, 404)
            HTTP.startwrite(stream)
        end
    end
    # HTTP.get("http://localhost:"*string(port)) # todo needed to reserve port?
end

function sleep(; path=".")
    filename = TOGAwaken.broadcastbrowserportfile(path=path)
    isfile(filename) && rm(filename)
end

end
