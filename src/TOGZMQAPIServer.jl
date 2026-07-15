module TOGZMQAPIServer

using ZMQ, Serialization
using LoopOS: @whiletrue
using TOGZMQ

struct APIData
    f::Symbol
    args::Any
end
const FUNCTIONS = Dict{Symbol,Function}()

# __init__() = atexit(sleep)
# sleep() = nothing
# sleep(socket) = TOGZMQ.sleep(socket)
function awaken(;socketlocation, functions)
    # @info "TOGZMQAPIServer.awaken", socketlocation
    socket = Socket(REP)
    bind(socket, socketlocation)
    socket, errormonitor(@async @whiletrue receive(socket, functions))
end
function receive(socket, functions)
    # @info "TOGZMQAPIServer.receive", socket, functions
    _, _, _, _, _, apidata = TOGZMQ.receive(socket)
    # @info "TOGZMQAPIServer.receive", apidata.f
    output = try
        f = functions[apidata.f]
        # @info typeof(apidata.args), typeof(f)
        isnothing(apidata.args) ? f() : f(apidata.args...)
        # functions[apidata.f](apidata.args...)
    catch e
        @info e
        bt = catch_backtrace()
        showerror(stdout, e, bt)
        string(e)
    end
    # @info "TOGZMQAPIServer.receive", typeof(output)
    TOGZMQ.send(socket, output)
end

end
