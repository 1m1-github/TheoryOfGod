module TOGZMQAPIServer

using ZMQ, Serialization
using LoopOS: @whiletrue
using TOGZMQ
using TOGZMQ: TOGMessage

struct APIData
    f::Symbol
    args::Any
end
const FUNCTIONS = Dict{Symbol,Function}()

# __init__() = atexit(sleep)
# sleep() = nothing
# sleep(socket) = TOGZMQ.sleep(socket)
function awaken(;socketlocation, functions)
    @info "TOGZMQAPIServer.awaken", socketlocation
    socket = Socket(REP)
    bind(socket, socketlocation)
    socket, errormonitor(@async @whiletrue receive(socket, functions))
end
function receive(socket, functions)
    # @info "TOGZMQAPIServer.receive", socket, functions
    message = TOGZMQ.receive(socket)
    @info "TOGZMQAPIServer.receive", message
    output = try
        f = functions[message.information.f]
        # @info typeof(message.information.args), typeof(f)
        isnothing(message.information.args) ? f() : f(message.information.args...)
        # functions[message.information.f](message.information.args...)
    catch e
        @info e
        bt = catch_backtrace()
        showerror(stdout, e, bt)
        string(e)
    end
    @info "TOGZMQAPIServer.receive", typeof(output)
    TOGZMQ.send(socket, TOGMessage(TOGZMQ.ID[], message.from, message.togroup, message.description, output))
end

end
