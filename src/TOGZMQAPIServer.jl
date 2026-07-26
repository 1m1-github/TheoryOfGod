module TOGZMQAPIServer

using ZMQ, Serialization
using LoopOS: @whiletrue
using TOGZMQ
using TOGZMQ: TOGMessage
using TOGLogging: exception_details

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
    message = TOGZMQ.receive(socket)
    output = try
        f = functions[message.information.f]
        isnothing(message.information.args) ? f() : f(message.information.args...)
    catch e
        @error exception_details(e, catch_backtrace())
        string(e)
    end
    TOGZMQ.send(socket, TOGMessage(TOGZMQ.ID[], message.from, message.togroup, message.description, output))
end

end
