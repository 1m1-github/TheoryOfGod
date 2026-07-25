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
    # @info "TOGZMQAPIServer.awaken", socketlocation
    socket = Socket(REP)
    bind(socket, socketlocation)
    socket, errormonitor(@async @whiletrue receive(socket, functions))
end
function exception_details(e, bt = nothing)
    estr = if e isa TaskFailedException
        "TaskFailedException:\n" * exception_details(e.task.exception)
    elseif e isa CompositeException
        join(exception_details.(e.exceptions), "\n")
    else
        string(e)
    end 
    estr * "\n" * sprint(Base.show_backtrace, catch_backtrace())
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
