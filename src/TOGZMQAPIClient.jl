module TOGZMQAPIClient

using ZMQ
using TOGZMQAPIServer, TOGZMQ

# __init__() = atexit(sleep)
# sleep() = nothing
# sleep(socket::Socket) = TOGZMQAPIServer.sleep(socket)
function awaken(socketlocation)
    # @show "TOGZMQAPIClient.awaken", socketlocation
    socket = Socket(REQ)
    connect(socket, socketlocation)
    socket
end
function call(socket::Socket, f::Symbol, x...)
    # @info "TOGZMQAPIClient.call", socket, f, typeof(x), length(x)
    TOGZMQ.send(socket, TOGZMQAPIServer.APIData(f, x))
    # @info "TOGZMQAPIClient.call, sent"
    _, _, _, _, _, information = TOGZMQ.receive(socket)
    # @info "TOGZMQAPIClient.call", typeof(information)
    information
end
call(socket::Socket, f::Symbol) = call(socket, f, nothing)

end
