module TOGZMQAPIClient

using ZMQ
using TOGZMQAPIServer, TOGZMQ
using TOGZMQ: TOGMessage

# __init__() = atexit(sleep)
# sleep() = nothing
# sleep(socket::Socket) = TOGZMQAPIServer.sleep(socket)
function awaken(socketlocation)
    @show "TOGZMQAPIClient.awaken", socketlocation
    socket = Socket(REQ)
    connect(socket, socketlocation)
    socket
end
function call(socket::Socket, f::Symbol, x...)
    TOGZMQ.send(socket, TOGMessage(TOGZMQ.ID[], "Ω", false, string(f), TOGZMQAPIServer.APIData(f, x)))
    message = TOGZMQ.receive(socket)
    message.information
end
call(socket::Socket, f::Symbol) = call(socket, f, nothing)

end
