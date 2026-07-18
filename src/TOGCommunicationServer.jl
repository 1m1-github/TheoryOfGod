module TOGCommunicationServer

using ZMQ
using LoopOS: @whiletrue
using TOGZMQ
using TOGZMQ: TOGMessage

const ROUTERSOCKET = Ref{Socket}()
const PUBSOCKET = Ref{Socket}()

# __init__() = atexit(sleep)
# function sleep()
#     TOGZMQ.sleep(ROUTERSOCKET[])
#     TOGZMQ.sleep(PUBSOCKET[])
# end

function awaken(; router, pub)
    @info "TOGCommunicationServer.jl, .awaken", router, pub
    ROUTERSOCKET[] = Socket(ROUTER)
    PUBSOCKET[] = Socket(PUB)
    bind(ROUTERSOCKET[], router)
    bind(PUBSOCKET[], pub)
    errormonitor(@async @whiletrue begin
        # @info "TOGCommunicationServer, waiting"
        message = TOGZMQ.receive(ROUTERSOCKET[])
        @info "TOGCommunicationServer.jl, received",  message.from, message.to, message.togroup, message.description, typeof(message.information)
        socket =  message.togroup ? PUBSOCKET[] : ROUTERSOCKET[]
        TOGZMQ.send(socket, message)
    end)
end

end
