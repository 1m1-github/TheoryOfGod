module TOGZMQ

export TOGMessage

using ZMQ, Serialization
using LoopOS: Peripheral
import TOGState: state

import Base: take!, put!

const ID = Ref{String}("")

struct TOGMessage
    from::String
    to::String
    togroup::Bool
    description::String
    information::Any
end
state(message::TOGMessage) = """$(message.from)->$(message.to): $(typeof(message.information)) [$(message.description)]"""

awaken(;name) = ID[] = name

# __init__() = atexit(sleep)
# sleep() = nothing
# function sleep(socket::Socket)
#     # rm(replace(socket.last_endpoint[1:end-1], r"^ipc://" => "")) # todo needed?
#     close(socket)
# end
function send(socket::Socket, message::TOGMessage)
    # isempty(ID[]) && error("Run TOGZMQ.awaken() first.")
    # @info "TOGZMQ.send", socket, message
    buffer = IOBuffer()
    # @info "TOGZMQ.send", sizeof(buffer)
    serialize(buffer, message.information)
    # @info "TOGZMQ.send, serialized", typeof(message.information), message.information
    information = Message(take!(buffer))
    # @info "TOGZMQ.send, buffer", typeof(information), information
    data = [message.to, message.togroup, message.from, message.description]
    # @info "TOGZMQ.send data1", data
    # @info "TOGZMQ.send socket.type", socket.type, socket.type == ROUTER
    socket.type == ROUTER && ( data = [message.to; data...] )
    # @info "TOGZMQ.send data2", data
    push!(data, information)
    # @info "TOGZMQ.send data3", data
    ZMQ.send_multipart(socket, data)
end
function receive(socket::Socket)
    # @info "TOGZMQ.receive", socket
    frames = ZMQ.recv_multipart(socket)
    # @info "TOGZMQ.receive", socket, length(frames)
    to = String(frames[end-4])
    # @info "TOGZMQ.receive, to", to
    togroup = Bool(only(frames[end-3]))
    # @info "TOGZMQ.receive, togroup", togroup
    from = String(frames[end-2])
    # @info "TOGZMQ.receive, from", from
    description = String(frames[end-1])
    # @info "TOGZMQ.receive, from", from
    buffer = IOBuffer(frames[end])
    # @info "TOGZMQ.receive, buffer", sizeof(buffer)
    information = try 
        deserialize(buffer) 
    catch e
        e
    end
    # @info "TOGZMQ.receive", socket, information
    TOGMessage(from, to, togroup, description, information)
end

end
