module TOGZMQ

using ZMQ, Serialization
using LoopOS: Peripheral

import Base: take!, put!
# struct Message <: Peripheral end
struct ReceiveChannel <: Peripheral
    channel::Channel
end
take!(a::ReceiveChannel) = take!(a.channel)
put!(a::ReceiveChannel, information) = put!(a.channel, information)
const RECEIVECHANNEL = ReceiveChannel(Channel{Any}(Inf))

# __init__() = atexit(sleep)
# sleep() = nothing
# function sleep(socket::Socket)
#     # rm(replace(socket.last_endpoint[1:end-1], r"^ipc://" => "")) # todo needed?
#     close(socket)
# end
send(socket::Socket, information) = send(socket, "", false, "", Symbol(""), "", information)
function send(socket::Socket, to::String, togroup::Bool, from::String, symbol::Symbol, description::String, information)
    # @info "TOGZMQ.send", socket, to, togroup, from, symbol, description, typeof(information)
    buffer = IOBuffer()
    # @info "TOGZMQ.send", sizeof(buffer)
    serialize(buffer, information)
    # @info "TOGZMQ.send, serialized"
    _information = Message(take!(buffer))
    # @info "TOGZMQ.send, got _information"
    ZMQ.send_multipart(socket, [to, togroup, from, string(symbol), description, _information])
end
# TOGgod.TOGCommunicationClient.send("∀", true, "i", :q, "qw", 12)
# send(to::String, togroup::Bool, from::String, symbol::Symbol, description::String, information)
# "i","∀",true,"i","q","qw",12
function receive(socket::Socket)
    # @info "TOGZMQ.receive", socket
    frames = ZMQ.recv_multipart(socket)
    # @info "TOGZMQ.receive", socket, length(frames)
    # n = socket.type == ROUTER ? 
    # frames1 = String(frames[1])
    # @info "TOGZMQ.receive", frames1
    to = length(frames) < 6 ? "" : String(frames[end-5])
    # @info "TOGZMQ.receive, to", to
    togroup = Bool(only(frames[end-4]))
    # @info "TOGZMQ.receive, togroup", togroup
    from = String(frames[end-3])
    # @info "TOGZMQ.receive, from", from
    symbol = Symbol(String(frames[end-2]))
    # @info "TOGZMQ.receive, symbol", symbol
    description = String(frames[end-1])
    # @info "TOGZMQ.receive, description", description
    buffer = IOBuffer(frames[end])
    # @info "TOGZMQ.receive, buffer", sizeof(buffer)
    information = try 
        deserialize(buffer) 
    catch e 
        @info e
        buffer
    end
    # @info "TOGZMQ.receive", typeof(information)
    to, togroup, from, symbol, description, information
end

end
