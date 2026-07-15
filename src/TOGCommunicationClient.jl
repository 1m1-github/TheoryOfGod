module TOGCommunicationClient

using ZMQ
using TOGZMQ
using LoopOS: @whiletrue, Peripheral
# import TOGZMQ: Message
import Base: put!, take!

const CACHE = Dict{Symbol, Tuple{String, Any}}()
const DEALERSOCKET = Ref{Socket}()
const SUBSOCKET = Ref{Socket}()

struct Message <: Peripheral end

put!(::Type{Message}, to::String, togroup::Bool, symbol::Symbol, description::String, information) = TOGZMQ.send(DEALERSOCKET[], to, togroup, DEALERSOCKET[].routing_id, symbol, description, information)
put!(::Type{Message}, to::String, togroup::Bool, information) = put!(Message, to, togroup, information)
put!(::Type{Message}, to::String, information) = put!(Message, to, false, information)

"Returns cached information given a symbol and removes it from CACHE"
function take!(::Type{Message}, symbol::Symbol)
    haskey(CACHE, symbol) || error("No information found for symbol :$symbol.")
    description, information = CACHE[symbol]
    delete!(CACHE, symbol)
    description, information
end

# __init__() = atexit(sleep)
# function sleep()
#     close(DEALERSOCKET[])
#     close(SUBSOCKET[])
# end
function awaken(;name, dealer, sub)
    # @show "TOGCommunicationClient.awaken", name, dealer, sub
    DEALERSOCKET[] = Socket(DEALER)
    # @show "TOGCommunicationClient.awaken", DEALERSOCKET[]
    # setproperty!(DEALERSOCKET[], :routing_id, name)
    DEALERSOCKET[].routing_id = name
    SUBSOCKET[] = Socket(SUB)
    # @show "TOGCommunicationClient.awaken", SUBSOCKET[]
    connect(DEALERSOCKET[], dealer)
    connect(SUBSOCKET[], sub)
    # @show "TOGCommunicationClient.awaken sockets connected"
    listen(group="∀")
    # @show "TOGCommunicationClient.awaken listening to group ∀"
    # LoopOS.listen(TOGZMQ.RECEIVECHANNEL)
    # @show "TOGCommunicationClient.awaken listening to TOGZMQ.RECEIVECHANNEL"
    dealertask = errormonitor(@async @whiletrue receive(DEALERSOCKET[]))
    subtask = errormonitor(@async @whiletrue receive(SUBSOCKET[]))
    dealertask, subtask
end

listen(;group) = subscribe(SUBSOCKET[], group)
ignore(;group) = unsubscribe(SUBSOCKET[], group)

"Stores received in CACHE"
function receive(socket)
    @info "TOGCommunicationClient.receive", socket
    to, togroup, from, symbol, description, information = TOGZMQ.receive(socket)
    @info "TOGCommunicationClient, received", to, togroup, from, symbol, description, typeof(information), information
    CACHE[symbol] = (description, information)
    imessage = string(information isa AbstractString ? information : description * " of type " * string(typeof(information)) * " (use x=TOGCommunicationClient.take!(:$symbol) to store the information in x)")
    socket.type == DEALER && ( to = ZMQ.getproperty(socket, :routing_id) )
    fromto = "$from>$(togroup ? "GROUP-" : "")$to>"
    message = fromto * imessage * "."
    @info "TOGCommunicationClient, message", message
    put!(TOGZMQ.RECEIVECHANNEL, message)
end

end
