module TOGCommunicationClient

export listentogroup, ignoregroup, CACHE, put!, take!, Messages

using ZMQ
using TOGZMQ
using TOGZMQ: TOGMessage, state
using LoopOS: @whiletrue, Peripheral, listen
import Base: put!, take!

const CACHE = TOGMessage[]
const DEALERSOCKET = Ref{Socket}()
const SUBSOCKET = Ref{Socket}()

struct Messages <: Peripheral
    channel::Channel{TOGMessage}
end
const MESSAGES = Messages(Channel{TOGMessage}())

# put!(::Messages, message::TOGMessage)
put!(::Messages, to::String, togroup::Bool, description::String, information) = put!(Messages, to, togroup, description, information)
# put!(::Type{Messages}, message::TOGMessage) = put!(MESSAGES.channel, message)
put!(::Type{Messages}, message::TOGMessage) = put!(Messages, message.to, message.togroup, message.description, message.information)
"""
Send any information to a specific god or to a group.
example: `put!(Messages, "∀", true, "greeting", "hi")`
"""
put!(::Type{Messages}, to::String, togroup::Bool, description::String, information) = TOGZMQ.send(DEALERSOCKET[], TOGMessage(TOGZMQ.ID[], to, togroup, description, information))

"""
Retrieve received message using its index in CACHE.
example: `take!(Messages, 1)`
"""
function take!(::Type{Messages}, i)
    @info "TOGCommunicationClient.jl, take!"
    message = CACHE[i]
    deleteat!(CACHE, i)
    message
end
# take!(::Messages, i) = take!(Messages, i)
take!(::Messages) = state(take!(MESSAGES.channel))
take!(::Type{Messages}) = take!(MESSAGES)

# __init__() = atexit(sleep)
# function sleep()
#     close(DEALERSOCKET[])
#     close(SUBSOCKET[])
# end
function awaken(;dealer, sub)
    @info "TOGCommunicationClient.jl, awaken"
    DEALERSOCKET[] = Socket(DEALER)
    DEALERSOCKET[].routing_id = TOGZMQ.ID[]
    # @info "TOGCommunicationClient.awaken", DEALERSOCKET[].routing_id
    SUBSOCKET[] = Socket(SUB)
    connect(DEALERSOCKET[], dealer)
    connect(SUBSOCKET[], sub)
    listentogroup("∀")
    # global MESSAGES
    listen(MESSAGES)
    dealertask = errormonitor(@async @whiletrue receive(DEALERSOCKET[]))
    subtask = errormonitor(@async @whiletrue receive(SUBSOCKET[]))
    dealertask, subtask
end

"""
Join a group.
example: `listentogroup("interestinggroup")`
"""
listentogroup(group) = subscribe(SUBSOCKET[], group)
"""
Leave a group.
example: `ignoregroup("interestinggroup")`
"""
ignoregroup(group) = unsubscribe(SUBSOCKET[], group)

function receive(socket)
    @info "TOGCommunicationClient.jl, receive", socket
    message = TOGZMQ.receive(socket)
    @info "TOGCommunicationClient.jl, received",  message.from, message.to, message.togroup, message.description, typeof(message.information)
    push!(CACHE, message)
    put!(MESSAGES.channel, message)
end

end
