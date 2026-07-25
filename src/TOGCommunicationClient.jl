module TOGCommunicationClient

export listentogroup, ignoregroup, CACHE, put!, readcache, Messages
export TOGMessage

using ZMQ
using TOGZMQ
using TOGZMQ: TOGMessage, state
using LoopOS: @whiletrue, Peripheral, listen
import Base: put!, take!

"""
Received messages are cached here. Retrieve them with `readcache`.
"""
const CACHE = TOGMessage[]
const DEALERSOCKET = Ref{Socket}()
const SUBSOCKET = Ref{Socket}()

struct Messages <: Peripheral
    channel::Channel{TOGMessage}
end
const MESSAGES = Messages(Channel{TOGMessage}())

put!(::Messages, to::String, togroup::Bool, description::String, information) = put!(Messages, to, togroup, description, information)
put!(::Type{Messages}, message::TOGMessage) = put!(Messages, message.to, message.togroup, message.description, message.information)
"""
Send any information to a specific god or to a group.
example: `put!(TOGCommunicationClient.Messages, "∀", true, "greeting", "hi")`
"""
put!(::Type{Messages}, to::String, togroup::Bool, description::String, information) = TOGZMQ.send(DEALERSOCKET[], TOGMessage(TOGZMQ.ID[], to, togroup, description, information))

"""
Retrieve received message using its index in CACHE.
example: `readcache(TOGCommunicationClient.Messages, index)`
"""
function readcache(::Type{Messages}, i)
    # @info "TOGCommunicationClient.jl, take!"
    message = CACHE[i]
    deleteat!(CACHE, i)
    message
end
take!(::Messages) = state(take!(MESSAGES.channel))
take!(::Type{Messages}) = take!(MESSAGES)

# __init__() = atexit(sleep)
# function sleep()
#     close(DEALERSOCKET[])
#     close(SUBSOCKET[])
# end
function awaken(;dealer, sub)
    # @info "TOGCommunicationClient.jl, awaken"
    DEALERSOCKET[] = Socket(DEALER)
    DEALERSOCKET[].routing_id = TOGZMQ.ID[]
    SUBSOCKET[] = Socket(SUB)
    connect(DEALERSOCKET[], dealer)
    connect(SUBSOCKET[], sub)
    listentogroup("∀")
    global MESSAGES
    listen(MESSAGES)
    dealertask = errormonitor(@async @whiletrue receive(DEALERSOCKET[]))
    subtask = errormonitor(@async @whiletrue receive(SUBSOCKET[]))
    dealertask, subtask
end

"""
Listen to messages from some group.
example: `listentogroup("interestinggroup")`
"""
listentogroup(group::String) = subscribe(SUBSOCKET[], group)
"""
Stop listening to messages from some group.
example: `ignoregroup("interestinggroup")`
"""
ignoregroup(group::String) = unsubscribe(SUBSOCKET[], group)

function receive(socket)
    message = TOGZMQ.receive(socket)
    push!(CACHE, message)
    put!(MESSAGES.channel, message)
end

end
