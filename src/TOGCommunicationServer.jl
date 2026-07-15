module TOGCommunicationServer

using ZMQ
using LoopOS: @whiletrue
using TOGZMQ

# const GROUPS = Dict{String, Vector{String}}("∀" => [""])
# const GROUP = Dict{String,String}()
# const GROUPS = Set{String}()
const ROUTERSOCKET = Ref{Socket}()
const PUBSOCKET = Ref{Socket}()

# creategroup(;group) = push!(GROUPS, group)
# rmgroup(;group) = delete!(GROUPS, group)
# function addclient(;group, name)
#     to ∈ GROUPS && error("$name already exists as a client")
#     push!(GROUPS, to)
#     # TOGZMQClient.start(group=group, router=ROUTERSOCKET[].last_endpoint, pub=PUBSOCKET[].last_endpoint, name=name)
# end
# function creategroup(;groupname, creatorname)
#     haskey(GROUPS, groupname) || return
#     GROUPS[groupname] = [creatorname]
# end
# function addtogroup(;groupname, user, membername)
#     haskey(GROUPS, groupname) || return
#     group = GROUPS[groupname]
#     user == group[1] || return
#     membername ∈ values(group) && return
#     push!(group, membername)
# end
# function rmgroup(;groupname, user)
#     haskey(GROUPS, groupname) || return
#     user == GROUPS[groupname][1] || return
#     delete!(GROUPS, groupname)
# end
# function rmfromgroup(;groupname, user, membername)
#     haskey(GROUPS, groupname) || return
#     group = GROUPS[groupname]
#     user == group[1] || return
#     user == membername && return
#     group = filter!(n -> n ≠ membername, group)
# end
# addtogroup(groupname="∀", user=name, membername=name)

# __init__() = atexit(sleep)
# function sleep()
#     TOGZMQ.sleep(ROUTERSOCKET[])
#     TOGZMQ.sleep(PUBSOCKET[])
# end

function awaken(; router, pub)
    ROUTERSOCKET[] = Socket(ROUTER)
    PUBSOCKET[] = Socket(PUB)
    bind(ROUTERSOCKET[], router)
    bind(PUBSOCKET[], pub)
    # listen(RECEIVE)
    errormonitor(@async @whiletrue begin
        @info "TOGCommunicationServer, waiting"
        to, togroup, from, symbol, description, information = TOGZMQ.receive(ROUTERSOCKET[])
        @info "TOGCommunicationServer, received", to, togroup, from, symbol, description, typeof(information), information
        # haskey(AGENT, from) || continue
        # TOGZMQ.send(PUBSOCKET[], to, from, information)
        # if to == "group"
        # if to ∈ GROUPS
        socket = togroup ? PUBSOCKET[] : ROUTERSOCKET[]
        # if togroup
        # group = AGENTGROUP[from]
        # TOGZMQ.send(PUBSOCKET[], [to, from, symbol, description, information])
        # send_multipart(PUBSOCKET[], [group, from, message])
        # elseif AGENTGROUP[from] == get(AGENTGROUP, to, "")
        # else
        # send_multipart(ROUTERSOCKET[], [to, from, message])
        # end
        TOGZMQ.send(socket, to, togroup, from, symbol, description, information)
    end)
end

# import Base: take!, put!
# send(socket, message, to) = ZMQ.send_multipart(socket, [to, Sys.username(), message])
# put!(::DirectMessage, message::String, to::String) = send(ROUTERSOCKET[], message, to)
# put!(::GroupMessage, message::String, to::String="∀") = send(PUBSOCKET[], message, to)

end
