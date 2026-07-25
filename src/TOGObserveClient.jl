module TOGObserveClient

export ∩

using ZMQ, TOGZMQAPIClient
using TOGExist: ∃

const SOCKET = Ref{Socket}()

# __init__() = atexit(sleep)
# sleep() = TOGZMQAPIClient.sleep(SOCKET[])
awaken(socketlocation) = SOCKET[] = TOGZMQAPIClient.awaken(socketlocation)
togtime() = TOGZMQAPIClient.call(SOCKET[], :time)
# togtype() = TOGZMQAPIClient.call(SOCKET[], :type)
"""
Find intersecting existences.
# Arguments
- `ϵ::∃`: Box to find intersecting existences in.
- `t::Number`: Max creation time of intersecting existences.
"""
Base.:∩(ϵ::∃, t::Number) = TOGZMQAPIClient.call(SOCKET[], :∩, ϵ, t)

end
