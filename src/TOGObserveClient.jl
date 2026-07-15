module TOGObserveClient

using ZMQ, TOGZMQAPIClient
using TOG: ∃

const SOCKET = Ref{Socket}()

# __init__() = atexit(sleep)
awaken(socketlocation) = SOCKET[] = TOGZMQAPIClient.awaken(socketlocation)
togtime() = TOGZMQAPIClient.call(SOCKET[], :time)
togtype() = TOGZMQAPIClient.call(SOCKET[], :type)
Base.:∩(ϵ::∃, Ο) = TOGZMQAPIClient.call(SOCKET[], :∩, ϵ, Ο)
# sleep() = begin
#     write("qq1", string(typeof(SOCKET[])))
#     write("qq2", string(SOCKET))
#     TOGZMQAPIClient.sleep(SOCKET[])
# end

end
