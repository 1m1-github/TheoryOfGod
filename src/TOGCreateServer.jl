module TOGCreateServer

using ZMQ
using TOGZMQAPIServer
using TOGExist: ∃, ∃!, 𝕋

const SOCKET = Ref{Socket}()
const TASK = Ref{Task}()

# __init__() = atexit(sleep)
# function sleep()
#     schedule(TASK[], InterruptException(), error=true)
#     TOGZMQAPIServer.sleep(SOCKET[])
# end
awaken(;socketlocation, ω::𝕋) = SOCKET[], TASK[] = TOGZMQAPIServer.awaken(socketlocation=socketlocation, functions=Dict(:create => create(ω)))
create(ω::𝕋) = (x...) -> create(x..., ω)
create(ϵ::∃, name::String, ω::𝕋) = ∃!(ϵ, name, ω)

end
