module TOGObserveServer

using ZMQ
using TOGZMQAPIServer
using TOG∃: t, 𝕋, ∩

const SOCKET = Ref{Socket}()
const TASK = Ref{Task}()

# __init__() = atexit(sleep)
# function sleep()
    # schedule(TASK[], InterruptException(), error=true)
    # TOGZMQAPIServer.sleep(SOCKET[])
# end
function awaken(;socketlocation, ω)
    SOCKET[], TASK[] = TOGZMQAPIServer.awaken(socketlocation=socketlocation, functions=Dict(
        :time => time(ω),
        :T => type(ω),
        :type => type(ω),
        :∩ => ∩(ω),
    ))
end

time(ω::𝕋) = (x...) -> t(ω)
type(ω::𝕋) = (x...) -> first(typeof(ω).parameters)
Base.:∩(ω::𝕋) = (x...) -> 
begin
    # @info typeof(x)
    # ∩(x[1][1], ω, x[1][2]) # put!
    ∩(x[1], ω, x[2]) # take!(octahedron)
end

end
