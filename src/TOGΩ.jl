module TOGΩ

# export TOGREPL

const T = Float64

using Serialization
using TOG∃: 𝕋
using TOGCommunicationServer, TOGAwaken, TOGLogging, TOGObserveServer, TOGCreateServer, TOGZMQ
using TOGREPL
using TOGMatrix, TOGOctahedron, TOGColor, Colors # DEBUG?

const Ωpath = joinpath(TOGAwaken.TOGDIR, "Ω")
const Ω = isfile(Ωpath) ? deserialize(Ωpath) : 𝕋(T)

__init__() = atexit(sleep)
function sleep(exitcode)
    @info "sleep", exitcode
    exitcode == TOGAwaken.ALREADYRUNNINGEXITCODE && return
    serialize(Ωpath, Ω)
    TOGAwaken.sleep()
#     TOGObserveServer.sleep()
#     TOGCreateServer.sleep()
#     TOGCommunicationServer.sleep()
    TOGREPL.sleep()
end
function awaken(; path=".", router=TOGAwaken.router(path=path), pub=TOGAwaken.pub(path=path), togobserve=TOGAwaken.togobserve(path=path), togcreate=TOGAwaken.togcreate(path=path))
    @info "TOGΩ, awaken"
    TOGLogging.awaken()
    TOGAwaken.awaken()
    TOGZMQ.awaken(name="Ω")
    TOGObserveServer.awaken(socketlocation=togobserve, ω=Ω)
    TOGCreateServer.awaken(socketlocation=togcreate, ω=Ω)
    TOGCommunicationServer.awaken(router=router, pub=pub)
    TOGREPL.awaken(name="Ω")
end

end
