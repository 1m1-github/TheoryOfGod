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
    # TOGREPL.sleep()
end
function awaken(; path=".", router=TOGAwaken.router(path=path), pub=TOGAwaken.pub(path=path), togobserve=TOGAwaken.togobserve(path=path), togcreate=TOGAwaken.togcreate(path=path))
    @info "TOGΩ, awaken"
    TOGLogging.awaken()
    TOGAwaken.awaken()
    @info "TOGΩ, awaken2"
    TOGZMQ.awaken(name="Ω")
    @info "TOGΩ, awaken3"
    TOGObserveServer.awaken(socketlocation=togobserve, ω=Ω)
    @info "TOGΩ, awaken4"
    TOGCreateServer.awaken(socketlocation=togcreate, ω=Ω)
    @info "TOGΩ, awaken5"
    TOGCommunicationServer.awaken(router=router, pub=pub)
    @info "TOGΩ, awaken6"
    # TOGREPL.awaken(name="Ω")
    wait(Condition())
    @info "TOGΩ, awaken7"
end

end
