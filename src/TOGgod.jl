module TOGgod

export learn, LoopOS, OCTAHEDRON, T, put!

using Pkg, Serialization
using LoopOS, TOGObserveClient, TOGCreateClient, TOGLearning, TOGCommunicationClient, TOGAwaken, TOGLogging, TOGREPL, TOGBroadcastBrowser, TOGOctahedronBrowser, TOGAudioAnalogToDigitalBrowser, TOGTextToAudioBrowser, TOGVisualAnalogToDigitalBrowser, TOGPort, TOGZMQ
using LoopOS: Peripheral
using TOGOctahedron: Octahedron
using TOGExist: ○
import Base: put!

"eltype of Ω."
const T = Ref{DataType}()
"Main octahedron."
const OCTAHEDRON = Ref{Octahedron}()
const ARGS = Ref{NamedTuple}((;))

__init__() = atexit((n)->begin
    @info "TOGgod, atexit"
    sleep(n)
end
)

struct Ω <: Peripheral end
put!(::Type{Ω}, args...) = TOGCreateClient.put!(OCTAHEDRON[], args...)

function sleep(exitcode)
    @info "TOGgod, sleep", exitcode
    exitcode == TOGAwaken.ALREADYRUNNINGEXITCODE && return
    isempty(ARGS[]) && return
    serialize(".tog/short", LoopOS.short())
    TOGAwaken.sleep(path=ARGS[][:path])
    #     # TOGObserveClient.sleep()
    #     # TOGCreateClient.sleep()
    #     # TOGCommunicationClient.sleep()
    TOGBroadcastBrowser.sleep(path=ARGS[][:path])
    # TOGREPL.sleep(path=ARGS[][:path])
end

function awaken(; args...)
    TOGLogging.awaken()
    # write(string(time()),"god"*string(args))
    @info "TOGgod.awaken", args
    ARGS[] = merge(ARGS[], args)
    ARGS[] = merge(ARGS[], [:path=>pwd()])
    remotereplport = get(args, :remotereplport, TOGPort.openport())
    broadcastbrowserport = get(args, :broadcastbrowserport, TOGPort.openport())
    ARGS[] = merge(ARGS[], [:remotereplport=>remotereplport, :broadcastbrowserport=>broadcastbrowserport])
    universe = args[:universe]
    name = basename(args[:path])
    intelligence = args[:intelligence]
    TOGAwaken.awaken()
    TOGZMQ.awaken(name=name)
    TOGObserveClient.awaken(TOGAwaken.togobserve(path=universe))
    TOGCreateClient.awaken(TOGAwaken.togcreate(path=universe))
    TOGCommunicationClient.awaken(dealer=TOGAwaken.router(path=universe), sub=TOGAwaken.pub(path=universe))
    T[] = TOGObserveClient.togtype()
    ϕ = MathConstants.golden
    OCTAHEDRON[] = Octahedron(
        t=TOGObserveClient.togtime(),
        d=[ϕ^-4, ϕ^-3, ϕ^-2, ϕ^-1],
        ẑeroμ=[zero(T[]), ○(T[]), ○(T[]), ○(T[])],
        ôneμ=[zero(T[]), ○(T[]), ○(T[]), ○(T[])+T[](0.1)],
        ρ=[T[](0.0), T[](0.1), T[](0.1), T[](0.0)],
        ♯=(1, 1))
    @info "TOGgod, broadcastbrowserport", broadcastbrowserport
    TOGBroadcastBrowser.awaken(root=browserconnect, port=broadcastbrowserport, functions=Dict(
        "/keypress"=>TOGOctahedronBrowser.keypress,
        # "/websocket"=>TOGAudioAnalogToDigitalBrowser.ws,
        "/audio"=>TOGAudioAnalogToDigitalBrowser.ws,
        # "/webcam"=>TOGVisualAnalogToDigitalBrowser.webcam,
    ))
    @info "TOGgod, after TOGBroadcastBrowser.awaken"
    # LoopOS.awaken(intelligence)
    # @info "TOGgod, after LoopOS.awaken"
    # TOGREPL.awaken(name=name, remotereplport=remotereplport)
    # @info "TOGgod, after TOGREPL.awaken"
    # write(string(time()),"godend"*string(remotereplport))
    # atreplinit(r -> begin
    #     @show "A", TOGAwaken.remotereplportfile(path=universe), isfile(TOGAwaken.remotereplportfile(path=universe))
    #     @show "B", TOGAwaken.readremotereplport(path=universe)
    #     isfile(TOGAwaken.remotereplportfile(path=universe)) && TOGREPL.connect(start_key="\\C-g", port=TOGAwaken.readremotereplport(path=universe))
    # end)
    isinteractive() ? nothing : wait(Condition())
    @info "TOGgod, after wait"
end

function browserconnect(port, browser)
    @info "TOGgod.browserconnect", port
    TOGOctahedronBrowser.awaken(octahedron=OCTAHEDRON[], browser=browser)
    # TOGAudioAnalogToDigitalBrowser.awaken()
end

"""
Learn by adding or removing Pkgs from yourself.
"""
function learn(; pkgs=String[], rmpkgs=String[])
    @info "TOGgod, learn"
    isempty(pkgs) || Pkg.add(pkgs)
    isempty(rmpkgs) || Pkg.rm(rmpkgs)
    sleep(0)
    @info "TOGgod, learn, after sleep"
    TOGAwaken.awakengod(; ARGS[]...)
    @info "TOGgod, learn, after awakengod"
    exit(0)
end

end
