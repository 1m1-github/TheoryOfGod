module TOGgod

export learn, LoopOS

using Pkg, Serialization
using LoopOS, TOGObserveClient, TOGCreateClient, TOGLearning, TOGCommunicationClient, TOGAwaken, TOGLogging, TOGREPL, TOGBroadcastBrowser, TOGOctahedronBrowser, TOGAudioAnalogToDigitalBrowser, TOGTextToAudioBrowser, TOGVisualAnalogToDigitalBrowser, TOGPort
using TOGOctahedron: Octahedron
using TOG: ○

const T = Ref{DataType}()
const OCTAHEDRON = Ref{Octahedron}()
const ARGS = Ref{NamedTuple}((;))

__init__() = atexit(sleep)

function sleep()
    isempty(ARGS[]) && return
    serialize(".tog/short", LoopOS.short())
    TOGAwaken.sleep(path=ARGS[][:path])
    #     # TOGObserveClient.sleep()
    #     # TOGCreateClient.sleep()
    #     # TOGCommunicationClient.sleep()
    TOGBroadcastBrowser.sleep(path=ARGS[][:path])
    TOGREPL.sleep(path=ARGS[][:path])
end

function awaken(; args...)
    # write(string(time()),"god"*string(args))
    # @info "TOGgod,awaken", args
    ARGS[] = merge(ARGS[], args)
    ARGS[] = merge(ARGS[], [:path=>pwd()])
    remotereplport = get(args, :remotereplport, TOGPort.openport())
    ARGS[] = merge(ARGS[], [:remotereplport=>remotereplport])
    universe = args[:universe]
    name = args[:name]
    intelligence = args[:intelligence]
    TOGLogging.awaken()
    TOGAwaken.awaken()
    TOGObserveClient.awaken(TOGAwaken.togobserve(path=universe))
    TOGCreateClient.awaken(TOGAwaken.togcreate(path=universe))
    TOGCommunicationClient.awaken(name=name, dealer=TOGAwaken.router(path=universe), sub=TOGAwaken.pub(path=universe))
    T[] = TOGObserveClient.togtype()
    ϕ = MathConstants.golden
    OCTAHEDRON[] = Octahedron(
        t=TOGObserveClient.togtime(),
        d=[ϕ^-4, ϕ^-3, ϕ^-2, ϕ^-1],
        ẑeroμ=[zero(T[]), ○(T[]), ○(T[]), ○(T[])],
        ôneμ=[zero(T[]), ○(T[]), ○(T[]), ○(T[])+T[](0.1)],
        ρ=[T[](0.0), T[](0.1), T[](0.1), T[](0.0)],
        ♯=(1, 1))
    TOGBroadcastBrowser.awaken(root=browserconnect, functions=Dict(
        "/keypress"=>TOGOctahedronBrowser.keypress,
        "/websocket"=>TOGAudioAnalogToDigitalBrowser.ws,
        "/webcam"=>TOGVisualAnalogToDigitalBrowser.webcam,
    ))
    LoopOS.awaken(intelligence)
    # TOGREPL.awaken(name=name, remotereplport=remotereplport)
    # write(string(time()),"godend"*string(remotereplport))
    # atreplinit(r -> begin
    #     @show "A", TOGAwaken.remotereplportfile(path=universe), isfile(TOGAwaken.remotereplportfile(path=universe))
    #     @show "B", TOGAwaken.readremotereplport(path=universe)
    #     isfile(TOGAwaken.remotereplportfile(path=universe)) && TOGREPL.connect(start_key="\\C-g", port=TOGAwaken.readremotereplport(path=universe))
    # end)
end

function browserconnect(port, browser)
    # @show "browserconnect", port
    TOGOctahedronBrowser.awaken(octahedron=OCTAHEDRON[], browser=browser)
    TOGAudioAnalogToDigitalBrowser.awaken()
end

function learn(; pkgs=String[], rmpkgs=String[])
    # @show "learn", pkgs, rmpkgs
    isempty(pkgs) || Pkg.add(pkgs)
    isempty(rmpkgs) || Pkg.rm(rmpkgs)
    # @show "learn, Pkg.update"
    Pkg.update(update_registry=false) # DEBUG update_registry=false
    # @show "learn, serialize"
    # serialize(".tog/short", LoopOS.short())
    # @show "learn, TOGAwaken.awakengod"
    TOGAwaken.awakengod(; ARGS[]..., self=true)
    # @show "learn, exit"
    exit(0)
end

end
