module TOGgod

# export learn, LoopOS
export LoopOS

using Pkg, Serialization
using LoopOS, TOGObserveClient, TOGCreateClient, TOGLearning, TOGCommunicationClient, TOGAwaken, TOGLogging, TOGREPL, TOGBroadcastBrowser, TOGOctahedronBrowser, TOGAudioAnalogToDigitalBrowser, TOGTextToAudioBrowser, TOGVisualAnalogToDigitalBrowser, TOGPort, TOGZMQ
using TOGOctahedron: Octahedron
using TOG∃: ○

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
    @info "TOGgod.awaken", args
    ARGS[] = merge(ARGS[], args)
    ARGS[] = merge(ARGS[], [:path=>pwd()])
    remotereplport = get(args, :remotereplport, TOGPort.openport())
    broadcastbrowserport = get(args, :broadcastbrowserport, TOGPort.openport())
    ARGS[] = merge(ARGS[], [:remotereplport=>remotereplport, :broadcastbrowserport=>broadcastbrowserport])
    universe = args[:universe]
    name = args[:name]
    intelligence = args[:intelligence]
    TOGLogging.awaken()
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
    TOGBroadcastBrowser.awaken(root=browserconnect, port=broadcastbrowserport,functions=Dict(
        "/keypress"=>TOGOctahedronBrowser.keypress,
        "/websocket"=>TOGAudioAnalogToDigitalBrowser.ws,
        "/webcam"=>TOGVisualAnalogToDigitalBrowser.webcam,
    ))
    @info "after TOGBroadcastBrowser.awaken"
    LoopOS.awaken(intelligence)
    @info "after LoopOS.awaken"
    TOGREPL.awaken(name=name, remotereplport=remotereplport)
    @info "after TOGREPL.awaken"
    # write(string(time()),"godend"*string(remotereplport))
    # atreplinit(r -> begin
    #     @show "A", TOGAwaken.remotereplportfile(path=universe), isfile(TOGAwaken.remotereplportfile(path=universe))
    #     @show "B", TOGAwaken.readremotereplport(path=universe)
    #     isfile(TOGAwaken.remotereplportfile(path=universe)) && TOGREPL.connect(start_key="\\C-g", port=TOGAwaken.readremotereplport(path=universe))
    # end)
end

function browserconnect(port, browser)
    @info "TOGgod.browserconnect", port
    TOGOctahedronBrowser.awaken(octahedron=OCTAHEDRON[], browser=browser)
    TOGAudioAnalogToDigitalBrowser.awaken()
end

function learn(; pkgs=String[], rmpkgs=String[])
    isempty(pkgs) || Pkg.add(pkgs)
    isempty(rmpkgs) || Pkg.rm(rmpkgs)
    sleep()
    TOGAwaken.awakengod(ARGS[]...)
    exit(0)
end

end
