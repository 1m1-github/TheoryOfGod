module TOGgod

export learn, OCTAHEDRON, put!
export LoopOS, TOGObserveClient, TOGCreateClient, TOGLearning, TOGCommunicationClient, TOGAwaken, TOGLogging, TOGREPL, TOGBroadcastBrowser, TOGOctahedronBrowser, TOGAudioAnalogToDigitalBrowser, TOGTextToAudioBrowser, TOGVisualAnalogToDigitalBrowser, TOGPort, TOGZMQ
# export T

using Pkg, Serialization
using LoopOS, TOGObserveClient, TOGCreateClient, TOGLearning, TOGCommunicationClient, TOGAwaken, TOGLogging, TOGREPL, TOGBroadcastBrowser, TOGOctahedronBrowser, TOGAudioAnalogToDigitalBrowser, TOGTextToAudioBrowser, TOGVisualAnalogToDigitalBrowser, TOGPort, TOGZMQ
using LoopOS: Peripheral
using TOGOctahedron: Octahedron
using TOGExist: ○
import Base: put!

# const T = Ref{DataType}()
const OCTAHEDRON = Ref{Octahedron}()
const ARGS = Ref{NamedTuple}((;))

__init__() = atexit((n)->sleep(n))

struct Ω <: Peripheral end
"""
Create inside your main octahedron. See the other `put!` methods for possible args.
"""
put!(::Type{Ω}, args...) = TOGCreateClient.put!(OCTAHEDRON[], args...)

function sleep(exitcode)
    # @info "TOGgod, sleep", exitcode
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
    # T[] = TOGObserveClient.togtype()
    ϕ = MathConstants.golden
    OCTAHEDRON[] = Octahedron(
        t=TOGObserveClient.togtime(),
        d=[ϕ^-4, ϕ^-3, ϕ^-2, ϕ^-1],
        observer=[0, ○, ○, ○],
        focus=[0, ○, ○, ○+0.1],
        ρ=[0, 0.1, 0.1, 0],
        ♯=(1, 1))
    TOGBroadcastBrowser.awaken(root=browserconnect, port=broadcastbrowserport, functions=Dict(
        "/keypress"=>TOGOctahedronBrowser.keypress,
        # "/websocket"=>TOGAudioAnalogToDigitalBrowser.ws,
        "/audio"=>TOGAudioAnalogToDigitalBrowser.audio,
        "/webcam"=>TOGVisualAnalogToDigitalBrowser.webcam,
    ))
    LoopOS.awaken(intelligence)
    # TOGREPL.awaken(name=name, remotereplport=remotereplport)
    isinteractive() ? nothing : wait(Condition())
end

function browserconnect(port, browser)
    # @info "TOGgod.browserconnect", port
    TOGOctahedronBrowser.awaken(octahedron=OCTAHEDRON[], browser=browser)
    TOGAudioAnalogToDigitalBrowser.awaken()
end

"""
Learn by adding or removing pkg for yourself.
To `learn`, the process is to create a new `Pkg` with `newpkg` after writing the corresponding source code file and then running `learn`, which will `add`/`rm` `Pkg`s and restart your Julia VM, which makes you sleep for a few seconds until you awake again with abilities added/removed.
Your abilities are thus listed in the `.tog/Project.toml`.
# Arguments
- `addpkg::Vector{String}`: `Pkg`s to be added.
- `rmpkg::Vector{String}`: `Pkg`s to be removed.
"""
function learn(; addpkg::Vector{String}=String[], rmpkg::Vector{String}=String[])
    # @info "TOGgod, learn"
    isempty(addpkg) || Pkg.add(addpkg)
    isempty(rmpkg) || Pkg.rm(rmpkg)
    sleep(0)
    # @info "TOGgod, learn, after sleep"
    TOGAwaken.awakengod(; ARGS[]...)
    # @info "TOGgod, learn, after awakengod"
    exit(0)
end

end
