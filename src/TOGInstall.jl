"""

"""
module TOGInstall

using Pkg, Revise, Git, LocalRegistry, TOML
include(joinpath("src", "TOGLearning.jl"))

function newpkg(; name, pkgs=String[])
    TOGLearning.newpkg(
        name=name,
        files=[joinpath("src", name * ".jl")],
        pkgs=pkgs,
        mvfiles=true)
end

function awaken()
    Pkg.add(["Revise", "LocalRegistry"])

    registrypath=joinpath(DEPOT_PATH[1], "registries", "TOGRegistry")
    create_registry(registrypath, registrypath)
    write(joinpath(DEPOT_PATH[1], "registries", "TOGRegistry", TOGLearning.GITIGNOREFILE), TOGLearning.GITIGNORE)
    TOGLearning.addcommit(path=joinpath(DEPOT_PATH[1], "registries", "TOGRegistry"), commitmessage=".gitignore")

    newpkg(name="LoopOS", pkgs=["Pkg", "Serialization"])
    # newpkg(name="TOGGPU", pkgs=["KernelAbstractions"])
    # newpkg(name="TOG", pkgs=["KernelAbstractions", "IntervalTrees", "TOGGPU"])
    # newpkg(name="TOGAwaken", pkgs=["Sockets"])
    # newpkg(name="TOGZMQ", pkgs=["ZMQ", "Serialization", "LoopOS"])
    # newpkg(name="TOGCommunicationServer", pkgs=["ZMQ", "LoopOS", "TOGZMQ"])
    # newpkg(name="TOGZMQAPIServer", pkgs=["ZMQ", "Serialization", "LoopOS", "TOGZMQ"])
    # newpkg(name="TOGObserveServer", pkgs=["ZMQ", "TOG", "TOGZMQAPIServer"])
    # newpkg(name="TOGPrivacy")
    # newpkg(name="TOGRay", pkgs=["FastGaussQuadrature", "KernelAbstractions", "TOG", "TOGGPU"])
    # newpkg(name="TOGZMQAPIClient", pkgs=["ZMQ", "TOGZMQ", "TOGZMQAPIServer"])
    # newpkg(name="TOGObserveClient", pkgs=["ZMQ", "TOG", "TOGZMQAPIClient"])
    # newpkg(name="TOGOctahedron", pkgs=["KernelAbstractions", "LinearAlgebra", "LoopOS", "TOG", "TOGPrivacy", "TOGRay", "TOGObserveClient", "TOGGPU"])
    # newpkg(name="TOGCreateServer", pkgs=["ZMQ", "TOG", "TOGZMQAPIServer", "TOGOctahedron"])
    # newpkg(name="TOGLogging", pkgs=["Logging"])
    # newpkg(name="TOGREPL", pkgs=["Sockets", "RemoteREPL", "ReplMaker", "LoopOS", "TOGAwaken"])
    # newpkg(name="TOGColor", pkgs=["Colors", "TOG"])
    # newpkg(name="TOGMatrix", pkgs=["Adapt", "TOGColor"])
    # newpkg(name="TOGΩ", pkgs=["Serialization", "TOG", "TOGAwaken", "Colors", "TOGColor", "TOGMatrix", "TOGCommunicationServer", "TOGObserveServer", "TOGCreateServer", "TOGREPL", "TOGOctahedron", "TOGLogging"])
    # newpkg(name="TOGCommunicationClient", pkgs=["ZMQ", "LoopOS", "TOGZMQ"])
    # newpkg(name="TOGLearning", pkgs=["Pkg", "TOML", "LocalRegistry", "Git"])
    # newpkg(name="TOGTypst", pkgs=["PNGFiles", "Typst_jll", "TOGMatrix"])
    # newpkg(name="TOGCreateClient", pkgs=["ZMQ", "PNGFiles", "TOGTypst", "TOGZMQAPIClient", "TOGOctahedron"])
    # newpkg(name="TOGBroadcastBrowser", pkgs=["HTTP", "URIs", "Sockets", "LoopOS", "TOGAwaken"])
    # newpkg(name="TOGMoveOctahedron", pkgs=["TOGOctahedron", "TOG"])
    # newpkg(name="TOGOctahedronBrowser", pkgs=["LoopOS", "TOG", "TOGBroadcastBrowser", "TOGOctahedron", "TOGColor", "TOGMoveOctahedron"])
    # newpkg(name="TOGXAI", pkgs=["HTTP", "JSON3", "Base64", "FileIO", "ImageIO", "Images", "ColorTypes", "TOGLogging"])
    # newpkg(name="TOGAudioAnalogToDigitalBrowser", pkgs=["HTTP", "LoopOS", "TOGBroadcastBrowser", "TOGXAI"])
    # newpkg(name="TOGVisualAnalogToDigitalBrowser", pkgs=["ColorTypes", "FixedPointNumbers", "JSON3", "Base64", "PNGFiles", "LoopOS", "TOGBroadcastBrowser"])
    # newpkg(name="TOGTextToAudioBrowser", pkgs=["HTTP", "LoopOS", "TOGBroadcastBrowser", "TOGXAI", "Base64"])
    # newpkg(name="TOGgod", pkgs=["Pkg", "Serialization", "LoopOS", "TOG", "TOGAwaken", "TOGCommunicationClient", "TOGOctahedron", "TOGLearning", "TOGObserveClient", "TOGCreateClient", "TOGREPL", "TOGLogging", "TOGBroadcastBrowser", "TOGOctahedronBrowser", "TOGAudioAnalogToDigitalBrowser", "TOGTextToAudioBrowser", "TOGVisualAnalogToDigitalBrowser"])
    # newpkg(name="TOGCaching", pkgs=["LoopOS"])
    # newpkg(name="TOGState", pkgs=["LoopOS", "TOGCaching"])
    # newpkg(name="TOGIntelligence", pkgs=["LoopOS", "TOGState", "TOGLogging"])
    # newpkg(name="TOGIntelligenceLocal", pkgs=["HTTP", "JSON3"])
    # newpkg(name="TOGIntelligenceHuman", pkgs=["TOGREPL", "TOGIntelligence"])
    # newpkg(name="TOGAnthropic", pkgs=["HTTP", "JSON3"])
    # newpkg(name="TOGHuman", pkgs=["TOGgod", "TOGIntelligenceHuman"])
    # newpkg(name="TOGi", pkgs=["TOGAwaken", "TOGHuman"])
    # newpkg(name="TOGAdvice")
    # newpkg(name="TOGPowerOfAttorney")
    # newpkg(name="TOGBasicTools", pkgs=["HTTP", "JSON3", "Base64", "Dates", "SMTPClient", "Serialization", "Gumbo", "Cascadia"])
    # newpkg(name="Dona", pkgs=["TOGgod", "TOGXAI", "TOGAdvice", "TOGPowerOfAttorney", "TOGIntelligence", "TOGAwaken", "TOGOctahedron", "TOGCommunicationClient", "TOGLearning", "TOGCreateClient", "TOGBroadcastBrowser", "TOGColor", "TOGBasicTools"])
end

end
