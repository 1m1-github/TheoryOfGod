"""
julia --optimize=3 --threads=auto --quiet --load tog.jl --eval 'tog.awaken()'
"""
module tog

using Pkg
include(joinpath("src", "TOGAwaken.jl"))

function awaken()
    awakenregistry()
    TOGAwaken.awakenΩ()
    awakenDona()
    # TOGAwaken.awakengod(name="Dona", universe=joinpath(pwd(),"Ω"))
end

function awakenDona()
    TOGAwaken.awakengod(path=joinpath(pwd(),"Dona"), universe=joinpath(pwd(),"Ω"), pkgs=[
        "LoopOS",
        "TOGgod", 
        "TOGXAI", 
        "TOGAdvice", 
        "TOGPowerOfAttorney", 
        "TOGIntelligence", 
        "TOGAwaken", 
        "TOGOctahedron", 
        "TOGCommunicationClient", 
        "TOGLearning", 
        "TOGCreateClient", 
        "TOGBroadcastBrowser", 
        "TOGColor", 
        "TOGBasicTools"
    ])
end

function awakenregistry()
    Pkg.add(["Revise", "LocalRegistry", "TOML", "Git"])
    currentdir = pwd()
    dir = joinpath(DEPOT_PATH[1], "dev")
    isdir(dir) || mkdir(dir)
    files(name) = """["$(joinpath("$currentdir", "src", name * ".jl"))"]"""
    TOGAwaken.julia(
    path=dir,
    project=dir,
    depot=DEPOT_PATH[1],
    code="""
    using Revise, LocalRegistry, TOML, Git
    include(joinpath("$currentdir", "tog.jl"))
    include(joinpath("$currentdir", "src", "TOGAwaken.jl"))
    include(joinpath("$currentdir", "src", "TOGLearning.jl"))
    if !isdir(TOGAwaken.REGISTRYPATH)
        create_registry(TOGAwaken.REGISTRYPATH, TOGAwaken.REGISTRYPATH)
        write(joinpath(TOGAwaken.REGISTRYPATH, TOGLearning.GITIGNOREFILE), TOGLearning.GITIGNORE)
        TOGLearning.addcommit(path=TOGAwaken.REGISTRYPATH, commitmessage=".gitignore")
    end
    TOGLearning.newpkg(name="LoopOS", pkgs=["Pkg", "Serialization"], files=$(files("LoopOS")))
    TOGLearning.newpkg(name="TOGGPU", pkgs=["KernelAbstractions"], files=$(files("TOGGPU")))
    TOGLearning.newpkg(name="TOG∃", pkgs=["KernelAbstractions", "IntervalTrees", "TOGGPU"], files=$(files("TOG∃")))
    TOGLearning.newpkg(name="TOGAwaken", pkgs=["Sockets"], files=$(files("TOGAwaken")))
    TOGLearning.newpkg(name="TOGCaching", pkgs=["LoopOS"], files=$(files("TOGCaching")))
    TOGLearning.newpkg(name="TOGState", pkgs=["LoopOS", "TOGCaching"], files=$(files("TOGState")))
    TOGLearning.newpkg(name="TOGZMQ", pkgs=["ZMQ", "Serialization", "LoopOS", "TOGState"], files=$(files("TOGZMQ")))
    TOGLearning.newpkg(name="TOGCommunicationServer", pkgs=["ZMQ", "LoopOS", "TOGZMQ"], files=$(files("TOGCommunicationServer")))
    TOGLearning.newpkg(name="TOGZMQAPIServer", pkgs=["ZMQ", "Serialization", "LoopOS", "TOGZMQ"], files=$(files("TOGZMQAPIServer")))
    TOGLearning.newpkg(name="TOGObserveServer", pkgs=["ZMQ", "TOG∃", "TOGZMQAPIServer"], files=$(files("TOGObserveServer")))
    TOGLearning.newpkg(name="TOGPrivacy", files=$(files("TOGPrivacy")))
    TOGLearning.newpkg(name="TOGRay", pkgs=["FastGaussQuadrature", "KernelAbstractions", "TOG∃", "TOGGPU"], files=$(files("TOGRay")))
    TOGLearning.newpkg(name="TOGZMQAPIClient", pkgs=["ZMQ", "TOGZMQ", "TOGZMQAPIServer"], files=$(files("TOGZMQAPIClient")))
    TOGLearning.newpkg(name="TOGObserveClient", pkgs=["ZMQ", "TOG∃", "TOGZMQAPIClient"], files=$(files("TOGObserveClient")))
    TOGLearning.newpkg(name="TOGOctahedron", pkgs=["KernelAbstractions", "LinearAlgebra", "LoopOS", "TOG∃", "TOGPrivacy", "TOGRay", "TOGObserveClient", "TOGGPU"], files=$(files("TOGOctahedron")))
    TOGLearning.newpkg(name="TOGCreateServer", pkgs=["ZMQ", "TOG∃", "TOGZMQAPIServer", "TOGOctahedron"], files=$(files("TOGCreateServer")))
    TOGLearning.newpkg(name="TOGLogging", pkgs=["Logging"], files=$(files("TOGLogging")))
    TOGLearning.newpkg(name="TOGPort", pkgs=["Sockets"], files=$(files("TOGPort")))
    TOGLearning.newpkg(name="TOGREPL", pkgs=["Sockets", "RemoteREPL", "ReplMaker", "LoopOS", "TOGAwaken", "TOGPort"], files=$(files("TOGREPL")))
    TOGLearning.newpkg(name="TOGColor", pkgs=["Colors", "TOG∃"], files=$(files("TOGColor")))
    TOGLearning.newpkg(name="TOGMatrix", pkgs=["Adapt", "TOGColor"], files=$(files("TOGMatrix")))
    TOGLearning.newpkg(name="TOGΩ", pkgs=["Serialization", "TOG∃", "TOGAwaken", "Colors", "TOGColor", "TOGMatrix", "TOGCommunicationServer", "TOGObserveServer", "TOGCreateServer", "TOGREPL", "TOGOctahedron", "TOGZMQ", "TOGLogging"], files=$(files("TOGΩ")))
    TOGLearning.newpkg(name="TOGCommunicationClient", pkgs=["ZMQ", "LoopOS", "TOGZMQ"], files=$(files("TOGCommunicationClient")))
    TOGLearning.newpkg(name="TOGLearning", pkgs=["Pkg", "TOML", "LocalRegistry", "Git"], files=$(files("TOGLearning")))
    TOGLearning.newpkg(name="TOGTypst", pkgs=["PNGFiles", "Typst_jll", "TOGMatrix"], files=$(files("TOGTypst")))
    TOGLearning.newpkg(name="TOGCreateClient", pkgs=["ZMQ", "PNGFiles", "TOGTypst", "TOGZMQAPIClient", "TOGOctahedron"], files=$(files("TOGCreateClient")))
    TOGLearning.newpkg(name="TOGBroadcastBrowser", pkgs=["HTTP", "URIs", "Sockets", "LoopOS", "TOGAwaken", "TOGPort"], files=$(files("TOGBroadcastBrowser")))
    TOGLearning.newpkg(name="TOGMoveOctahedron", pkgs=["TOGOctahedron", "TOG∃"], files=$(files("TOGMoveOctahedron")))
    TOGLearning.newpkg(name="TOGOctahedronBrowser", pkgs=["LoopOS", "TOG∃", "TOGBroadcastBrowser", "TOGOctahedron", "TOGColor", "TOGMoveOctahedron"], files=$(files("TOGOctahedronBrowser")))
    TOGLearning.newpkg(name="TOGXAI", pkgs=["HTTP", "JSON3", "Base64", "FileIO", "ImageIO", "Images", "ColorTypes", "TOGLogging"], files=$(files("TOGXAI")))
    TOGLearning.newpkg(name="TOGAudioAnalogToDigitalBrowser", pkgs=["HTTP", "LoopOS", "TOGBroadcastBrowser", "TOGXAI"], files=$(files("TOGAudioAnalogToDigitalBrowser")))
    TOGLearning.newpkg(name="TOGVisualAnalogToDigitalBrowser", pkgs=["ColorTypes", "FixedPointNumbers", "JSON3", "Base64", "PNGFiles", "LoopOS", "TOGBroadcastBrowser"], files=$(files("TOGVisualAnalogToDigitalBrowser")))
    TOGLearning.newpkg(name="TOGTextToAudioBrowser", pkgs=["HTTP", "LoopOS", "TOGBroadcastBrowser", "TOGXAI", "Base64"], files=$(files("TOGTextToAudioBrowser")))
    TOGLearning.newpkg(name="TOGgod", pkgs=["Pkg", "Serialization", "LoopOS", "TOG∃", "TOGAwaken", "TOGCommunicationClient", "TOGOctahedron", "TOGLearning", "TOGObserveClient", "TOGCreateClient", "TOGREPL", "TOGLogging", "TOGBroadcastBrowser", "TOGOctahedronBrowser", "TOGAudioAnalogToDigitalBrowser", "TOGTextToAudioBrowser", "TOGVisualAnalogToDigitalBrowser", "TOGPort", "TOGZMQ"], files=$(files("TOGgod")))
    TOGLearning.newpkg(name="TOGIntelligence", pkgs=["LoopOS", "TOGState", "TOGLogging"], files=$(files("TOGIntelligence")))
    TOGLearning.newpkg(name="TOGIntelligenceLocal", pkgs=["HTTP", "JSON3"], files=$(files("TOGIntelligenceLocal")))
    TOGLearning.newpkg(name="TOGIntelligenceHuman", pkgs=["TOGREPL", "TOGIntelligence"], files=$(files("TOGIntelligenceHuman")))
    TOGLearning.newpkg(name="TOGAnthropic", pkgs=["HTTP", "JSON3"], files=$(files("TOGAnthropic")))
    TOGLearning.newpkg(name="TOGHuman", pkgs=["TOGgod", "TOGIntelligenceHuman"], files=$(files("TOGHuman")))
    TOGLearning.newpkg(name="TOGi", pkgs=["TOGAwaken", "TOGHuman"], files=$(files("TOGi")))
    TOGLearning.newpkg(name="TOGAdvice", files=$(files("TOGAdvice")))
    TOGLearning.newpkg(name="TOGPowerOfAttorney", files=$(files("TOGPowerOfAttorney")))
    TOGLearning.newpkg(name="TOGBasicTools", pkgs=["HTTP", "JSON3", "Base64", "Dates", "SMTPClient", "Serialization", "Gumbo", "Cascadia"], files=$(files("TOGBasicTools")))
    """)
    # TOGLearning.newpkg(name="TOGIntelligenceJanet", pkgs=["LoopOS", "TOGState", "TOGLogging"], files=$(files("TOGIntelligenceJanet")))
    # TOGLearning.newpkg(name="Janet", pkgs=["TOGgod", "TOGXAI", "TOGAdvice", "TOGPowerOfAttorney", "TOGIntelligenceJanet", "TOGAwaken", "TOGOctahedron", "TOGCommunicationClient", "TOGLearning", "TOGCreateClient", "TOGBroadcastBrowser", "TOGColor", "TOGBasicTools"], files=$(files("Janet")))
    # TOGLearning.newpkg(name="Dona", pkgs=["TOGgod", "TOGXAI", "TOGAdvice", "TOGPowerOfAttorney", "TOGIntelligence", "TOGAwaken", "TOGOctahedron", "TOGCommunicationClient", "TOGLearning", "TOGCreateClient", "TOGBroadcastBrowser", "TOGColor", "TOGBasicTools"], files=$(files("Dona")))
end

end
