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
    # awakenJanet()
    # TOGAwaken.awakengod(name="Dona", universe=joinpath(pwd(),"Ω"))
end

function awakenDona()
    TOGAwaken.awakengod(path=joinpath(pwd(),"Dona"), universe=joinpath(pwd(),"Ω"), pkg=[
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
        "TOGBasicTools",
        "TOGTextToAudioBrowser",
        "TOGVisualAnalogToDigitalBrowser",
        "TOGColor",
    ])
end

function awakenJanet()
    TOGAwaken.awakengod(path=joinpath(pwd(),"Janet"), universe=joinpath(pwd(),"Ω"), pkg=[
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
        "TOGBasicTools",
        "TOGTextToAudioBrowser",
        "TOGVisualAnalogToDigitalBrowser",
        "TOGColor",
    ])
end

function awakenregistry()
    Pkg.add(["Revise", "LocalRegistry", "TOML", "Git"])
    currentdir = pwd()
    dir = joinpath(DEPOT_PATH[1], "dev")
    isdir(dir) || mkdir(dir)
    file(name) = """["$(joinpath("$currentdir", "src", name * ".jl"))"]"""
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
    TOGLearning.newpkg(name="LoopOS", pkg=["Pkg", "Serialization"], file=$(file("LoopOS")))
    TOGLearning.newpkg(name="TOGGPU", pkg=["KernelAbstractions"], file=$(file("TOGGPU")))
    TOGLearning.newpkg(name="TOGExist", pkg=["KernelAbstractions", "IntervalTrees", "TOGGPU"], file=$(file("TOGExist")))
    TOGLearning.newpkg(name="TOGAwaken", pkg=["Sockets"], file=$(file("TOGAwaken")))
    TOGLearning.newpkg(name="TOGCaching", pkg=["LoopOS"], file=$(file("TOGCaching")))
    TOGLearning.newpkg(name="TOGState", pkg=["Revise", "LoopOS", "TOGCaching"], file=$(file("TOGState")))
    TOGLearning.newpkg(name="TOGZMQ", pkg=["ZMQ", "Serialization", "LoopOS", "TOGState"], file=$(file("TOGZMQ")))
    TOGLearning.newpkg(name="TOGCommunicationServer", pkg=["ZMQ", "LoopOS", "TOGZMQ"], file=$(file("TOGCommunicationServer")))
    TOGLearning.newpkg(name="TOGZMQAPIServer", pkg=["ZMQ", "Serialization", "LoopOS", "TOGZMQ", "TOGLogging"], file=$(file("TOGZMQAPIServer")))
    TOGLearning.newpkg(name="TOGObserveServer", pkg=["ZMQ", "TOGExist", "TOGZMQAPIServer"], file=$(file("TOGObserveServer")))
    TOGLearning.newpkg(name="TOGPrivacy", file=$(file("TOGPrivacy")))
    TOGLearning.newpkg(name="TOGRay", pkg=["FastGaussQuadrature", "KernelAbstractions", "TOGExist", "TOGGPU"], file=$(file("TOGRay")))
    TOGLearning.newpkg(name="TOGZMQAPIClient", pkg=["ZMQ", "TOGZMQ", "TOGZMQAPIServer"], file=$(file("TOGZMQAPIClient")))
    TOGLearning.newpkg(name="TOGObserveClient", pkg=["ZMQ", "TOGExist", "TOGZMQAPIClient"], file=$(file("TOGObserveClient")))
    TOGLearning.newpkg(name="TOGOctahedron", pkg=["KernelAbstractions", "LinearAlgebra", "Colors", "FileIO", "ImageIO", "LoopOS", "TOGColor", "TOGXAI", "TOGExist", "TOGPrivacy", "TOGRay", "TOGObserveClient", "TOGGPU"], file=$(file("TOGOctahedron")))
    TOGLearning.newpkg(name="TOGCreateServer", pkg=["ZMQ", "TOGExist", "TOGZMQAPIServer"], file=$(file("TOGCreateServer")))
    TOGLearning.newpkg(name="TOGLogging", pkg=["Logging"], file=$(file("TOGLogging")))
    TOGLearning.newpkg(name="TOGPort", pkg=["Sockets"], file=$(file("TOGPort")))
    TOGLearning.newpkg(name="TOGREPL", pkg=["Sockets", "RemoteREPL", "ReplMaker", "LoopOS", "TOGAwaken", "TOGPort"], file=$(file("TOGREPL")))
    TOGLearning.newpkg(name="TOGColor", pkg=["Colors", "ColorTypes", "TOGExist"], file=$(file("TOGColor")))
    TOGLearning.newpkg(name="TOGMatrix", pkg=["Adapt", "TOGColor"], file=$(file("TOGMatrix")))
    TOGLearning.newpkg(name="TOGOmega", pkg=["Serialization", "Colors", "ColorTypes", "TOGExist", "TOGAwaken", "TOGColor", "TOGMatrix", "TOGCommunicationServer", "TOGObserveServer", "TOGCreateServer", "TOGREPL", "TOGOctahedron", "TOGZMQ", "TOGLogging"], file=$(file("TOGOmega")))
    TOGLearning.newpkg(name="TOGCommunicationClient", pkg=["ZMQ", "LoopOS", "TOGZMQ"], file=$(file("TOGCommunicationClient")))
    TOGLearning.newpkg(name="TOGLearning", pkg=["Pkg", "TOML", "LocalRegistry", "Git"], file=$(file("TOGLearning")))
    TOGLearning.newpkg(name="TOGTypst", pkg=["PNGfile", "Typst_jll", "TOGMatrix"], file=$(file("TOGTypst")))
    TOGLearning.newpkg(name="TOGCreateClient", pkg=["ZMQ", "PNGfile", "TOGExist", "TOGMatrix", "TOGTypst", "TOGZMQAPIClient", "TOGOctahedron"], file=$(file("TOGCreateClient")))
    TOGLearning.newpkg(name="TOGBroadcastBrowser", pkg=["HTTP", "URIs", "Sockets", "LoopOS", "TOGAwaken", "TOGPort"], file=$(file("TOGBroadcastBrowser")))
    TOGLearning.newpkg(name="TOGMoveOctahedron", pkg=["TOGOctahedron", "TOGExist"], file=$(file("TOGMoveOctahedron")))
    TOGLearning.newpkg(name="TOGOctahedronBrowser", pkg=["Colors", "LoopOS", "TOGExist", "TOGBroadcastBrowser", "TOGOctahedron", "TOGColor", "TOGMoveOctahedron"], file=$(file("TOGOctahedronBrowser")))
    TOGLearning.newpkg(name="TOGXAI", pkg=["HTTP", "JSON3", "Base64", "FileIO", "ImageIO", "Images", "ColorTypes", "TOGLogging"], file=$(file("TOGXAI")))
    TOGLearning.newpkg(name="TOGAudioAnalogToDigitalBrowser", pkg=["HTTP", "LoopOS", "TOGBroadcastBrowser", "TOGXAI"], file=$(file("TOGAudioAnalogToDigitalBrowser")))
    TOGLearning.newpkg(name="TOGVisualAnalogToDigitalBrowser", pkg=["ColorTypes", "FixedPointNumbers", "JSON3", "Base64", "PNGfile", "LoopOS", "TOGXAI", "TOGBroadcastBrowser"], file=$(file("TOGVisualAnalogToDigitalBrowser")))
    TOGLearning.newpkg(name="TOGTextToAudioBrowser", pkg=["HTTP", "LoopOS", "TOGBroadcastBrowser", "TOGXAI", "Base64"], file=$(file("TOGTextToAudioBrowser")))
    TOGLearning.newpkg(name="TOGgod", pkg=["Pkg", "Serialization", "LoopOS", "TOGExist", "TOGAwaken", "TOGCommunicationClient", "TOGOctahedron", "TOGLearning", "TOGObserveClient", "TOGCreateClient", "TOGREPL", "TOGLogging", "TOGBroadcastBrowser", "TOGOctahedronBrowser", "TOGAudioAnalogToDigitalBrowser", "TOGTextToAudioBrowser", "TOGVisualAnalogToDigitalBrowser", "TOGPort", "TOGZMQ"], file=$(file("TOGgod")))
    TOGLearning.newpkg(name="TOGIntelligence", pkg=["LoopOS", "TOGState", "TOGLogging"], file=$(file("TOGIntelligence")))
    TOGLearning.newpkg(name="TOGIntelligenceLocal", pkg=["HTTP", "JSON3"], file=$(file("TOGIntelligenceLocal")))
    TOGLearning.newpkg(name="TOGIntelligenceHuman", pkg=["TOGREPL", "TOGIntelligence"], file=$(file("TOGIntelligenceHuman")))
    TOGLearning.newpkg(name="TOGAnthropic", pkg=["HTTP", "JSON3"], file=$(file("TOGAnthropic")))
    TOGLearning.newpkg(name="TOGHuman", pkg=["TOGgod", "TOGIntelligenceHuman"], file=$(file("TOGHuman")))
    TOGLearning.newpkg(name="TOGi", pkg=["TOGAwaken", "TOGHuman"], file=$(file("TOGi")))
    TOGLearning.newpkg(name="TOGAdvice", file=$(file("TOGAdvice")))
    TOGLearning.newpkg(name="TOGPowerOfAttorney", file=$(file("TOGPowerOfAttorney")))
    TOGLearning.newpkg(name="TOGBasicTools", pkg=["HTTP", "JSON3", "Base64", "Dates", "SMTPClient", "Serialization", "Gumbo", "Cascadia"], file=$(file("TOGBasicTools")))
    """)
end

end
