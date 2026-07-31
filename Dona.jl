module Dona

# export TOGObserveClient, TOGCreateClient, TOGREPL, TOGBroadcastBrowser, TOGOctahedronBrowser, TOGTextToAudioBrowser, TOGVisualAnalogToDigitalBrowser, TOGPort,  TOGOctahedron, TOGColor, RGBA
# export LoopOS, TOGLearning, TOGCommunicationClient, TOGAwaken, TOGLogging, TOGZMQ, TOGAdvice, TOGPowerOfAttorney, TOGBasicTools, TOGgod, TOGMessage
export LoopOS, TOGLearning, TOGCommunicationClient, TOGAwaken, TOGAdvice, TOGPowerOfAttorney, TOGBasicTools, TOGgod, TOGIntelligenceLocal, TOGIntelligence, TOGMessage, TOGBroadcastBrowser
export NAME

# using TOGObserveClient, TOGCreateClient, TOGBroadcastBrowser, TOGOctahedronBrowser, TOGTextToAudioBrowser, TOGVisualAnalogToDigitalBrowser, TOGREPL, TOGPort, TOGZMQ, TOGLogging, TOGOctahedron, TOGColor,
using TOGgod: LoopOS, TOGLearning, TOGCommunicationClient, TOGAwaken, TOGBroadcastBrowser
using TOGAdvice, TOGPowerOfAttorney, TOGBasicTools, TOGgod, TOGIntelligenceLocal, TOGIntelligence
using TOGCommunicationClient: TOGMessage
# using TOGColor: RGBA

const NAME = """
Let your name be Dona.
"""

intelligence(;
    self,
    history,
    longmemory,
    shortmemory,
    newinput,
    complexity,
) = begin
    TOGIntelligence.intelligence(
    self=self,
    history=history,
    longmemory=longmemory,
    shortmemory=shortmemory,
    newinput=newinput,
    complexity=complexity,
    state_post="",
    model=TOGIntelligenceLocal.intelligence)
end

function awaken(; args...)
    _args = merge(NamedTuple(;), args)
    _args = merge(_args, [
        :intelligence => intelligence,
    ])
    TOGgod.awaken(; _args...)
end

end
