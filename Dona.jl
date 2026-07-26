module Dona

export LoopOS, TOGObserveClient, TOGCreateClient, TOGLearning, TOGCommunicationClient, TOGAwaken, TOGLogging, TOGREPL, TOGBroadcastBrowser, TOGOctahedronBrowser, TOGAudioAnalogToDigitalBrowser, TOGTextToAudioBrowser, TOGVisualAnalogToDigitalBrowser, TOGPort, TOGZMQ
export TOGAdvice, TOGPowerOfAttorney, TOGOctahedron, TOGColor, TOGBasicTools, TOGgod, TOGXAI
export TOGMessage, RGBA

using TOGgod: LoopOS, TOGObserveClient, TOGCreateClient, TOGLearning, TOGCommunicationClient, TOGAwaken, TOGLogging, TOGREPL, TOGBroadcastBrowser, TOGOctahedronBrowser, TOGAudioAnalogToDigitalBrowser, TOGTextToAudioBrowser, TOGVisualAnalogToDigitalBrowser, TOGPort, TOGZMQ
using TOGAdvice, TOGPowerOfAttorney, TOGOctahedron, TOGColor, TOGBasicTools, TOGgod, TOGXAI, TOGIntelligence
using TOGCommunicationClient: TOGMessage
using TOGColor: RGBA

intelligence(;
    self,
    history,
    longmemory,
    shortmemory,
    newinput,
    complexity,
) = begin
    @info "Dona, intelligence"
    TOGIntelligence.intelligence(
    self=self,
    history=history,
    longmemory=longmemory,
    shortmemory=shortmemory,
    newinput=newinput,
    complexity=complexity,
    state_post="",
    model=TOGXAI.intelligence)
end

function awaken(; args...)
    _args = merge(NamedTuple(;), args)
    _args = merge(_args, [
        :intelligence => intelligence,
    ])
    TOGgod.awaken(; _args...)
end

end
