module Dona

export LoopOS
import TOGgod: LoopOS

export TOGAdvice, TOGPowerOfAttorney, TOGAwaken, TOGOctahedron, TOGCommunicationClient, TOGLearning, TOGCreateClient, TOGBroadcastBrowser, TOGColor, TOGBasicTools
using TOGAdvice, TOGPowerOfAttorney, TOGAwaken, TOGOctahedron, TOGCommunicationClient, TOGLearning, TOGCreateClient, TOGBroadcastBrowser, TOGColor, TOGBasicTools
# export TOGCommunicationClient: Messages

using TOGIntelligence, TOGXAI, TOGgod

intelligence(;
    self,
    history,
    longmemory,
    shortmemory,
    newinput,
    complexity,
) = TOGIntelligence.intelligence(
    self=self,
    history=history,
    longmemory=longmemory,
    shortmemory=shortmemory,
    newinput=newinput,
    complexity=complexity,
    state_post="",
    model=TOGXAI.intelligence)

function awaken(; args...)
    @info "awaken Dona"
    # write(string(time()),"Dona1")
    _args = merge(NamedTuple(;), args)
    _args = merge(_args, [
        :name => "Dona",
        :intelligence => intelligence,
        # :pkg => string(@__MODULE__),
    ])
    # write(string(time()),"Dona2")
    TOGgod.awaken(; _args...)
    # write(string(time()),"Dona3")
end

end
