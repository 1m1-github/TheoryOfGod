module TOGHuman

export LoopOS
import TOGgod: LoopOS

using TOGgod, TOGIntelligenceHuman

function awaken(; args...)
    # write(string(time()),"human"*string(args))
    # write("TOGHuman", string(args))
    # name = args[:name]
    # universe = args[:universe]
    # TOGgod.awaken(path=pwd(),pkg=string(@__MODULE__), intelligence=TOGIntelligenceHuman.intelligence, name=name, universe=universe)
    _args = merge(NamedTuple(;), args)
    _args = merge(_args, [
        :intelligence=>TOGIntelligenceHuman.intelligence,
        :pkg => string(@__MODULE__),
    ])
    # write(string(time()),"human2"*string(_args))
    TOGgod.awaken(; _args...)
end

end
