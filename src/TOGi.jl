module TOGi

export LoopOS
import TOGHuman: LoopOS

using TOGAwaken, TOGHuman

function awaken(; universe="..")
    # TOGΩ.awaken(path=universe)
    TOGHuman.awaken(name="i", universe=universe)
    args = merge(NamedTuple(;), [
        # :name => "Dona",
        # :intelligence => TOGXAI.intelligence,
        :pkg => "Dona",
        :path => "../Dona",
        :universe => universe,
    ])
    # TOGAwaken.awakengod(;args...) # DEBUG
end

end
