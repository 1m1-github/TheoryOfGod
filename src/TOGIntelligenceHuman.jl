module TOGIntelligenceHuman

using TOGREPL, TOGIntelligence

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
    model=_intelligence)

_intelligence = complexity -> __intelligence
__intelligence(
    input_system,
    input_user,
    max_output_tokens,
    temperature,
) = begin
    # @show "TOGIntelligenceHuman.intelligence", input_system, input_user, max_output_tokens, temperature
    sleep(5)
    # output = take!(TOGREPL.repl)
    output = """nothing"""
    # @show "TOGIntelligenceHuman.intelligence", output
    output, zero(typeof(temperature))
end

end
