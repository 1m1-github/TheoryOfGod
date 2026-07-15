module TOGState

using LoopOS
using LoopOS: Loop, LOOP, Action, TrackedSymbol, Input, Peripheral
using TOGCaching: cache!

# todo add methods to functions based on execution to avoid typing in code

function os_time(timestamp)
    ΔT = timestamp - LOOP.boottime
    isinf(ΔT) && return "[∞s]"
    "[$(round(Int, ΔT))s]"
end

function state(
    self::String,
    history::Vector{Action},
    long_memory::Vector{String},
    short_memory::Vector{TrackedSymbol},
    input::Vector{Input},
    state_post::String,
)
    outputperipheral = [t.value for t in short_memory if t.value isa Type && hasmethod(Base.put!, (t.value, Vararg{Any}))]
    cached, volatile = cache!(short_memory)
    cached = filter(c -> !(c.m === Main && c.sym == :Main && c.value === Main), cached)
    for (i, action) = enumerate(history)
        push!(volatile, TrackedSymbol(LoopOS, Symbol("history[][$i].input"), action.input, action.timestamp))
        if istaskfailed(action.task)
            push!(volatile, TrackedSymbol(LoopOS, Symbol("history[][$i].task"), action.task, action.timestamp))
            push!(volatile, TrackedSymbol(LoopOS, Symbol("history[][$i].output"), action.output, action.timestamp))
        end
    end
    push!(volatile, TrackedSymbol(LoopOS, :LOOP, LOOP, Inf))
    cached_sections = [self, state("SHORT MEMORY", cached)]
    volatile_section = [state("LONG MEMORY", long_memory), state("HISTORY ∪ SHORT MEMORY", volatile), state("OUTPUT PERIPHERALS", outputperipheral), state("INPUTS", input), state_post]
    join(cached_sections, "\n\n"), join(volatile_section, "\n\n")
end
state(x) = string(x) # Use `dump` if you need to see more of anything but careful, it could be a lot
state(description::String, value::Any) = description * " === BEGIN" * "\n\n" * state(value) *  "\n\n" * description * " === END"
state(T::DataType) = strip(sprint(dump, T)) * " end"
state(r::Ref) = state(r[])
function state(timestamp::Float64, m::Module)
    name = string(nameof(m))
    # if m ∉ MODULES
        m ∈ [Base, Core] && return ""
        # return os_time(timestamp) * name * "::Module (`export`ed symbols not shown, use `add_module_to_state` if you need to)"
    # end
    _state = String[]
    for name = names(m)
        f = getfield(m, name)
        f isa Module && continue
        push!(_state, state(TrackedSymbol(m, name, f, timestamp)))
    end
    join(_state, '\n')
end
state(s::String) = "\"$s\""
state(v::Vector) = "[" * join(state.(v), ",\n") * "]"
state(v::Vector{T}) where T <: Number = "[" * join(string.(v), ", ") * "]"
state(i::Input) = "LoopOS.Input($(os_time(i.timestamp)), $(state(i.source)), $(state(i.input)))"
state(p::Peripheral) = state(typeof(p))
function state(a::Action)
    _state = "inputs=$(state(a.input))"
    _state *= "\n$(state(a.task))"
    istaskfailed(a.task) && ( _state *= "\noutput=$(a.output)" )
    _state
end
state(::Loop) = "LoopOS.LOOP"
function state(t::Task)
    _state = ["$(repr(f)):$(f(t))" for f = [istaskstarted, istaskdone, istaskfailed]]
    exception = istaskfailed(t) ? ",exception:$(state(t.exception))" : ""
    "Task(" * join(_state, ",") * exception * ")"
end
function state(x::Exception)
    x isa TaskFailedException && return state(x.task.exception)
    sprint(showerror, x)
end
function state(method::Method)
    sig = method.sig
    sig isa UnionAll && (sig = Base.unwrap_unionall(sig))
    params = sig.parameters[2:end]
    m = method.module
    sig_str = split(string(method), " @")[1]
    sig_str = replace(sig_str, "__source__::LineNumberNode, __module__::Module, " => "")
    binding = Docs.Binding(m, method.name)
    meta = Docs.meta(m)
    haskey(meta, binding) || return ""
    multidoc = meta[binding]
    paramstuple = Tuple{params...}
    docs = []
    if haskey(multidoc.docs, paramstuple)
        docs = multidoc.docs[paramstuple].text
    else
        for (t, s) in multidoc.docs
            if paramstuple <: t
                docs = s.text
                break
            end
        end
    end
    isempty(docs) && return ""
    doc_str = "\"" * strip(string(join(docs))) * "\" "
    doc_str * sig_str
end
function state(v::TrackedSymbol)
    # v.m ∉ MODULES && return ""
    m = v.m == Main ? "" : string(v.m) * "."
    value = v.value
    T = typeof(value)
    ref = ""
    if T <: Ref
        ref = "[]"
        value = v.value[]
        T = typeof(value)
    end
    if T <: Function
        return join(state.([TrackedSymbol(v.m, v.sym, method, v.timestamp) for method = methods(value, v.m)]), '\n')
    elseif T <: Method
        return os_time(v.timestamp) * state(value)
    end
    T_str = T ∈ [DataType, Method] ? "" : string(T)
    _sizeofvalue = value isa Type ? 0 : sizeof(value)
    sizeofvalue = iszero(_sizeofvalue) ? "" : "(sizeof=" * string(sizeof(value)) * ")"
    if T == Module
        state(v.timestamp, value)
    else
        _state = if value === LOOP && isinf(v.timestamp)
            _s = strip(sprint(dump, value))
            replace(_s, r": (\w+) " => s"::\1=") * " end"
        else
            state(value)
        end
        os_time(v.timestamp) * m * string(v.sym) * ref * "::" * T_str * sizeofvalue * "=" * _state
    end
end
function state(_state::Vector{TrackedSymbol})
    sort!(_state, lt = (s, _s) -> s.timestamp == _s.timestamp ? s.value isa Action : s.timestamp < _s.timestamp)
    replace(join(filter(!isempty, state.(_state)), '\n'), "Main." => "")
end

end
