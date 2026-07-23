module TOGState

export ignore, unignore

using Revise
using LoopOS
using LoopOS: Loop, LOOP, Action, TrackedSymbol, Input, Peripheral, LISTENING
using TOGCaching
using TOGCaching: cache!
import Base: take!, put!

const IGNORE = Set{Module}([Main, Base, Core, Revise])
ignore(m::Module) = push!(IGNORE, m)
unignore(m::Module) = delete!(IGNORE, m)

# todo add methods to functions based on execution to avoid typing in code

function os_time(timestamp)
    ΔT = timestamp - LOOP.birthtime
    isinf(ΔT) && return "[∞s]"
    "[$(round(Int, ΔT))s]"
end

function expand!(short::Vector{TrackedSymbol})
    modules = [TrackedSymbol(ts.value, sym, f, ts.timestamp) for ts = short if ts.value isa Module for sym = names(ts.value) if sym ≠ ts.sym for f = (getfield(ts.value, sym),) if !(f isa Module)]
    push!(short, modules...)
    _methods = [TrackedSymbol(ts.m, ts.sym, method, ts.timestamp) for ts = short if ts.value isa Function for method = methods(ts.value, ts.m)]
    push!(short, _methods...)
    filter!(ts -> !(ts.value isa Function), short)
end
function dichotomy(f, v)
    in = eltype(v)[]
    out = eltype(v)[]
    for x in v
        push!(f(x) ? in : out, x)
    end
    in, out
end
isoutputperipheralmethod(sym::Symbol, method::Method) = sym == :put! && isperipheralmethod(method)
isinputperipheralmethod(sym::Symbol, method::Method) = sym == :take! && isperipheralmethod(method)
function isperipheralmethod(method::Method)
    t = fieldtypes(method.sig)
    length(t) == 1 && return false
    arg = t[2]
    arg <: Peripheral || arg <: Type{<:Peripheral}
end
function methodperipheraltype(method::Method)
    arg = fieldtypes(method.sig)[2]
    arg <: Peripheral ? arg : arg.parameters[1]
end

function categorize(short::Vector{TrackedSymbol})
    # expand!(short)
    methodshort, short = dichotomy(ts -> ts.value isa Method, short)
    typeshort, short = dichotomy(ts -> ts.value isa Type, short) # todo DataType?
    moduleshort, short = dichotomy(ts -> ts.value isa Module, short)
    outputperipheralmethodshort, methodshort = dichotomy(ts -> isoutputperipheralmethod(ts.sym, ts.value), methodshort)
    inputperipheralmethodshort, methodshort = dichotomy(ts -> isinputperipheralmethod(ts.sym, ts.value), methodshort)
    outputperipheralmethodtypeshort = unique([TrackedSymbol(ts.m, nameof(t), t, ts.timestamp) for ts = outputperipheralmethodshort for t = (methodperipheraltype(ts.value),)])
    inputperipheralmethodtypeshort = unique([TrackedSymbol(ts.m, nameof(t), t, ts.timestamp) for ts = inputperipheralmethodshort for t = (methodperipheraltype(ts.value),)])
    outputperipheralmethodtypeshortvector = map(ts -> ts.value, outputperipheralmethodtypeshort)
    inputperipheralmethodtypeshortvector = map(ts -> ts.value, inputperipheralmethodtypeshort)
    peripheralshort, short = dichotomy(ts -> ts.value isa Peripheral, short)
    outputperipheralshort = filter(ts -> typeof(ts.value) ∈ outputperipheralmethodtypeshortvector, peripheralshort)
    inputperipheralshort = filter(ts -> typeof(ts.value) ∈ inputperipheralmethodtypeshortvector, peripheralshort)
    filter!(ts -> ts.value ∉ outputperipheralmethodtypeshortvector, typeshort)
    filter!(ts -> ts.value ∉ inputperipheralmethodtypeshortvector, typeshort)
    short, outputperipheralshort, inputperipheralshort, methodshort, outputperipheralmethodshort, inputperipheralmethodshort, typeshort, outputperipheralmethodtypeshort, inputperipheralmethodtypeshort, moduleshort
end

function shortsections(short, outputperipheralshort, inputperipheralshort, methodshort, outputperipheralmethodshort, inputperipheralmethodshort, typeshort, outputperipheralmethodtypeshort, inputperipheralmethodtypeshort, moduleshort)
    sections = String[]
    !isempty(short) && push!(sections, state("SHORT MEMORY", short))
    !isempty(short) && push!(sections, state("SHORT MEMORY OUTPUT PERIPHERALS", outputperipheralshort))
    !isempty(short) && push!(sections, state("SHORT MEMORY INTPUT PERIPHERALS", inputperipheralshort))
    !isempty(short) && push!(sections, state("SHORT MEMORY METHODS", methodshort))
    !isempty(short) && push!(sections, state("SHORT MEMORY OUTPUT PERIPHERAL METHODS", outputperipheralmethodshort))
    !isempty(short) && push!(sections, state("SHORT MEMORY INPUT PERIPHERAL METHODS", inputperipheralmethodshort))
    !isempty(short) && push!(sections, state("SHORT MEMORY TYPES", typeshort))
    !isempty(short) && push!(sections, state("SHORT MEMORY OUTPUT PERIPHERAL TYPES", outputperipheralmethodtypeshort))
    !isempty(short) && push!(sections, state("SHORT MEMORY INPUT PERIPHERAL TYPES", inputperipheralmethodtypeshort))
    !isempty(short) && push!(sections, state("SHORT MEMORY MODULES", moduleshort))
    sections
end

function state(
    self::String,
    history::Vector{Action},
    long::Vector{String},
    short::Vector{TrackedSymbol},
    input::Vector{Input},
    state_post::String,
)
    filter!(ts -> ts.value ∉ IGNORE, short)
    expand!(short)
    # cached=short # DEBUG
    # volatile = TrackedSymbol[] # DEBUG
    cached, volatile = cache!(short)
    # short, outputperipheralshort, inputperipheralshort, methodshort, outputperipheralmethodshort, inputperipheralmethodshort, typeshort, outputperipheralmethodtypeshort, inputperipheralmethodtypeshort, moduleshort = categorize(short)
    # cached, volatile = cache!(short)
    # cached, volatile = cache!(outputperipheralshort)
    # cached, volatile = cache!(inputperipheralshort)
    # cached, volatile = cache!(methodshort)
    # cached, volatile = cache!(outputperipheralmethodshort)
    # cached, volatile = cache!(short)
    # cached, volatile = cache!(short)
    # cached, volatile = cache!(short)
    # cached, volatile = cache!(short)
    cached, outputperipheralcached, inputperipheralcached, methodcached, outputperipheralmethodcached, inputperipheralmethodcached, typecached, outputperipheralmethodtypecached, inputperipheralmethodtypecached, modulecached = categorize(cached)
    volatile, outputperipheralvolatile, inputperipheralvolatile, methodvolatile, outputperipheralmethodvolatile, inputperipheralmethodvolatile, typevolatile, outputperipheralmethodtypevolatile, inputperipheralmethodtypevolatile, modulevolatile = categorize(volatile)
    historyvolatile = TrackedSymbol[]
    for (i, action) = enumerate(history)
        if istaskfailed(action.task)
            # length(action.input) == 1 && only(action.input).source isa Loop && continue
            push!(historyvolatile, TrackedSymbol(LoopOS, Symbol("history[][$i].input"), action.input, action.timestamp))
            push!(historyvolatile, TrackedSymbol(LoopOS, Symbol("history[][$i].output"), action.output, action.timestamp))
            push!(historyvolatile, TrackedSymbol(LoopOS, Symbol("history[][$i].task"), action.task, action.timestamp))
        end
    end
    push!(historyvolatile, TrackedSymbol(LoopOS, :LOOP, LOOP, Inf))
    sectionscached = shortsections(cached, outputperipheralcached, inputperipheralcached, methodcached, outputperipheralmethodcached, inputperipheralmethodcached, typecached, outputperipheralmethodtypecached, inputperipheralmethodtypecached, modulecached)
    sectionsvolatile = shortsections(volatile, outputperipheralvolatile, inputperipheralvolatile, methodvolatile, outputperipheralmethodvolatile, inputperipheralmethodvolatile, typevolatile, outputperipheralmethodtypevolatile, inputperipheralmethodtypevolatile, modulevolatile)
    cached_sections = [self, sectionscached...]
    volatile_section = [
        state("LONG MEMORY", long),
        sectionsvolatile...,
        state("HISTORY", historyvolatile),
        state("INPUTS", input),
        state_post]
    join(cached_sections, "\n\n"), join(volatile_section, "\n\n")
end
state(x) = string(x) # Use `dump` if you need to see more of anything but careful, it could be a lot
state(description::String, value::Vector) = isempty(value) ? "" : description * " === BEGIN" * "\n\n" * state(value) * "\n\n" * description * " === END"
state(T::DataType) = strip(sprint(dump, T)) * " end"
state(r::Ref) = state(r[])
state(m::Module) = string(nameof(m)) * "::Module"
state(s::String) = "\"$s\""
state(v::Vector) = "[" * join(state.(v), ",\n") * "]"
state(v::Vector{T}) where T<:Number = "[" * join(string.(v), ", ") * "]"
state(i::Input) = "LoopOS.Input($(os_time(i.timestamp)), $(state(i.source)), $(state(i.input)))"
state(p::Peripheral) = state(typeof(p)) * "(" * (p ∈ LISTENING ? "" : "not") * "listening)"
function state(a::Action)
    _state = "inputs=$(state(a.input))"
    _state *= "\n$(state(a.task))"
    istaskfailed(a.task) && (_state *= "\noutput=$(a.output)")
    _state
end
# state(::Loop) = "LoopOS.LOOP"
function state(t::Task)
    _state = ["$(repr(f)):$(f(t))" for f = [istaskstarted, istaskdone, istaskfailed]]
    exception = istaskfailed(t) ? ",exception:$(state(t.exception))" : ""
    "Task(" * join(_state, ",") * exception * ")"
end
function state(e::Exception)
    e isa TaskFailedException && return state(e.task.exception)
    string(e) * string(catch_backtrace())
end
function docs(method::Method, multidoc::Docs.MultiDoc)
    sig = method.sig
    sig isa UnionAll && (sig = Base.unwrap_unionall(sig))
    params = sig.parameters[2:end]
    paramstuple = Tuple{params...}
    if haskey(multidoc.docs, paramstuple)
        multidoc.docs[paramstuple].text
    else
        for (t, s) in multidoc.docs
            if paramstuple <: t
                return s.text
            end
        end
        ""
    end
end
function docs(ts::TrackedSymbol)
    binding = Docs.Binding(ts.m, ts.sym)
    meta = Docs.meta(ts.m)
    haskey(meta, binding) || return ""
    multidoc = meta[binding]
    d = length(multidoc.docs) == 1 ? only(multidoc.docs)[2].text : docs(ts.value, multidoc)
    isempty(d) ? "" : "\"" * strip(string(join(d))) * "\"\n"
end
function state(method::Method)
    sig_str = split(string(method), " @")[1]
    replace(sig_str, "__source__::LineNumberNode, __module__::Module, " => "")
end
function state(v::TrackedSymbol)
    value = v.value
    T = typeof(value)
    ref = if T <: Ref
        ref = "[]"
        value = v.value[]
        T = typeof(value)
    else
        ""
    end
    T_str = T ∈ [Type, Method] ? "" : "::" * string(T)
    # @info v, T
    _sizeofvalue = T ∈ [DataType, Type, Method, Module, UnionAll] ? 0 : sizeof(value)
    sizeofvalue = iszero(_sizeofvalue) ? "" : "(sizeof=" * string(sizeof(value)) * ")"
    _state = if value === LOOP && isinf(v.timestamp)
        _s = strip(sprint(dump, value))
        replace(_s, r": (\w+) " => s"::\1=") * " end"
    else
        state(value)
    end
    m = v.m == Main ? "" : string(v.m) * "."
    os_time(v.timestamp) * "\n" * docs(v) * m * string(v.sym) * ref * T_str * sizeofvalue * "=" * _state
end
function state(_state::Vector{TrackedSymbol})
    sort!(_state, lt=(s, _s) -> s.timestamp == _s.timestamp ? s.value isa Action : s.timestamp < _s.timestamp)
    replace(join(filter(!isempty, state.(_state)), '\n'), "Main." => "")
end

end
