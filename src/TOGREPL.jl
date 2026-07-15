module TOGREPL

# using REPL
using ReplMaker, RemoteREPL, Sockets
using TOGAwaken, TOGPort
# using LoopOS: listen, Peripheral
# import Base.take!, Base.put!

# struct repl <: Peripheral
#     c::Channel{String}
# end
# const REPLINSTANCE = repl(Channel{String}(Inf))
# take!(::Type{repl}) = begin
#     # @show "TOGREPL.take!", REPLINSTANCE.c
#     output = take!(REPLINSTANCE.c)
#     # @show "TOGREPL.take!", output
#     output
# end
# put!(::Type{repl}, a) = begin
#     # @show "TOGREPL.put!", a
#     println(stdout, a)
# end
# state(::Type{repl}) = "TOGREPL.REPL"

# old code
# struct REPLInput <: Peripheral
#     c::Channel{String}
# end
# take!(::REPLInput) = take!(REPL.c)
# const REPLINPUT = REPLInput(Channel{String}(Inf))
# struct REPLOutput <: Peripheral end
# put!(::REPLOutput, a) = println(stdout, a)
# old code

# repl_parse(s) = begin
#     # @show "TOGREPL.repl_parse", s
#     put!(REPLINSTANCE.c, string(strip("""$s""")))
# end
# repl_parse(s) = println(s)
const REMOTE_REPL_TASK = Ref{Task}()
const NAME = Ref{String}()
# function awaken(GOD)
function connect(; start_key, port, host=Sockets.localhost)
    @show "connect", start_key, port, host
    connection = RemoteREPL.connect_remote(host, port)
    @show "connect", connection
    name = @remote(connection, TOGΩ.TOGREPL.NAME[])
    @show "connect", name
    initrepl(
        s -> RemoteREPL.remotecmd(connection, stdout, s),
        # start_key=(start_key),
        mode_name=Symbol(name),
        prompt_text="$name> ",
        # startup_text=false,
    )
end
function awaken(; name, remotereplport=TOGPort.openport(), path=".")
    NAME[] = name
    # atreplinit(r -> begin
    #     ReplMaker.initrepl(
    #         repl_parse,
    #         repl=r,
    #         prompt_text="> ",
    #         prompt_color=:light_cyan,
    #         start_key="\\C-G", # "\x07",
    #         mode_name="GOD",
    #     )
    #     listen(REPLINSTANCE)
    # end
    # )
    # TOGAwaken.writeremotereplport(port=remotereplport)
    TOGAwaken.writeremotereplport(path=path, port=remotereplport)
    # REMOTE_REPL_TASK[] = serve_repl(remotereplport)
    if isinteractive()
        REMOTE_REPL_TASK[] = @async serve_repl(remotereplport)
    else
        serve_repl(remotereplport)
    end

    # remotereplport
    # write(stdin.buffer, "\x07")
    # GOD && write(stdin.buffer, "\x07")
end
# __init__() = atexit(sleep)
function sleep(;path=".")
    filename = TOGAwaken.remotereplportfile(path=path)
    isfile(filename) && rm(filename)
end

# serve_repl([address=Sockets.localhost,] port=27754; [on_client_connect=nothing])
# initrepl(parser::Function;
#                     prompt_text = "myrepl> ",
#                     prompt_color = :blue,
#                     start_key = ')',
#                     repl = Base.active_repl,
#                     mode_name = :mylang,
#                     show_function = nothing,
#                     show_function_io = stdout,
#                     valid_input_checker::Function = (s -> true),
#                     keymap::Dict = REPL.LineEdit.default_keymap_dict,
#                     completion_provider = REPL.REPLCompletionProvider(),
#                     sticky = true,
#                     startup_text=true
#                     )


end
