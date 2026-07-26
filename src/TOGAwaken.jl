module TOGAwaken

export awakengod

const TOGDIR = ".tog"
const REGISTRYPATH = joinpath(DEPOT_PATH[1], "registries", "TOGRegistry")
const JULIA_ARG = ["--optimize=3", "--threads=auto", "--quiet"]
const JULIA_ARG_INTERACTIVE = [JULIA_ARG..., "--interactive"]
const ALREADYRUNNINGEXITCODE = 8888

globalpath(; path=".") = joinpath(pwd(), path)
router(; path=".") = "ipc://$(globalpath(path=path))/$TOGDIR/router.ipc" # change to tcp if on separate machines
pub(; path=".") = "ipc://$(globalpath(path=path))/$TOGDIR/pub.ipc" # change to tcp if on separate machines
togobserve(; path=".") = "ipc://$(globalpath(path=path))/$TOGDIR/togobserve.ipc" # change to tcp if on separate machines
togcreate(; path=".") = "ipc://$(globalpath(path=path))/$TOGDIR/togcreate.ipc" # change to tcp if on separate machines
pidfile(; path=".") = joinpath(path, TOGDIR, "pid")
readpid(; file=pidfile()) = parse(Int, read(file, String))
readremotereplport(; path=".", file=remotereplportfile(path=path)) = parse(Int, read(file, String))
writepid(; path=".") = write(pidfile(path=path), string(getpid()))
remotereplportfile(; path=".") = joinpath(path, TOGDIR, "remotereplport")
writeremotereplport(; port, path=".") = write(remotereplportfile(path=path), string(port))
broadcastbrowserportfile(; path=".") = joinpath(path, TOGDIR, "broadcastbrowserport")
writebroadcastbrowserport(; port, path=".") = write(broadcastbrowserportfile(path=path), string(port))
function rmpid(; path=".")
    file = pidfile(path=path)
    isfile(file) && rm(file)
end
# __init__() = atexit(rmpid)
sleep = rmpid
function isawake(; path=".")
    file = pidfile(path=path)
    isfile(file) || return false
    pid = readpid(file=file)
    success(`kill -0 $pid`)
end
function awaken(; path=".")
    if isawake(path=path)
        @error "Already awake at $path."
        exit(ALREADYRUNNINGEXITCODE)
    end
    writepid(path=path)
end

function julia(; code, path=".", project=joinpath(path, TOGDIR), dev=joinpath(project, "julia", "dev"), depot=joinpath(project, "julia"), arg=JULIA_ARG, wait=true, sandbox=false)
    # @info "TOGAwaken.jl, julia"
    cd(path) do
        depot *= ":" * join(DEPOT_PATH, ":")
        cmd = addenv(
            `julia $arg -e $code`,
            "JULIA_PKG_DEVDIR" => dev,
            "JULIA_PROJECT" => project,
            "JULIA_DEPOT_PATH" => depot,
        )
        sandbox && (cmd = sandboxcmd(cmd=cmd))
        run(cmd, wait=wait)
    end
end

function installΩ(;path)
    # @info "TOGAwaken.jl, installΩ", path
    project=TOGAwaken.TOGDIR
    isdir(joinpath(path, project)) && return
    isdir(path) || mkpath(path)
    julia(
        path=path,
        project=project,
        dev=joinpath(ENV["HOME"], ".julia", "dev"),
        depot=DEPOT_PATH[1],
        code="""using Pkg;Pkg.add("TOGOmega")""")
end
function updateΩ(;path)
    # @info "TOGAwaken.jl, updateΩ", path
    julia(
        path=path,
        project=TOGAwaken.TOGDIR,
        dev=joinpath(ENV["HOME"], ".julia", "dev"),
        depot=DEPOT_PATH[1],
        code="""using Pkg;Pkg.update()""")
end
function awakenΩ()
    # @info "TOGAwaken.jl, awakenΩ"
    path=joinpath(pwd(), "Ω")
    installΩ(path=path)
    updateΩ(path=path)
    julia(
        path=path,
        project=TOGAwaken.TOGDIR,
        dev=joinpath(ENV["HOME"], ".julia", "dev"),
        depot=DEPOT_PATH[1],
        code="""using TOGOmega;TOGOmega.awaken()""", wait=false)
end

function installgod(; path, pkg)
    # @info "TOGAwaken.jl, installgod"
    project=TOGAwaken.TOGDIR
    isdir(joinpath(path, project)) && return
    isdir(path) || mkpath(path)
    addpkg = isempty(pkg) ? "" : "Pkg.add([" * join(map(pkg->""""$pkg\"""", pkg), ',') * "])"
    currentdir = pwd()
    name = basename(path)
    julia(
        path=path,
        project=project,
        code="""
        using Pkg
        Pkg.Registry.add()
        Pkg.Registry.add(path="$REGISTRYPATH")
        cp(joinpath("$currentdir", "$name.jl"), joinpath(".tog", "$name.jl"))
        $addpkg
        """)
end
function updategod(; path)
    # @info "TOGAwaken.jl, updategod"
    julia(
        path=path,
        project=TOGAwaken.TOGDIR,
        code="""using Pkg;Pkg.update()""")
end
"""
Awakens a `god`. If new, creates a folder, connects the `god` to a `universe` and adds `Pkg`s to it.
Awaken `god`s to help you. These can be clones of yourself with differing abilities.
example to connect to the same universe: `awakengod(path="Anna", pkg=["Dates"], universe=TOGgod.ARG[][:universe])`
"""
function awakengod(; arg...)
    # @info "TOGAwaken.jl, awakengod"
    path=arg[:path]
    pkg=get(arg, :pkg, String[])
    installgod(path=path, pkg=pkg)
    updategod(path=path)
    argpart = ["$k = $(repr(v))" for (k, v) in pairs(arg)]
    argstring = join(argpart, ",")
    name = basename(path)
    godfile = joinpath(path, TOGDIR, "$name.jl")
    julia(
        path=path,
        project=TOGAwaken.TOGDIR,
        code="""include("$godfile");using .$name;$name.awaken($argstring)""", wait=false)
end

function sandboxcmd(; path=TOGDIR, cmd)
    arg = if Sys.isapple()
        profile = joinpath(path, ".sb")
        write(
            profile,
            """
            (version 1)
            (allow default)
            (deny file-write*)
            (allow file-write* (subpath "$path"))
            """
        )
        [
            "sandbox-exec",
            "-f",
            profile,
            "sh",
            "-c",
        ]
    elseif Sys.islinux()
        [
            "bwrap",
            "--ro-bind", "/", "/",
            "--bind", name, name,
            "--proc", "/proc",
            "--dev", "/dev",
            "--die-with-parent",
        ]
    else
        error("Currently only MacOS and Linux")
    end
    prependargs(cmd, arg)
end

function prependargs(cmd::Cmd, arg...)
    str = string.(arg)
    exec = cmd.exec
    new_exec = vcat(exec[1:1], str, exec[2:end])
    Cmd(new_exec)
end

end
