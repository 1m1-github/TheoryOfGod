module TOGAwaken

export awakengod

const TOGDIR = ".tog"
const REGISTRYPATH = joinpath(DEPOT_PATH[1], "registries", "TOGRegistry")
const JULIA_ARGS = ["--optimize=3", "--threads=auto", "--quiet"]
const JULIA_ARGS_INTERACTIVE = [JULIA_ARGS..., "--interactive"]
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

function julia(; code, path=".", project=joinpath(path, TOGDIR), dev=joinpath(project, "julia", "dev"), depot=joinpath(project, "julia"), args=JULIA_ARGS, wait=true, sandbox=false)
    @info "TOGAwaken.jl, julia"
    cd(path) do
        depot *= ":" * join(DEPOT_PATH, ":")
        cmd = addenv(
            `julia $args -e $code`,
            "JULIA_PKG_DEVDIR" => dev,
            "JULIA_PROJECT" => project,
            "JULIA_DEPOT_PATH" => depot,
        )
        sandbox && (cmd = sandboxcmd(cmd=cmd))
        run(cmd, wait=wait)
    end
end

function installΩ(;path)
    @info "TOGAwaken.jl, installΩ", path
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
    @info "TOGAwaken.jl, updateΩ", path
    julia(
        path=path,
        project=TOGAwaken.TOGDIR,
        dev=joinpath(ENV["HOME"], ".julia", "dev"),
        depot=DEPOT_PATH[1],
        code="""using Pkg;Pkg.update()""")
end
function awakenΩ()
    @info "TOGAwaken.jl, awakenΩ"
    path=joinpath(pwd(), "Ω")
    installΩ(path=path)
    updateΩ(path=path)
    julia(
        path=path,
        project=TOGAwaken.TOGDIR,
        dev=joinpath(ENV["HOME"], ".julia", "dev"),
        depot=DEPOT_PATH[1],
        # args=TOGAwaken.JULIA_ARGS_INTERACTIVE,
        code="""using TOGOmega;TOGOmega.awaken()""", wait=false)
end

function installgod(; path, pkgs)
    @info "TOGAwaken.jl, installgod"
    project=TOGAwaken.TOGDIR
    isdir(joinpath(path, project)) && return
    isdir(path) || mkpath(path)
    pkgadd = isempty(pkgs) ? "" : "Pkg.add([" * join(map(pkg->""""$pkg\"""", pkgs), ',') * "])"
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
        $pkgadd
        """)
end
function updategod(; path)
    @info "TOGAwaken.jl, updategod"
    julia(
        path=path,
        project=TOGAwaken.TOGDIR,
        code="""using Pkg;Pkg.update()""")
end
"""
Awakens a god 
example: `awakengod(path="Anna", pkgs=["Dates"], universe=TOGgod.ARGS[][:universe])`
"""
function awakengod(; args...)
    @info "TOGAwaken.jl, awakengod"
    path=args[:path]
    pkgs=get(args, :pkgs, String[])
    installgod(path=path, pkgs=pkgs)
    updategod(path=path)
    argsparts = ["$k = $(repr(v))" for (k, v) in pairs(args)]
    argsstring = join(argsparts, ",")
    name = basename(path)
    godfile = joinpath(path, TOGDIR, "$name.jl")
    @info "TOGAwaken.awakengod", argsstring
    julia(
        path=path,
        project=TOGAwaken.TOGDIR,
        # args=TOGAwaken.JULIA_ARGS_INTERACTIVE,
        code="""include("$godfile");using .$name;$name.awaken($argsstring)""", wait=false)
    # code="""using $name;$name.awaken(universe="$universe")""") # DEBUG
end

function sandboxcmd(; path=TOGDIR, cmd)
    args = if Sys.isapple()
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
    prependargs(cmd, args)
end

function prependargs(cmd::Cmd, args...)
    strs = string.(args)
    exec = cmd.exec
    new_exec = vcat(exec[1:1], strs, exec[2:end])
    Cmd(new_exec)
end

end
