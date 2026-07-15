module TOGAwaken

const TOGDIR = ".tog"
const REGISTRYPATH = joinpath(DEPOT_PATH[1], "registries", "TOGRegistry")
const JULIA_ARGS = ["--optimize=3", "--threads=auto", "--quiet"]
const JULIA_ARGS_INTERACTIVE = [JULIA_ARGS..., "--interactive"]

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
    isawake(path=path) && error("TOGgod at $path is already awake.")
    writepid(path=path)
end

function julia(; path=".", project=joinpath(path, TOGDIR), depot=joinpath(project, "julia"),args=JULIA_ARGS, wait=true, code)
    cd(path) do
        dev = joinpath(ENV["HOME"], ".julia", "dev")
        depot *= ":" * join(DEPOT_PATH, ":")
        cmd = addenv(
                `julia $args -e $code`,
                "JULIA_PKG_DEVDIR" => dev,
                "JULIA_PROJECT" => project,
                "JULIA_DEPOT_PATH" => depot,
            )
        # run(sandbox(cmd), wait=wait)
        run(cmd, wait=wait)
    end
end

function installΩ()
    dir = "Ω"
    project=TOGAwaken.TOGDIR
    isdir(joinpath(dir, project)) && return
    isdir(dir) || mkdir(dir)
    TOGAwaken.julia(
    path=dir,
    project=project,
    depot=DEPOT_PATH[1],
    code="""
    using Pkg
    Pkg.add("TOGΩ")
    """)
end
function updateΩ()
    TOGAwaken.julia(
    path="Ω",
    project=TOGAwaken.TOGDIR,
    depot=DEPOT_PATH[1],
    args=TOGAwaken.JULIA_ARGS_INTERACTIVE,
    code="""using Pkg;Pkg.update()""", wait=false)
end
function awakenΩ()
    installΩ()
    updateΩ()
    TOGAwaken.julia(
    path="Ω",
    project=TOGAwaken.TOGDIR,
    depot=DEPOT_PATH[1],
    args=TOGAwaken.JULIA_ARGS_INTERACTIVE,
    code="""using TOGΩ;TOGΩ.awaken()""", wait=false)
end

function installgod(;name)
    dir = name
    project=TOGAwaken.TOGDIR
    isdir(joinpath(dir, project)) && return
    isdir(dir) || mkdir(dir)
    TOGAwaken.julia(
    path=dir,
    project=project,
    code="""
    using Pkg
    Pkg.Registry.add()
    Pkg.Registry.add(path="$REGISTRYPATH")
    Pkg.add("$name")
    """)
end
function updategod(;name)
    TOGAwaken.julia(
    path=name,
    project=TOGAwaken.TOGDIR,
    code="""using Pkg;Pkg.update()""", wait=false)
end
function awakengod(;name)
    installgod(name=name)
    updategod(name=name)
    TOGAwaken.julia(
    path=name,
    project=TOGAwaken.TOGDIR,
    args=TOGAwaken.JULIA_ARGS_INTERACTIVE,
    code="""using $name;$name.awaken(universe="../Ω")""", wait=false)
end

function sandbox(;path, cmd)
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
