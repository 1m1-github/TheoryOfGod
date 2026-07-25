module TOGLearning

# todo handle fails

export newpkg, updatepkg

using Pkg, TOML, LocalRegistry, Git
using Pkg.Types: PackageSpec, Context

const REGISTRYNAME = "TOGRegistry"
const JULIACODEPATH = joinpath(DEPOT_PATH[1], "dev")

const LICENSEFILE = "LICENSE"
const LICENSE = """
Study it, use it, enjoy it.
Any one deriving value from this should share a fair value >= 0.
Fair value is derived via negotiation xor ultimately a court.
"""
const READMEFILE = "README.md"
const README(name) = "# $name"
const GITIGNOREFILE = ".gitignore"
const GITIGNORE = """
.*
**/.*
!**/.gitignore
Manifest.toml
tmp*
"""

pkgdir(; name) = joinpath(JULIACODEPATH, name)
projecttoml(; name) = joinpath(pkgdir(name=name), "Project.toml")

"""
Each of your abilities lives modularly in a `Pkg`.
`newpkg` creates a new `Pkg` inside your Julia dev dir.
The process of learning a new ability is to write a new module source file, then create a `newpkg`, then `learn` to have the ability present.
`newpkg` will `generate` a new `Pkg`, `cp`/`mv` files, `add` `Pkg`s, `git commit`, `register` to the `LocalRegistry` and `rm` all on error.
# Arguments
- `name::String`: The name of the `Pkg`.
- `file=String[]`: The `file`s to be copied xor moved to the `src` folder of the `Pkg`.
- `pkg=String[]`: The `Pkg`s to be `add`ed to this `Pkg`.
- `mvfile::Bool=false`: Whether to copy xor move the `file`s.
"""
function newpkg(; name::String, file=String[], pkg=String[], mvfile=false)
    path = pkgdir(name=name)
    isdir(path) && return
    # @info "TOGLearning.jl, newpkg"
    Pkg.generate(path)
    try
        changefiles(name=name, addfile=file, rmfile=String[], cpmv=mvfile ? mv : cp, init=true)
        changepkgs(name=name, addpkg=pkg, rmpkg=String[])
        version = initversion(name=name)
        addcommit(path=path, commitmessage=version)
        # newremoterepo(path=path)
        registerpkg(name=name)
    catch e
        rm(path, force=true, recursive=true)
        rethrow(e)
    end
end

"""
Update an existing `Pkg` in your Julia dev dir by force copying `file`s to the `src` dir, removing `rmfile`s from the `src` dir, `add`ing `addpkg`s, removing `rmpkg`s.
`updatepkg` will `cp`/`mv` files, `add`/`rm` `Pkg`s, increase the `version`, `git commit`, `register` to the `LocalRegistry`.
# Arguments
- `name::String`: The name of the `Pkg`.
- `addfile=String[]`: The `file`s to be copied xor moved to the `src` folder of the `Pkg` by force.
- `rmfile=String[]`: The `file`s to be removed from the `src` folder of the `Pkg`.
- `addpkg=String[]`: The `Pkg`s to be `add`ed to this `Pkg`.
- `rmpkg=String[]`: The `Pkg`s to be `add`ed to this `Pkg`.
- `mvfile::Bool=false`: Whether to copy xor move the `file`s.
"""
function updatepkg(; name::String, addfile=String[], addpkg=String[], rmfile=String[], rmpkg=String[], mvfile=false)
    # @info "TOGLearning.jl, updatepkg"
    changefiles(name=name, addfile=addfile, rmfile=rmfile, cpmv=mvfile ? mv : cp)
    changepkgs(name=name, addpkg=addpkg, rmpkg=rmpkg)
    path = pkgdir(name=name)
    isdirty(path=path) || return
    version = updateversion(name=name)
    addcommit(path=path, commitmessage=version)
    # if hasremote(path=path)
    #     pushremote(path=path)
    # else
    #     newremoterepo(path=path)
    # end
    registerpkg(name=name)
end
# function rmpkg(; name::String, pushregistry=false, githubuser=get(ENV, "GITHUB_USER", ""), githubauth=get(ENV, "GITHUB_AUTH", ""))
# rmdir(joinpath(JULIACODEPATH, name))
# path = registrytoml()
# registry = TOML.parsefile(path)
# pkgkeys = filter(k -> registry["packages"][k]["name"] == name, keys(registry["packages"]))
# if !isempty(pkgkeys)
# pkgkey = only(pkgkeys)
# rmdir(joinpath(LOOPOSREGISTRYPATH, registry["packages"][pkgkey]["path"]))
# delete!(registry["packages"], pkgkey)
# open(path, "w") do io
# TOML.print(io, registry)
# end
# addcommitpush(LOOPOSREGISTRYPATH, push=pushregistry)
# end
# rmrepo(name, githubuser, githubauth)
# end

# """
# Convenience method to copy an existing pkg 
# """
# function cppkg(; name::String, newname::String)
#     @info "TOGLearning.jl, cppkg"
#     files = readdir(joinpath(pkgdir(name=name), "src"), join=true)
#     project = TOML.parsefile(projecttoml(name=name))
#     pkgs = haskey(project, "deps") ? collect(keys(project["deps"])) : String[]
#     newpkg(name=newname, files=files, pkgs=pkgs)
# end
# function mvpkg(; name::String, newname::String, pushregistry=false, githubuser=get(ENV, "GITHUB_USER", ""), githubauth=get(ENV, "GITHUB_AUTH", ""))
#     cppkg(name=name, newname=newname, pushregistry=pushregistry, githubuser=githubuser, githubauth=githubauth)
#     rmpkg(name=name)
# end
# function checkpkg(;name)
#     # cd(pkgdir(name=name)) do
#         oldenv = Base.active_project()
#         # Pkg.activate(".")
#         Pkg.activate(pkgdir(name=name))
#         Pkg.precompile()
#         @show "reactivate", oldenv
#         Pkg.activate(oldenv)
#         @show "reactivated oldenv"
#     # end
# end

# rmdir(path) = isdir(path) && rm(path, recursive=true)
function addfile(; name, file, content)
    file = joinpath(pkgdir(name=name), file)
    !isfile(file) && write(file, content)
end
srcfile(; name, file) = joinpath(pkgdir(name=name), "src", basename(file))
function changefiles(; name, addfile=String[], rmfile=String[], cpmv=cp, init=false)
    if init
        addfile(name=name, file=LICENSEFILE, content=LICENSE)
        addfile(name=name, file=GITIGNOREFILE, content=GITIGNORE)
        addfile(name=name, file=READMEFILE, content=README(name))
    end
    for file = addfile
        cpmv(file, srcfile(name=name, file=file), force=true)
    end
    for file = rmfile
        rm(srcfile(name=name, file=file))
    end
end

function changepkg(pkg, f)
    if startswith(pkg, "http")
        f(url=pkg)
    elseif ispath(pkg)
        f(path=pkg)
    else
        f(pkg)
    end
end
function changepkgs(; name, addpkg=String[], rmpkg=String[])
    # cd(pkgdir(name=name)) do
    oldenv = Base.active_project()
    # Pkg.activate(".")
    Pkg.activate(pkgdir(name=name))
    # Pkg.Registry.add("General")
    # Pkg.Registry.add(path="/Users/1m1/.julia/registries/TOGRegistry")
    isempty(rmpkg) || Pkg.rm(rmpkg)
    isempty(addpkg) || Pkg.add(addpkg)
    Pkg.update(update_registry=false) # DEBUG
    # Pkg.update()
    # isempty(pkgs) || Pkg.develop(pkgs)
    Pkg.precompile()
    Pkg.activate(oldenv)
    # end
end

function registerpkg(; name)
    register(
        pkgdir(name=name),
        registry=REGISTRYNAME,
        repo=remoteurllocal(name=name),
        push=false
    )
end

function resolve(; name)
    cd(pkgdir(name=name)) do
        oldenv = Base.active_project()
        Pkg.activate(".")
        Pkg.resolve()
        Pkg.activate(oldenv)
    end
end
function changeversion(; name, newversion)
    path = projecttoml(name=name)
    project = TOML.parsefile(path)
    version = VersionNumber(project["version"])
    project["version"] = string(newversion(version))
    delete!(project, "compat")
    open(path, "w") do file
        TOML.print(file, project)
    end
    resolve(name=name)
    project["version"]
end
initversion(; name) = changeversion(name=name, newversion=_ -> v"1")
updateversion(; name) = changeversion(name=name, newversion=v -> VersionNumber(v.major + 1))

function git1m1()
    g = Sys.which("git")
    g === nothing && return nothing
    chomp(read(`$g config user.name`, String)) == "1m1" ? g : nothing
end
function gitcmd()
    g = git1m1()
    isnothing(g) && return git()
    g
end

isdirty(; path=".") =
    cd(path) do
        !isempty(read(`$(gitcmd()) status --porcelain`))
    end
remoteurllocal(; name) = joinpath(JULIACODEPATH, name)
hasremote(; path=".") =
    cd(path) do
        !isempty(readlines(`$(gitcmd()) remote`))
    end
getremoteurl(; path=".") =
    cd(path) do
        readlines(`$(gitcmd()) remote get-url origin`) |> only
    end
addsetremote(; path, addset) =
    cd(path) do
        run(`$(gitcmd()) remote $addset origin $(remoteurllocal(name=basename(path)))`)
    end
addremote(; path) = addsetremote(path=path, addset="add")
setremote(; path) = addsetremote(path=path, addset="set-url")
pushremote(; path=".") =
    cd(path) do
        g = git1m1()
        !isnothing(g) && run(`$g push -f -u origin main`)
    end
function addcommit(; path, commitmessage=".")
    cd(path) do
        isnew = !isdir(".git")
        isnew && run(`$(gitcmd()) init`)
        run(`$(gitcmd()) add .`)
        !isnew && !isdirty(path=".") && return
        run(`$(gitcmd()) commit -m $commitmessage`)
    end
end

end
