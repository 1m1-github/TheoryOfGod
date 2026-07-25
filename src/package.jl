# full process
using Pkg
# Pkg.offline(true)
Pkg.update()
Pkg.add(["Pkg","Revise", "Git", "LocalRegistry"])
using Revise, Git, LocalRegistry, TOML
# Pkg.Registry.add()
# Pkg.Registry.add("General")
# Pkg.Registry.add(url="https://github.com/1m1-github/TOGRegistry.git")
# Pkg.Registry.add(path=joinpath(DEPOT_PATH[1],"registries","TOGRegistry"))
# Pkg.add(["TOGLearning"])
# Pkg.develop("TOGLearning")
using TOGLearning
# githubuser=get(ENV, "GITHUB_USER", "")
# githubauth=get(ENV, "GITHUB_AUTH", "")
# create_repo(
#     GitHub.owner(githubuser),
#     "TOGRegistry",
#     auth=authenticate(githubauth),
#     )
# create_registry("TOGRegistry", "https://github.com/1m1-github/TOGRegistry.git", push=true)

# allfile(; name) = joinpath(ENV["HOME"], "Documents", "TOGAll", name * ".jl")
# allfile(; name) = joinpath(ENV["HOME"], "dev", name, "src", name * ".jl")
allfile(; name) = joinpath(ENV["HOME"], "TheoryOfGod/dev", name, "src", name * ".jl")
# allfile(; name) = joinpath(DEPOT_PATH[1], "dev", name, "src", name * ".jl")
function newpkg(; name, pkgs=String[])
    setglobal!(Base.MainInclude, :err, nothing)
    TOGLearning.newpkg(
        name=name,
        files=[allfile(name=name)],
        pkgs=pkgs,
        mvfiles=false,
        githubuser="",
        )
    # path = TOGLearning.pkgdir(name=name)
    # version = TOGLearning.initversion(name=name)
    # TOGLearning.addcommit(path=path, commitmessage=version)
    # TOGLearning.registerpkg(name=name)
    println((isnothing(err) ? "OK" : string(err))* " " * name)
end



# registry = TOML.parsefile(joinpath(DEPOT_PATH[1], "registries", "TOGRegistry", "Registry.toml"))
# name="TOGCreateServer"
# for (_,pkg) = registry["packages"]
#     name = pkg["name"]
#     @show name
#     path=joinpath(DEPOT_PATH[1],"dev",name)
#     TOGLearning.hasremote(path=path) || TOGLearning.addremote(path=path)
#     exists, _ = TOGLearning.remoterepoexists(name=name)
#     exists || TOGLearning.newremoterepo(path=path)
#     TOGLearning.pushremote(path=path)
# end

# develop all registered pkgs
using Pkg, TOML
registry = TOML.parsefile(joinpath(DEPOT_PATH[1], "registries", "TOGRegistry", "Registry.toml"))
for (_,pkg) = registry["packages"]
    @show pkg["name"]
    Pkg.develop(pkg["name"])
end


TOGLearning.updatepkg(name="TOGOctahedron", addpkg=["TOGColor"])
files(name) = [joinpath(ENV["HOME"], "TheoryOfGod", "src", name * ".jl")]
TOGLearning.newpkg(name="TOGIntelligenceJanet", pkgs=["LoopOS", "TOGState", "TOGLogging"], files=files("TOGIntelligenceJanet"))

cd(joinpath(DEPOT_PATH[1], "dev"))

# cp from src to dev
for filename = readdir(joinpath(ENV["HOME"], "TheoryOfGod", "src"), join=true)
    endswith(filename, ".jl") || continue
    name = join(split(basename(filename),'.')[1:end-1],'.')
    target = joinpath(DEPOT_PATH[1],"dev", name, "src", name * ".jl")
    isfile(target) || continue
    hash(read(filename, String)) == hash(read(target, String)) && continue
    @info filename, target
    cp(filename, target, force=true)
end

# update all registered pkgs
for pkg = readdir(joinpath(DEPOT_PATH[1], "dev"), join=true)
    isdir(pkg) || continue
    TOGLearning.isdirty(path=pkg) || continue
    @show pkg
    TOGLearning.updatepkg(name=basename(pkg))
end
Pkg.update()

for pkg = readdir(joinpath(DEPOT_PATH[1], "dev"), join=true)
    # isdir(pkg) || continue
    isfile(joinpath(pkg, "Manifest.toml")) || continue
    @show pkg
    rm(joinpath(pkg, "Manifest.toml"))
    # Pkg.activate(pkg)
    # Pkg.resolve()
end

# cp readme
for pkg = readdir(joinpath(ENV["HOME"],"Downloads/dev 2"),join=true)
    filename = joinpath(pkg, "README.md")
    isfile(filename) || continue
    @show pkg
    cp(filename, joinpath(DEPOT_PATH[1], "dev", basename(pkg), "README.md"), force=true)
end

TOGLearning.getremoteurl(path="LoopOS")

Pkg.Registry.add(url="https://github.com/1m1-github/TOGRegistry.git")
TOGLearning.rmremoterepo(name="TOGgod")

registry = TOML.parsefile(joinpath(DEPOT_PATH[1], "registries", "TOGRegistry", "Registry.toml"))
for (_,pkg) = registry["packages"]
    @show pkg["name"]
    register(
                TOGLearning.pkgdir(name=pkg["name"]),
                repo=TOGLearning.remoteurl(name=pkg["name"]),
                registry=TOGLearning.REGISTRYNAME,
                push=true
            )
end

for (_,pkg) = registry["packages"]
    @show pkg["name"]
    # TOGLearning.updateversion(name=pkg["name"])
    path = TOGLearning.pkgdir(name=pkg["name"])
    TOGLearning.addcommit(path=path)
    if TOGLearning.hasremote(path=path)
        TOGLearning.pushremote(path=path)
    else
        TOGLearning.newremoterepo(path=path, githubuser=githubuser, githubauth=githubauth)
    end
end

# path = readdir(joinpath(DEPOT_PATH[1], "dev"), join=true)[2]
for path = readdir(joinpath(DEPOT_PATH[1], "dev"), join=true)
    isdir(path) || continue
    @show path
    TOGLearning.setremote(path=path)
end

Pkg.activate(joinpath("dev", "TOGΩ"))
Pkg.update()
Pkg.precompile()
for pkg = readdir(joinpath(DEPOT_PATH[1], "dev"), join=true)
    isdir(pkg) || continue
    @show pkg
    Pkg.activate(pkg)
    Pkg.update()
    Pkg.precompile()
end

pwd()
for pkg = readdir()
    startswith(pkg, "TOG") || continue
    pkg == "TOGAll" && continue
    mv(joinpath(pkg, "src", pkg * ".jl"), joinpath("TOGAll", pkg * ".jl"))
end

Pkg.Apps.add(path=joinpath(DEPOT_PATH[1], "dev", "TOGΩ"))
Pkg.Apps.update("tog")
Pkg.Apps.update(name)
Pkg.Apps.rm("tog")
Pkg.update()
Pkg.activate(name)
Pkg.precompile()
Pkg.activate()
setglobal!(Base.MainInclude, :err, nothing)

Pkg.Apps.rm("tog")
Pkg.Apps.add(path=joinpath(DEPOT_PATH[1], "dev", "TOGΩ"))

# ====

name="TOGΩ"
pkgs=["Git", "GitHub"]
pkgs=String[]
TOGLearning.newpkg(
        name=name,
        files=[name * ".jl"],
        pkgs=pkgs,
        mvfiles=false,
    )
Base.active_project()
Pkg.status()
path=joinpath(DEPOT_PATH[1], "dev", name)
TOGGit.addcommitpush(path=path)
TOGGit.setremote(path=path)
TOGGit.newremoterepo(path=path)
TOGGit.push(path=path)
LocalRegistry.register(
        path,
        registry=TOGLocalRegistry.REGISTRYNAME,
        push=true,
    )

project = TOML.parsefile(joinpath(path, "Project.toml"))
version = VersionNumber(project["version"])
project["version"] = string(newversion(version))
project["apps"]
delete!(project, "compat")
open(path, "w") do file
    TOML.print(file, project)
end
Pkg.Apps.add(path=path)
name="TOGLocalRegistry"
name="TOGLearning"


for pkg = readdir(joinpath(DEPOT_PATH[1], "dev"), join=true)
    isdir(pkg) || continue
    name=basename(pkg)
    devhash=hash(read(joinpath(pkg, "src", name*".jl")))
    togallhash=hash(read(joinpath(ENV["HOME"], "Documents", "TOGAll", name*".jl")))
    devhash == togallhash || @show name
end


for file = readdir(joinpath(ENV["HOME"], "Documents", "LoopOSAll"))
    endswith(file, ".jl") || continue
    startswith(file, "LoopOS") || continue
    file == "LoopOS.jl" && continue§
    name=file[1:end-3]
    @show name
    TOGLearning.rmremoterepo(name=name)
end
TOGLearning.rmremoterepo(name="LoopOSRegistry")

using Pkg, TOML
registry = TOML.parsefile(joinpath(DEPOT_PATH[1], "registries", "TOGRegistry", "Registry.toml"))
name="TOGAwaken"
using GitHub
for (_,pkg) = registry["packages"]
    name = pkg["name"]
    startswith(name, "TOG") || continue
    @show name
    try delete_repo(
            repo("$githubuser/$name"),
            auth=authenticate(githubauth),
        )
    catch e @show e end
end
TOGLearning.rmremoterepo(name="LoopOS")

remoterepoexists(; name, githubuser=get(ENV, "GITHUB_USER", "")) =
    try
        true, repo("$githubuser/$name")
    catch
        false, nothing
    end
function rmremoterepo(; name, githubuser=get(ENV, "GITHUB_USER", ""), githubauth=get(ENV, "GITHUB_AUTH", ""))
    if !isempty(githubuser) && !isempty(githubauth)
        exists, repo = remoterepoexists(name=name, githubuser=githubuser)
        exists && delete_repo(
            repo,
            auth=authenticate(githubauth),
        )
    end
end