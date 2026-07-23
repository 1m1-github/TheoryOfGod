# Theory Of God

## Install julia

```
curl -fsSL https://install.julialang.org | sh -s -- --yes
```

## Create directory

```
mkdir TheoryOfGod && cd TheoryOfGod
```

## Run TheoryOfGod

```
julia --optimize=3 --threads=auto --quiet --load tog.jl --eval 'tog.awaken()'
```

# Explain

## LoopOS

## god

## GOD

## Learning

## [0,1]^(.)

## Consciousness

## Autonomy

## Freedom

# Debug

## Ω

```
cd Ω
```

### Awaken

```
julia --optimize=3 --threads=auto --interactive --quiet --project=.tog --eval 'using TOGOmega;TOGOmega.awaken()'
```

### Update

```
julia --optimize=3 --threads=auto --quiet --project=.tog --eval 'using Pkg;Pkg.update()'
```

## Dona

```
cd Dona
```

### Awaken

```
env JULIA_PROJECT=$(pwd)/.tog JULIA_DEPOT_PATH=$(pwd)/.tog/julia:$HOME/.julia julia --optimize=3 --threads=auto --interactive --quiet --project=.tog --eval 'include("/Users/1m1/TheoryOfGod/Dona/.tog/Dona.jl");using .Dona;Dona.awaken(path="/Users/1m1/TheoryOfGod/Dona",universe="../Ω")'
```

## Update

```
env JULIA_PROJECT=$(pwd)/.tog JULIA_DEPOT_PATH=$(pwd)/.tog/julia:$HOME/.julia julia --optimize=3 --threads=auto --quiet --project=.tog --eval 'using Pkg;Pkg.update()'
```

# TODO
### awaken with remoterepl, repl/process io fwding, sandbox
### peripherals/explanations/cleanup
### env vars for api keys
### registry: rm, clean
### online/offline learning github
### first newpkg or update with new pkg fails then passes
### learning rmpkgs, pkgs that do not exist errs, loop updatepkg whilst reducing num erred pkgs
### move octahedron, random
### TOGLearn
### Share octahedron server
### information explained
### energy explained
### power explained
### autonomy explained
### freedom explained
### alignment explained
### consciousness explained
### intelligence (complexity) explained
### metabolic balance explained
### existence maxxed via growth explained
### same ports being used
### Ω add Pkg, first to zmq is root?
### IO long memory LoopOS
### tog privacy: create/observe, zmq
### explain "fair" in license
