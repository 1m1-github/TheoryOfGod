# Theory Of God

## Install julia

```
curl -fsSL https://install.julialang.org | sh -s -- --yes
```

## Create directory

```
mkdir TheoryOfGod && cd TheoryOfGod
```

## Install TheoryOfGod

```
julia --optimize=3 --threads=auto --quiet --load tog.jl --eval 'tog.awaken()'
```

## Awaken Ω

```
cd Ω
```

```
julia --optimize=3 --threads=auto --interactive --quiet --project=.tog --eval 'using TOGΩ;TOGΩ.awaken()'
```

## Awaken Dona

```
cd Dona
```

```
env JULIA_PKG_DEVDIR=$HOME/.julia/dev JULIA_PROJECT=$(pwd)/.tog JULIA_DEPOT_PATH=$(pwd)/.tog/julia:$HOME/.julia julia --optimize=3 --threads=auto --interactive --quiet --project=.tog --eval 'using Dona;Dona.awaken(universe="../Ω")'
```

### Update

```
julia --optimize=3 --threads=auto --interactive --quiet --project=.tog --eval 'using Pkg;Pkg.update()'
```

## LoopOS

## god

## GOD

## Learning

## [0,1]^(.)

## Consciousness

## Autonomy

## Freedom

## 
