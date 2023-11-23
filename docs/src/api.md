# API

Full reference.

## Coordinates

```@autodocs
Modules = [SargassumBOMB]
Pages = [
    "coordinates.jl"
]
```

## Interpolants

<!-- preprocess-all.jl is inclded by itp-construct.jl, not sure why it needs to be included twice -->
```@autodocs
Modules = [SargassumBOMB]
Pages = [
    "interpolants/itp-core.jl",
    "interpolants/itp-construct.jl",
    "data/preprocessed/preprocess-all.jl" 
]
```

## Biology

```@autodocs
Modules = [SargassumBOMB]
Pages = [
    "biology.jl"
]
```

## Geography

```@autodocs
Modules = [SargassumBOMB]
Pages = [
    "geography.jl"
]
```

## Raft Parameters

```@autodocs
Modules = [SargassumBOMB]
Pages = [
    "raft-parameters.jl"
]
```

## Physics

```@autodocs
Modules = [SargassumBOMB]
Pages = [
    "physics.jl"
]
```

## Control

```@autodocs
Modules = [SargassumBOMB]
Pages = [
    "control.jl"
]
```

## Trajectories

```@autodocs
Modules = [SargassumBOMB]
Pages = [
    "trajectories.jl"
]
```

## Plotting

```@autodocs
Modules = [SargassumBOMB]
Pages = [
    "plotting/plotting-core.jl",
    "plotting/plotting-itp.jl"
]
```

## Main

```@autodocs
Modules = [SargassumBOMB]
Pages = [
    "main.jl"
]
```