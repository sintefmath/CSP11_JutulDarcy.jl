# CSP11_JutulDarcy.jl

Simulate the [CSP11](https://www.spe.org/en/csp/) CO2 storage benchmark using [JutulDarcy.jl](https://github.com/sintefmath/JutulDarcy.jl) on meshes made in [MRST](https://mrst.no/). Currently contains scripts for case B and C only.

## First time setup

To use, check out this repository, do `git lfs pull` and open a Julia 1.9+ session in the local environment:

```bash
julia --project=/path/to/CSP11_JutulDarcy.jl
```

Add the local module with various helper functions:

```julia
]           # Enter package mode
dev ./CSP11 # Add local package
instantiate # Add dependencies and precompile
```

## Running cases

You can now run the example:

```julia
include("run_mrst_grid_spe11.jl")
```
