# CSP11_JutulDarcy.jl

## First time setup

To use, check out this repository, do `git lfs pull` and open a Julia session in the local environment:

```bash
julia --project=/path/to/CSP11_JutulDarcy.jl
```

Add the local module with various helper functions:

```julia
]           # Enter package mode
dev ./CSP11 # Add local package
instantiate # Add dependencies
```

## Running cases

You can now run the example:

```julia
include("run_mrst_grid_spe11.jl")
```
