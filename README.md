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

## Running in docker

To run it inside a [Docker](https://docs.docker.com/desktop/) image, either build the image described uner _docker/jutul-spe11csp-juliaimg.Dockerfile_ or pull it from public [docker repository](https://hub.docker.com)

```bash
docker build --build-arg VERSION=v0 -t jutul-juliaimg-spe11csp -f jutul-spe11csp-juliaimg.Dockerfile .
```

or

```bash
docker push jafranc/jutul-juliaimg-spe11csp:latest
```

Then, as docker environement cannot easily be used with Cairo-GUI, _Project.toml_ has to be adjusted discarding inaccessible deps. It is done under the folder _docker/Project.toml_.

Eventually you can run it

```bash
docker run -it --rm -v /path/to/output:/tmp image_name julia -e "import Pkg;Pkg.add([\"Jutul\",\"JutulDarcy\",\"HYPRE\"]);Pkg.develop(path=\"/opt/spe11csp/CSP11/docker\");include(\"docker_run_mrst_grid_spe11.jl\")"
```

where _image_name_ is to be replaced either by name taken during local build step (suggested above _jutul-juliaimg-spe11csp_) or by pulled image name _jafranc/jutul-juliaimg-spe11csp_. _/path/to/output_ to be replaced by the path to a local folder. Note that generated results will be root-owned as per vanilla docker logic. 
