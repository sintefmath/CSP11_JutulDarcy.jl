# CSP11_JutulDarcy.jl

Simulate the [SPE11 CO₂ storage benchmark](https://www.spe.org/csp/spe11/)
using [JutulDarcy.jl](https://github.com/sintefmath/JutulDarcy.jl).
The CSP11 module supports cases A, B and C on Cartesian or supplied reservoir
domains. Facies come from the bundled [official geometry](https://github.com/Simulation-Benchmarks/11thSPE-CSP/tree/main/geometries).

## Setup

Open Julia in this repository's environment, then develop the local module:

```julia
using Pkg
Pkg.develop(path = "CSP11")
Pkg.instantiate()
```

## Running cases

```julia
using CSP11, JutulDarcy

# Lab scale: 2.8 × 0.01 × 1.2 m, isothermal at 20 °C, atmospheric top boundary
case_a, name_a = setup_spe11_case((280, 120); case = :a)

# Field scale: one metre nominal thickness in the x/depth plane
case_b, name_b = setup_spe11_case((400, 60); case = :b)

# Three dimensions, warped into the physical C geometry
case_c, name_c = setup_spe11_case((85, 50, 40); case = :c)

result = simulate_reservoir(case_a)
```

Injection uses well and facility models. Case A injects at 1.7e-7 kg/s per
well: the first well operates for five hours, and the second starts after
2.5 hours. Migration continues to 120 hours. A has no initialization period;
B/C have a 1000-year initialization period. Sparse data is reported every ten
minutes for A and every 0.1 year for B/C. Dense data uses hourly snapshots for A,
five-year snapshots for B, and the 26 times listed in Section 4.7 for C.
These settings follow the [official benchmark statement](https://github.com/Simulation-Benchmarks/11thSPE-CSP/blob/main/description/spe_csp11_description.pdf).

An existing `DataDomain` must contain permeability and porosity. If `:satnum`
is absent, facies are inferred from the geometry:

```julia
case_b, name_b = setup_spe11_case(domain; case = :b)
```

`run_spe11.jl` is a complete simulation example with sparse reporting. The
delivery scripts run multiple Cartesian resolutions. Dense reporting accepts
explicit `(report_cell, simulation_cell, weight)` tuples via the `weights`
keyword to `CSP11.write_reporting_grid`. Reporting dimensions use `(nx, ny, nz)`;
indices vary fastest in x, then y, then upward elevation. Weights are normalized
separately for intensive fields and total CO₂ mass. Include the injection-time
state when writing dense data, or pass explicit `times` in seconds since injection.
The sparse hook reads the converged state before Jutul commits the ministep;
both writers interpolate linearly between available times.

## Benchmark comparison and numerical choices

Case definitions live in `CSP11/src/specification.jl`, with equation/table
references. The following describes the default setup, including the native
JutulDarcy equations it uses; custom domains and keyword overrides can change it.

| Requirement | Implementation or approximation |
| --- | --- |
| Geometry, facies, wells and POPs (Sections 2.2, 3.1, 4.2) | Bundled official polygons and physical C arch. Facies, boxes and POPs use cell-centre/nearest-cell sampling. Both POP indices are retained on coarse grids. |
| C permeability (Eqs. 4.3–4.5) | Full tensor push-forward in depth-positive coordinates; signed cross terms are preserved. |
| Rock, relative permeability and capillary pressure (Tables 1, 2, 4, 5; Eqs. 2.6–2.9) | Specified values; capillary pressure is adaptively tabulated. Facies 7 has exactly zero permeability. |
| Facies-7 porosity | Regularized to 0.001 to retain nonsingular fluid equations; the statement specifies zero. Volume diagnostics include this regularization. |
| Gravity and initial pressure (Eq. 2.1, Section 5.5) | 9.81 m/s²; native hydrostatic integration with water-table density and the specified physical pressure datum, independent of grid/well-cell centres. |
| Lateral buffer storage (Sections 3.3, 4.3) | Adds 50,000 m times each eligible face area, including corner faces; restricted to facies 2–5. Section 4.3 inherits B's rules; facies 6 on C's new lateral faces has no explicit assignment, and receives no buffer. |
| Temperature boundaries (Sections 3.3, 4.3) | Geothermal temperatures at physical top/bottom faces. Active injectors impose 10 °C through native thermal boundaries; disabled wells conduct no heat. Conductance uses the native well-index approximation to the circular well boundary. |
| Injection (Sections 2.4, 3.5, 4.5) | Prescribed total mass rates and periods. Native well volumes are scaled by `1e-8` to make their artificial storage negligible; inactive perforations are closed. C defaults to independently controlled perforations weighted by physical trajectory length. `divide_c_wells=false` instead uses two common-pressure well models, whose perforation rates need not be uniform per length. |
| Phase properties and partitioning (Section 2.1) | Native CSP11 NIST/Spycher tables and Garcia water-density mixing. Native PVT/flash evaluation uses the common reservoir pressure, rather than each phase's capillary-adjusted pressure. |
| Molecular diffusion (Eq. 2.3) | Specified phase diffusion coefficients, zero in facies 7. Native JutulDarcy diffusion uses mole-fraction gradients; the statement specifies mass-fraction gradients. |
| Mechanical dispersion (Eq. 2.3) | Not implemented by the native solver. The specified dispersivities, 0.01 m and 10 m, are recorded in the case definitions but do not add a dispersive flux. |
| Fluid thermal conductivity (Section 3.2) | Fixed at 0.6 and 0.088 W/(m K), rather than the stated pressure/temperature-dependent NIST values. |
| Fluid internal energy (Section 3.2) | Native enthalpy mixes pure-component tables by phase composition and derives internal energy as `h - p/ρ`. The statement instead specifies phase internal energy independent of composition. |
| Thermal energy (Eq. 3.1) | Native internal-energy storage, enthalpy advection and conduction; the statement's gravitational potential-energy terms are omitted. Rock heat capacity/density match Table 5. |
| Sparse/dense data (Sections 2.6, 3.7, 4.7) | Includes time zero, correct cadence and A/B/C column sets, and conservative mass remapping. Sparse boundary CO₂ counts the added buffer fraction of each cell, rather than its entire mass. |
| Convection proxy | Intentionally skipped; `M_C` remains `NaN`. |
| Performance data (Section 2.6.3) | Detailed time-series writer uses actual model degrees of freedom when supplied; otherwise emits `n/a`. |

The native transport/property/energy differences above are model limitations,
not adjustments to the benchmark definition. They matter when interpreting
results or claiming equation-level agreement. No separate CSP11 variable types
are introduced to replace those native equations.

## Tests

```julia
using Pkg
Pkg.test("CSP11")
```

## Docker

Build `docker/jutul-spe11csp-juliaimg.Dockerfile` with the desired Julia version,
then run `run_spe11.jl` in the container after developing `CSP11` and
instantiating its environment.
