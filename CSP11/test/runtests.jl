using CSP11
using Jutul
using JutulDarcy
using Test

@testset "SPE11 geometry" begin
    @test CSP11.spe11_facies(2700.0, 300.0) in 1:7
    @test CSP11.spe11_facies(5100.0, 700.0) in 1:7
    @test CSP11.spe11c_elevation_offset(0.0) ≈ 0.0
    @test CSP11.spe11c_elevation_offset(2500.0) ≈ 155.0
    @test CSP11.spe11c_elevation_offset(5000.0) ≈ 10.0
end

@testset "MAT-independent Cartesian domains" begin
    domain_b = CSP11.setup_spe11_domain((40, 12); case = :b)
    @test number_of_cells(domain_b) == 40*12
    @test extrema(domain_b[:satnum]) == (1, 7)
    @test !isempty(domain_b[:boundary])
    @test length(CSP11.setup_spe11_wells!(domain_b, :b)) == 2

    domain_c = CSP11.setup_spe11_domain((12, 10, 8); case = :c)
    wells_c = CSP11.setup_spe11_wells!(domain_c, :c)
    n1, n2 = domain_c[:num_well_cells]
    @test length(wells_c) == n1 + n2
    @test n1 > 0 && n2 > 0
    @test sum(domain_c[:well_rates][1:n1]) ≈ 50.0
    @test sum(domain_c[:well_rates][n1+1:end]) ≈ 50.0
    @test minimum(domain_c[:permeability][5, :]) < 0
    @test maximum(domain_c[:permeability][5, :]) > 0
end

@testset "Direct SPE11 saturation functions" begin
    @test CSP11.spe11_erf(0.0) == 0.0
    @test CSP11.spe11_erf(0.5) ≈ 0.5204998778130465 atol = 3.1e-8
    @test CSP11.spe11_erf(1.0) ≈ 0.8427007929497149 atol = 3.1e-8
    @test CSP11.spe11_erf(-2.0) ≈ -0.9953222650189527 atol = 3.1e-8

    satnum = collect(1:7)
    kr, pc = CSP11.spe11_saturation_functions(satnum)
    @test kr.regions == satnum
    @test pc.regions == satnum

    expected_half = 0.5^1.5
    @test kr.krg[1](0.55) ≈ expected_half rtol = 5e-4
    @test kr.krog[2](0.57) ≈ expected_half rtol = 5e-4
    @test kr.krg[1](0.10) ≈ 0.0 atol = 1e-15
    @test kr.krog[1](0.32) ≈ 0.0 atol = 1e-15

    entry_pressure_facies_2 = 6.12e-3*sqrt(0.20/1e-13)
    @test pc.pc[1][2](0.0) ≈ entry_pressure_facies_2 rtol = 1e-7
    @test pc.pc[1][2](0.86) ≈ 3e7
    @test sum(length(table.X) for table in pc.pc[1]) < 2000
    @test all(length(table.k.X) < 100 for table in kr.krg)
    @test all(length(table.k.X) < 100 for table in kr.krog)

    for region in 1:6
        immobile = CSP11.spe11_wetting_immobile_saturation[region]
        entry_pressure = 6.12e-3*sqrt(
            (0.10, 0.20, 0.20, 0.20, 0.25, 0.35)[region]/
            (1e-16, 1e-13, 2e-13, 5e-13, 1e-12, 2e-12)[region]
        )
        table = pc.pc[1][region]
        for sg in range(0.0, 1.0; length = 1001)
            exact = CSP11.spe11_capillary_pressure(sg, immobile, entry_pressure)
            @test table(sg) ≈ exact rtol = 1.1e-3 atol = 1.0
        end
    end
end

@testset "Complete Cartesian case" begin
    case_b, name = CSP11.setup_spe11_case((20, 10);
        case = :b,
        thermal = false,
        nstep_initialization = 0,
        nstep_injection1 = 1,
        nstep_injection2 = 1,
        nstep_migration = 1,
        use_reporting_steps = false
    )
    @test length(case_b.dt) == 3
    @test occursin("cartesian_20x10", name)
    @test !isempty(JutulDarcy.well_symbols(case_b.model))
end

@testset "SPE11 reservoir sources" begin
    schedule = (nstep_initialization = 1, nstep_injection1 = 1,
        nstep_injection2 = 1, nstep_migration = 1, use_reporting_steps = false)
    for (specase, dims, rate) in ((:b, (20, 10), 0.035), (:c, (12, 6, 8), 50.0)),
            thermal in (false, true)
        case, name = CSP11.setup_spe11_case(dims; case = specase, thermal,
            use_wells = false, schedule...)
        @test isempty(JutulDarcy.well_symbols(case.model))
        @test keys(case.model.models) == (:Reservoir,)
        @test endswith(name, "_sources")
        @test length(case.dt) == 4
        @test isnothing(case.forces[1][:Reservoir].sources)
        @test isnothing(case.forces[4][:Reservoir].sources)
        first = case.forces[2][:Reservoir].sources
        both = case.forces[3][:Reservoir].sources
        @test sum(s.value for s in first) ≈ rate
        @test sum(s.value for s in both) ≈ 2rate
        @test all(s.temperature == 283.15 for s in both)
        @test allunique(s.cell for s in both)
        domain = reservoir_domain(case.model)
        cells, weights = domain[:spe11_source_cells], domain[:spe11_source_weights]
        @test Set(s.cell for s in first) == Set(cells[1])
        @test [s.value for s in first] ≈ [rate*weights[1][findfirst(==(s.cell), cells[1])] for s in first]
        if thermal
            model = reservoir_model(case.model)
            nc = number_of_cells(domain)
            state = (Pressure = case.state0[:Reservoir][:Pressure],
                ComponentHeatCapacity = case.parameters[:Reservoir][:ComponentHeatCapacity])
            storage = (state = state,)
            mass, energy = zeros(2, nc), zeros(1, nc)
            Jutul.apply_forces_to_equation!(mass, storage, model,
                model.equations[:mass_conservation], nothing, both, 0.0)
            Jutul.apply_forces_to_equation!(energy, storage, model,
                model.equations[:energy_conservation], nothing, both, 0.0)
            @test all(iszero, mass[1, :])
            @test sum(mass[2, :]) ≈ -2rate
            for source in both
                cell, temperature = source.cell, source.temperature
                pressure = state.Pressure[cell]
                density = model.secondary_variables[:PhaseMassDensities].tab(pressure, temperature)[2]
                control = CSP11.rate_to_injection_control(source.value,
                    last(JutulDarcy.reference_densities(model.system)), temperature)
                well_state = (Pressure = state.Pressure,
                    Saturations = repeat([0.0, 1.0], 1, nc),
                    PhaseMassDensities = repeat([1.0, density], 1, nc),
                    ComponentHeatCapacity = state.ComponentHeatCapacity)
                enthalpy = JutulDarcy.well_top_node_enthalpy(control, model, well_state, temperature, cell)
                @test energy[1, cell] ≈ -source.value*enthalpy
            end
        end
    end
    # Dividing the C wells must leave the physical per-cell sources unchanged.
    domain = CSP11.setup_spe11_domain((12, 6, 8); case = :c)
    options = (case = :c, thermal = false, use_wells = false, schedule...)
    unsplit, _ = CSP11.setup_spe11_case(domain; options...)
    divided, _ = CSP11.setup_spe11_case(domain; options..., divide_c_wells = true)
    for step in (2, 3)
        a, b = unsplit.forces[step][:Reservoir].sources, divided.forces[step][:Reservoir].sources
        @test [(s.cell, s.temperature) for s in a] == [(s.cell, s.temperature) for s in b]
        @test [s.value for s in a] ≈ [s.value for s in b]
    end
    # Duplicate cells can arise on coarse grids or overlapping trajectories.
    domain[:spe11_source_cells, nothing] = [[1, 2], [2, 3]]
    domain[:spe11_source_weights, nothing] = [[0.5, 0.5], [0.25, 0.75]]
    combined = CSP11.spe11_source_forces(divided.model, [50.0, 50.0], 283.15)[:Reservoir].sources
    @test [s.cell for s in combined] == [1, 2, 3]
    @test [s.value for s in combined] == [25.0, 37.5, 37.5]
    @test_throws ArgumentError CSP11.spe11_source_forces(divided.model, [-1.0, 0.0], 283.15)
    # MAT/custom wells are DataDomains; source mode extracts their physical
    # perforations while retaining prescribed rates and removing well models.
    domain = CSP11.setup_spe11_domain((20, 10); case = :b)
    domain[:well_cells, nothing] = [20, 30]
    wells = [setup_well(domain, [c]; name = Symbol(:INJ, i-1))
        for (i, c) in enumerate(domain[:well_cells])]
    custom, _ = CSP11.setup_spe11_case(domain; case = :b, wells,
        thermal = false, use_wells = false, schedule...)
    @test keys(custom.model.models) == (:Reservoir,)
    @test [s.cell for s in custom.forces[3][:Reservoir].sources] == [20, 30]
    @test sum(s.value for s in custom.forces[3][:Reservoir].sources) ≈ 0.07
end
