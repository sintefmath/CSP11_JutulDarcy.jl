@testset "Benchmark statement parameters and schedules" begin
    @test CSP11.spe11_reporting_times(:a)[end] == 432000.0
    @test length(CSP11.spe11_reporting_times(:a)) == 121
    @test length(CSP11.spe11_reporting_times(:b)) == 201
    expected_c = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 75, 100,
        150, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000]
    @test CSP11.spe11_reporting_times(:c)./31536000 == expected_c
    for case in (:b, :c)
        times = CSP11.spe11_reporting_times(case; kind = :sparse)
        @test length(times) == 10001
        @test all(==(3153600.0), diff(times))
        spec = CSP11.spe11_case_specification(case)
        @test spec.permeability == (1e-16, 1e-13, 2e-13, 5e-13, 1e-12, 2e-12, 0.0)
        @test spec.porosity[1:6] == (0.10, 0.20, 0.20, 0.20, 0.25, 0.35)
        @test spec.rock_conductivity == (1.9, 1.25, 1.25, 1.25, 0.92, 0.26, 2.0)
        @test spec.diffusion == (1e-9, 2e-8)
        @test spec.buffer_facies == 2:5
        @test spec.dispersivity == 10.0
    end
    c, = setup_spe11_case((12, 6, 8); case = :c, thermal = false)
    @test length(c.dt) == 26
    @test c.dt[1] == 1000*31536000.0
    @test cumsum(c.dt)[2:end]./31536000 .- 1000 == expected_c[2:end]
    domain = reservoir_domain(c.model)
    @test all(iszero, domain[:permeability][:, domain[:satnum] .== 7])
    @test all(iszero, domain[:spe11_buffer_volumes][domain[:satnum] .== 6])
    @test length(JutulDarcy.well_symbols(c.model)) > 2
    @test sum(domain[:well_rates]) ≈ 100.0
    @test_throws ArgumentError CSP11.setup_reservoir_forces_and_timesteps_csp11(c.model, :c; nstep_injection1 = 0)
    @test_throws ArgumentError CSP11.prepare_spe11_domain!(domain, :b)
    # Both POPs remain usable if a coarse grid maps them to the same cell.
    coarse = setup_spe11_domain((1, 1); case = :a)
    @test coarse[:spe11_observation_cells] == (1, 1)
    # Eq. 4.4 is a tensor push-forward, not a rotation or entrywise magnitude.
    diagonal = repeat([2e-13, 3e-13, 4e-14], 1, 3)
    for (i, y) in enumerate((0.0, 2500.0, 5000.0))
        transformed = CSP11.transform_spe11c_permeability(diagonal[:, i:i], [y])
        slope = -3/25*((y - 2500.0)/2500.0) + 1/500
        f = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 -slope 1.0]
        k = [2e-13 0.0 0.0; 0.0 3e-13 0.0; 0.0 0.0 4e-14]
        @test Jutul.expand_perm(transformed[:, 1], 3) ≈ f*k*f'
    end
end

@testset "Specified gravity and hydrostatic pressure datum" begin
    for (specase, dims) in ((:b, (20, 10)), (:b, (40, 24)), (:c, (12, 6, 8)))
        domain = setup_spe11_domain(dims; case = specase)
        model, parameters = setup_reservoir_model_csp11(domain)
        table = CSP11.spe11_hydrostatic_pressure_table(model, specase)
        @test table(900.0) ≈ 3e7
        rho = reservoir_model(model).secondary_variables[:PhaseMassDensities].tab(3e7, 335.65)[1]
        @test (table(901.0) - table(899.0))/2 ≈ 9.81*rho rtol = 1e-3
        z = domain[:cell_centroids][3, :]
        left, right = eachrow(domain[:neighbors])
        @test parameters[:Reservoir][:TwoPointGravityDifference] ≈ -9.81.*(z[right] - z[left])
        @test domain[:temperature] ≈ 343.15 .- 0.025.*(1200 .- z)
    end
    @test CSP11.spe11_geothermal_temperature(-155.0) - 273.15 ≈ 36.125
end

@testset "Conservative and validated reporting maps" begin
    state = Dict(:Pressure => [1.0, 3.0], :Temperature => [293.15, 293.15],
        :Saturations => [1.0 1.0; 0.0 0.0],
        :LiquidMassFractions => [1.0 1.0; 0.0 0.0], :VaporMassFractions => [0.0 0.0; 1.0 1.0],
        :PhaseMassDensities => [1000.0 1000.0; 2.0 2.0], :TotalMasses => [0.0 0.0; 6.0 12.0])
    weights = [(1, 1, 10.0), (2, 1, 5.0), (2, 2, 5.0), (3, 2, 10.0)]
    mapping = CSP11.reporting_mapping(weights, 2, (3, 1, 1))
    result = CSP11.map_to_reporting_grid(state, mapping)
    @test result[:p] == [1.0, 2.0, 3.0]
    @test sum(result[:co2_mass]) == 18.0
    @test_throws ArgumentError CSP11.reporting_mapping([(1, 1, -1.0)], 1, (1, 1, 1))
    @test_throws ArgumentError CSP11.reporting_mapping([(2, 1, 1.0)], 1, (1, 1, 1))
    @test_throws ArgumentError CSP11.reporting_mapping([(1, 1, 1.0)], 2, (1, 1, 1))
    mktempdir(@__DIR__) do path
        case = (model = Dict(:Reservoir => (domain = CartesianMesh((2, 1, 1)),)), dt = [3600.0])
        @test_throws ArgumentError CSP11.write_reporting_grid(case, [state], path, :a;
            weights = weights, dims = (3, 1, 1), time_offset = 0.0)
        CSP11.write_reporting_grid(case, [state, state], path, :a;
            weights = weights, dims = (3, 1, 1), time_offset = 0.0)
        @test isfile(joinpath(path, "spe11a_spatial_map_0h.csv"))
        @test length(split(readlines(joinpath(path, "spe11a_spatial_map_1h.csv"))[2], ',')) == 9
        # C includes y even for a single y layer, and uses all 26 specified times.
        c = (model = case.model, dt = [1000*31536000.0])
        map_c = [(1, 1, 1.0), (2, 2, 1.0)]
        CSP11.write_reporting_grid(c, [state, state], path, :c;
            weights = map_c, dims = (2, 1, 1), time_offset = 0.0)
        files = filter(file -> startswith(file, "spe11c_spatial_map"), readdir(path))
        @test length(files) == 26
        @test !isfile(joinpath(path, "spe11c_spatial_map_55y.csv"))
        @test length(split(readlines(joinpath(path, "spe11c_spatial_map_75y.csv"))[2], ',')) == 11
    end
end

@testset "Thermal injection with prescribed well temperature" begin
    case, = setup_spe11_case((20, 10); case = :b, thermal = true,
        nstep_initialization = 0, nstep_injection1 = 1, nstep_injection2 = 1,
        nstep_migration = 1, use_reporting_steps = false)
    smoke = JutulCase(case.model, fill(3600.0, 3), case.forces;
        state0 = case.state0, parameters = case.parameters)
    @test case.forces[1][:INJ0].mask.values == [1.0]
    @test case.forces[1][:INJ1].mask.values == [0.0]
    @test case.forces[end][:INJ0].mask.values == [0.0]
    @test case.forces[end][:INJ1].mask.values == [0.0]
    mktempdir(@__DIR__) do path
        domain = reservoir_domain(case.model)
        hook, csv = CSP11.get_reporting_hook(path, domain; interval = 3600.0)
        sim, config = setup_reservoir_simulator(smoke; info_level = -1,
            max_timestep = 3600.0, max_nonlinear_iterations = 20, post_ministep_hook = hook)
        result = simulate_reservoir(smoke; simulator = sim, config = config, validate = false)
        @test length(result.states) == 3
        @test all(state -> all(isfinite, state[:Temperature]), result.states)
        @test Jutul.value(sim.storage.state0.INJ0.Temperature[1]) ≈ 283.15 atol = 1e-3
        lines = readlines(csv)
        @test length(lines) == 5
        @test parse(Float64, first(split(lines[end], ','))) == 10800.0
        pop1 = domain[:spe11_observation_cells][1]
        @test parse(Float64, split(lines[end], ',')[2]) ≈ last(result.states)[:Pressure][pop1] rtol = 0.0 atol = 1e-8
    end
end
