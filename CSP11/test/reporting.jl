@testset "Reporting without grid-file metadata" begin
    state = Dict(
        :Pressure => [1.0, 3.0], :Temperature => [293.15, 293.15],
        :Saturations => [0.75 0.5; 0.25 0.5],
        :LiquidMassFractions => [1.0 1.0; 0.0 0.0],
        :VaporMassFractions => [0.0 0.0; 1.0 1.0],
        :PhaseMassDensities => [1000.0 1000.0; 2.0 2.0],
        :TotalMasses => [10.0 20.0; 6.0 12.0]
    )
    weights = [(1, 1, 1.0), (2, 1, 0.5), (2, 2, 0.5), (3, 2, 1.0)]
    mapping = CSP11.reporting_mapping(weights, 2, (3, 1, 1))
    mapped = CSP11.map_to_reporting_grid(state, mapping)
    @test mapped[:p] == [1.0, 2.0, 3.0]
    @test mapped[:co2_mass] ≈ [4.0, 6.0, 8.0]
    @test sum(mapped[:co2_mass]) ≈ sum(state[:TotalMasses][2, :])

    domain = setup_spe11_domain((40, 24); case = :a)
    model, = setup_reservoir_model_csp11(domain)
    nc = number_of_cells(domain)
    report_state = (
        Pressure = fill(1.1e5, nc), RelativePermeabilities = zeros(2, nc),
        PhaseMassDensities = repeat([1000.0, 2.0], 1, nc),
        LiquidMassFractions = repeat([1.0, 0.0], 1, nc),
        VaporMassFractions = repeat([0.0, 1.0], 1, nc),
        Saturations = repeat([1.0, 0.0], 1, nc), FluidVolume = pore_volume(domain),
        TotalMasses = zeros(2, nc)
    )
    mktempdir(@__DIR__) do path
        grid_case = (model = Dict(:Reservoir => (domain = CartesianMesh((2, 1, 1)),)), dt = [3600.0])
        CSP11.write_reporting_grid(grid_case, [state, state], path, :a;
            weights = weights, dims = (3, 1, 1), time_offset = 0.0)
        dense = readlines(joinpath(path, "spe11a_spatial_map_1h.csv"))
        @test length(dense) == 4
        @test sum(parse(Float64, split(line, ',')[9]) for line in dense[2:end]) ≈ 18.0
        hook, csv, = CSP11.get_reporting_hook(path, domain)
        current = merge(report_state, (Pressure = report_state.Pressure .+ 600.0,))
        simulator = (model = model, storage = (state0 = (Reservoir = report_state,), state = (Reservoir = current,)))
        hook(false, Dict(:success => true, :dt => 600.0), simulator, 600.0, nothing, 1, nothing)
        lines = readlines(csv)
        @test basename(csv) == "spe11a_time_series.csv"
        @test length(lines) == 3
        @test parse(Float64, first(split(lines[2], ','))) == 0.0
        @test parse(Float64, first(split(lines[3], ','))) == 600.0
        @test parse(Float64, split(lines[3], ',')[2]) == 1.1e5 + 600.0
        @test length(split(lines[1], ',')) == length(split(lines[2], ',')) == 13
        # A failed ministep must neither advance the reporting clock nor write.
        hook(false, Dict(:success => false, :dt => 750.0), simulator, 750.0, nothing, 1, nothing)
        @test length(readlines(csv)) == 3
        later = merge(report_state, (Pressure = report_state.Pressure .+ 1350.0,))
        simulator = (storage = (state0 = (Reservoir = current,), state = (Reservoir = later,)),)
        hook(false, Dict(:success => true, :dt => 750.0), simulator, 750.0, nothing, 1, nothing)
        last_line = split(last(readlines(csv)), ',')
        @test parse(Float64, last_line[1]) == 1200.0
        @test parse(Float64, last_line[2]) ≈ 1.1e5 + 1200.0
        # Only effective buffer storage belongs in boundary CO₂, not the full
        # mass of a coarse cell that also contains ordinary physical pore volume.
        bdomain = setup_spe11_domain((20, 10); case = :b)
        nb = number_of_cells(bdomain)
        bstate = (Pressure = fill(3e7, nb), RelativePermeabilities = zeros(2, nb),
            PhaseMassDensities = ones(2, nb), LiquidMassFractions = zeros(2, nb),
            VaporMassFractions = zeros(2, nb), Saturations = zeros(2, nb),
            FluidVolume = ones(nb), TotalMasses = [zeros(1, nb); ones(1, nb)])
        helper = CSP11.CSP11ReportHelper(bdomain)
        values = CSP11.reporting_values(helper, bstate)
        @test values[end] ≈ sum(bdomain[:spe11_buffer_volumes]./bdomain[:volumes])
        @test values[end] < length(bdomain[:boundary])
        hook_b, csv_b = CSP11.get_reporting_hook(path, bdomain; time_offset = 0.0)
        sim_b = (storage = (state0 = (Reservoir = bstate,), state = (Reservoir = bstate,)),)
        hook_b(false, Dict(:success => true, :dt => 31536000.0), sim_b, 31536000.0, nothing, 1, nothing)
        lines_b = readlines(csv_b)
        @test length(lines_b) == 12 # header plus 0, 0.1, ..., 1.0 years
        @test parse(Float64, first(split(lines_b[3], ','))) == 3153600.0
        stats = (newtons = 1, linearizations = 2, linear_iterations = 3,
            linear_setup = 0.1, linear_solve = 0.2)
        report = Dict(:total_time => 1.0, :ministeps => [
            (success = true, dt = 10.0, stats = stats),
            (success = false, dt = 5.0, stats = stats)])
        CSP11.print_performance(path, [report_state], [report], :a; model = model)
        performance = split(last(readlines(joinpath(path, "spe11a_performance_time_series_detailed.csv"))), ',')
        @test parse(Float64, performance[1]) == 10.0
        @test parse(Float64, performance[3]) == 1.0
        @test parse(Float64, performance[5]) == Jutul.number_of_degrees_of_freedom(model)
        @test parse(Float64, performance[6]) == 2.0
        @test parse(Float64, performance[8]) == 6.0
    end
end
