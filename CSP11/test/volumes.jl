@testset "SPE11 buffer storage" begin
    # Non-unit cell widths expose a volume-multiplier implementation. Corner
    # cells have one eligible x face in B and both x and y faces in C.
    for specase in (:b, :c)
        dims = specase == :b ? (2, 1, 2) : (2, 2, 2)
        sizes = specase == :b ? (8.0, 1.0, 4.0) : (8.0, 6.0, 4.0)
        nc = prod(dims)
        tags = fill(2, nc)
        tags[end-1:end] = [1, 7]
        domain = reservoir_domain(CartesianMesh(dims, sizes);
            satnum = tags, permeability = 1e-13, porosity = 0.2)
        CSP11.prepare_spe11_domain!(domain, specase)
        volume = specase == :b ? 8.0 : 24.0
        buffer = specase == :b ? 100000.0 : 700000.0
        @test domain[:spe11_physical_volumes] ≈ fill(volume, nc)
        @test domain[:volumes][1:end-2] ≈ fill(volume + buffer, nc - 2)
        @test domain[:volumes][end-1:end] ≈ fill(volume, 2)
        @test domain[:boundary] == collect(1:nc-2)
        @test domain[:porosity][end] == 0.001
        @test all(==(0.2), domain[:porosity][1:end-1])
        @test pore_volume(domain)[end] ≈ 0.001*volume
        original_volumes = copy(domain[:volumes])
        original_heat_capacity = copy(domain[:rock_heat_capacity])
        CSP11.prepare_spe11_domain!(domain, specase)
        @test domain[:volumes] == original_volumes
        @test domain[:rock_heat_capacity] == original_heat_capacity
        # Existing domains from the old setup may carry inflated volumes and
        # a boundary marker. Recover physical volumes from the mesh first.
        legacy = reservoir_domain(CartesianMesh(dims, sizes);
            satnum = tags, permeability = 1e-13, porosity = 0.05)
        legacy[:boundary, nothing] = collect(1:nc)
        legacy[:volumes] .*= 50001.0
        CSP11.prepare_spe11_domain!(legacy, specase)
        @test legacy[:volumes] ≈ original_volumes
        @test legacy[:porosity][end] == 0.001
    end
end

@testset "SPE11 integrated volume references" begin
    b = CSP11.spe11_reference_volumes(:b)
    c = CSP11.spe11_reference_volumes(:c)
    a = CSP11.spe11_reference_volumes(:a)
    @test a.physical_bulk ≈ 2.8*0.01*1.2
    @test a.bulk == a.physical_bulk
    @test a.pore == a.physical_pore
    @test a.buffer_bulk == a.buffer_pore == 0.0
    @test b.physical_bulk ≈ 8400.0*1200.0
    @test c.physical_bulk ≈ 8400.0*5000.0*1200.0
    @test c.physical_pore ≈ 5000.0*b.physical_pore
    # Independent reference totals from integration of the published polygons.
    @test b.bulk ≈ 9.383552519e7
    @test b.pore ≈ 2.051114471635e7
    @test c.bulk ≈ 1.1615512576601074e12
    @test c.pore ≈ 2.6039960975358762e11
    @test_throws ArgumentError CSP11.spe11_reference_volumes(:d)
end

@testset "SPE11 volume diagnostics" begin
    schedule = (thermal = false, nstep_initialization = 0,
        nstep_injection1 = 1, nstep_injection2 = 1, nstep_migration = 1,
        use_reporting_steps = false)
    for (specase, dims) in ((:b, (20, 10)), (:c, (12, 6, 8)))
        case, _ = @test_logs min_level = Logging.Warn CSP11.setup_spe11_case(dims;
            case = specase, schedule...)
        domain = reservoir_domain(case.model)
        @test all(==(0.001), domain[:porosity][domain[:satnum] .== 7])
        report = @test_logs min_level = Logging.Warn CSP11.check_spe11_volumes(
            case.model, case.parameters, specase)
        @test report.physical ≈ (specase == :b ? 1.008e7 : 5.04e10)
        @test report.bulk ≈ sum(domain[:volumes])
        @test report.pore ≈ sum(pore_volume(domain))
        # A real setup inconsistency must still warn; coarse facies sampling
        # alone must not be reported as an incorrect volume.
        domain[:spe11_physical_volumes] .*= 1.1
        @test_logs (:warn, r"physical bulk volume") (:warn, r"effective bulk volume") min_level = Logging.Warn CSP11.check_spe11_volumes(
            case.model, case.parameters, specase)
        domain[:spe11_physical_volumes] ./= 1.1
        domain[:volumes] .*= 1.1
        @test_logs (:warn, r"effective bulk volume") (:warn, r"model pore volume") min_level = Logging.Warn CSP11.check_spe11_volumes(
            case.model, case.parameters, specase)
        domain[:volumes] ./= 1.1
        case.parameters[:Reservoir][:FluidVolume] .*= 1.1
        @test_logs (:warn, r"model pore volume") min_level = Logging.Warn CSP11.check_spe11_volumes(
            case.model, case.parameters, specase)
    end
end
