using CSP11
using Jutul
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
end
