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

@testset "Complete Cartesian case" begin
    case_b, name = CSP11.setup_spe11_case((20, 10);
        case = :b,
        thermal = false,
        include_satfun = false,
        nstep_initialization = 0,
        nstep_injection1 = 1,
        nstep_injection2 = 1,
        nstep_migration = 1,
        use_reporting_steps = false
    )
    @test length(case_b.dt) == 3
    @test occursin("cartesian_20x10", name)
end
