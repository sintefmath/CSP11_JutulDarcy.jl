@testset "SPE11A geometry and properties" begin
    domain = CSP11.setup_spe11_domain((40, 24); case = :a)
    @test number_of_cells(domain) == 960
    @test sum(domain[:volumes]) ≈ 2.8*0.01*1.2
    @test domain[:volumes] == domain[:spe11_physical_volumes]
    @test isempty(domain[:boundary])
    @test extrema(domain[:satnum]) == (1, 7)
    @test all(==(293.15), domain[:temperature])
    k, phi = CSP11.rock_props_from_satnum(collect(1:7), :a)
    @test k[1, 1:6] == [4e-11, 5e-10, 1e-9, 2e-9, 4e-9, 1e-8]
    @test k[1, :] == k[2, :] == k[3, :]
    @test phi[1:6] == [0.44, 0.43, 0.44, 0.45, 0.43, 0.46]
    @test domain[:porosity] == phi[domain[:satnum]]
    @test all(==(1.6e-5), domain[:diffusion][2, domain[:satnum] .!= 7])
    @test all(iszero, domain[:diffusion][:, domain[:satnum] .== 7])
    @test domain[:satnum] == CSP11.setup_spe11_domain((40, 1, 24); case = :a)[:satnum]
    for (x, elevation) in ((0.9, 0.3), (1.7, 0.7), (1.5, 0.5))
        @test spe11_facies(x, elevation; case = :a) == spe11_facies(3000x, 1000elevation)
        @test spe11_facies(x, 0.005, 1.2 - elevation; case = :a) == spe11_facies(x, elevation; case = :a)
    end
    cc = domain[:cell_centroids]
    @test findall(==(1.0), domain[:A]) == findall(i -> cc[1, i] >= 1.1 && cc[3, i] >= 0.6, axes(cc, 2))
    @test findall(==(1.0), domain[:B]) == findall(i -> cc[1, i] <= 1.1 && cc[3, i] <= 0.6, axes(cc, 2))
    @test count(>(0), domain[:observation_points]) == 2
    @test_throws ArgumentError CSP11.setup_spe11_domain((10, 2, 10); case = :a)
    @test_throws ArgumentError CSP11.setup_spe11_case((10, 10); case = :a, thermal = true)
end

@testset "SPE11A injection" begin
    case, = setup_spe11_case((28, 12); case = :a, use_reporting_steps = false,
        nstep_injection1 = 1, nstep_injection2 = 1, nstep_migration = 1)
    smoke = JutulCase(case.model, [60.0], [case.forces[1]];
        state0 = case.state0, parameters = case.parameters)
    result = simulate_reservoir(smoke; info_level = -1, validate = false,
        max_timestep = 60.0, max_nonlinear_iterations = 20)
    @test length(result.states) == 1
    state = only(result.states)
    @test all(isfinite, state[:Pressure])
    @test all(s -> -1e-8 <= s <= 1 + 1e-8, state[:Saturations])
    @test sum(state[:TotalMasses][2, :]) > 0
end

@testset "SPE11A saturation functions" begin
    kr, pc = CSP11.spe11_saturation_functions(collect(1:7); case = :a)
    @test kr.krg[1](0.55) ≈ 0.25
    @test kr.krog[2](0.57) ≈ 0.25
    for (i, entry_pressure) in enumerate((1500.0, 300.0, 100.0, 25.0, 10.0, 1.0))
        @test pc.pc[1][i](0.0) ≈ entry_pressure rtol = 1e-4
        @test pc.pc[1][i](1.0 - CSP11.spe11_wetting_immobile_saturation[i]) ≈ 9.5e4
    end
end

@testset "Complete SPE11A case" begin
    case, name = @test_logs min_level = Logging.Warn setup_spe11_case((40, 24); case = :a)
    domain = reservoir_domain(case.model)
    @test occursin("spe11a_cartesian_40x24_isothermal", name)
    @test !JutulDarcy.model_is_thermal(reservoir_model(case.model))
    @test length(case.dt) == 720
    @test all(==(600.0), case.dt)
    @test sum(case.dt) == 432000.0
    @test Set(JutulDarcy.well_symbols(case.model)) == Set((:INJ0, :INJ1))
    state = case.state0[:Reservoir]
    @test all(==(293.15), case.parameters[:Reservoir][:Temperature])
    @test all(==(1.0), state[:OverallMoleFractions][1, :])
    pressure = state[:Pressure]
    density = reservoir_model(case.model).secondary_variables[:PhaseMassDensities].tab(1.1e5, 293.15)[1]
    @test pressure ≈ 1.1e5 .+ density*9.81.*domain[:cell_centroids][3, :] rtol = 1e-6
    _, rho_co2 = JutulDarcy.reference_densities(reservoir_model(case.model).system)
    controls(step) = case.forces[step][:Facility].control
    @test controls(1)[:INJ0].target.value*rho_co2 ≈ 1.7e-7
    @test controls(15)[:INJ1] isa DisabledControl
    @test controls(16)[:INJ1].target.value*rho_co2 ≈ 1.7e-7
    @test controls(30)[:INJ0].temperature == 293.15
    @test all(c isa DisabledControl for c in values(controls(31)))
    faces = domain[:spe11_pressure_boundary_faces]
    @test length(faces) == 40
    @test all(iszero, domain[:boundary_centroids][3, faces])
    trans = compute_boundary_trans(domain)
    for forces in (case.forces[1], case.forces[16], case.forces[end])
        bc = forces[:Reservoir].bc
        @test [b.cell for b in bc] == domain[:boundary_neighbors][faces]
        @test all(b.pressure ≈ pressure[b.cell] for b in bc)
        @test all(b.temperature == 293.15 && b.fractional_flow == (1.0, 0.0) && b.density == density for b in bc)
        @test [b.trans_flow for b in bc] ≈ trans[faces]
    end
    # The atmospheric top boundary must preserve an initially resting column.
    nc = number_of_cells(domain)
    boundary_state = (Pressure = pressure, PhaseMobilities = ones(2, nc),
        PhaseMassDensities = repeat([density, 1.0], 1, nc))
    rmodel = reservoir_model(case.model)
    for bc in case.forces[1][:Reservoir].bc
        @test all(iszero, JutulDarcy.compute_bc_mass_fluxes(rmodel.system, bc,
            Jutul.global_map(rmodel), boundary_state))
    end
    custom = reservoir_domain(UnstructuredMesh(CartesianMesh((40, 1, 24), (2.8, 0.01, 1.2)));
        permeability = domain[:permeability], porosity = domain[:porosity])
    short, = setup_spe11_case(custom; case = :a, use_reporting_steps = false,
        nstep_injection1 = 1, nstep_injection2 = 1, nstep_migration = 1)
    @test short.dt == [9000.0, 9000.0, 414000.0]
    @test reservoir_domain(short.model)[:satnum] == domain[:satnum]
end
