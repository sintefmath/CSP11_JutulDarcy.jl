@testset "SPE11 fixed face temperatures" begin
    for (specase, dims) in ((:b, (20, 10)), (:c, (12, 6, 8)))
        for use_wells in (true, false)
            case, _ = CSP11.setup_spe11_case(dims; case = specase,
                use_wells, thermal = true, nstep_initialization = 1,
                nstep_injection1 = 1, nstep_injection2 = 1,
                nstep_migration = 1, use_reporting_steps = false)
            domain = reservoir_domain(case.model)
            @test all(==(850.0), domain[:rock_heat_capacity])
            @test all(==(850.0), case.parameters[:Reservoir][:RockHeatCapacity])
            faces = domain[:spe11_thermal_boundary_faces]
            cells = domain[:boundary_neighbors][faces]
            for forces in case.forces
                bc = forces[:Reservoir].bc
                @test Set(b.cell for b in bc) == Set(cells)
                @test allunique(b.cell for b in bc)
                @test all(b.trans_flow == 0.0 && b.trans_thermal > 0.0 for b in bc)
                for b in bc
                    face = faces[findfirst(==(b.cell), cells)]
                    depth = domain[:boundary_centroids][3, face]
                    @test b.temperature ≈ 343.15 - 0.025*(1200.0 - depth)
                    # The fixed temperature belongs at the face, not at the
                    # centre of a cell with an artificially large heat capacity.
                    @test b.temperature != case.state0[:Reservoir][:Temperature][b.cell]
                end
            end
        end
        # Isothermal cases carry no thermal forcing.
        case, _ = CSP11.setup_spe11_case(dims; case = specase,
            use_wells = false, thermal = false, nstep_initialization = 0,
            nstep_injection1 = 1, nstep_injection2 = 1,
            nstep_migration = 1, use_reporting_steps = false)
        @test all(isnothing(f[:Reservoir].bc) for f in case.forces)
    end
    # Top and bottom may belong to the same cell. Sum their conductances and
    # use the conductance-weighted temperature, avoiding parallel write races.
    domain = reservoir_domain(CartesianMesh((2, 1, 1), (8.0, 1.0, 4.0));
        satnum = [7, 7], permeability = 1e-20, porosity = 0.001,
        rock_thermal_conductivity = 2.0, rock_heat_capacity = 850.0)
    CSP11.prepare_spe11_domain!(domain, :b)
    model, = CSP11.setup_reservoir_model_csp11(domain; thermal = true)
    bc = CSP11._spe11_thermal_boundary_conditions(model)
    @test length(bc) == 2
    @test all(b.temperature ≈ 313.20 for b in bc)
    @test all(b.trans_thermal ≈ 2*4/2*((1 - 0.001)*2.0 + 0.001*0.6) for b in bc)
    # Legacy prepared domains are restored once, preserving physical capacity.
    domain[:rock_heat_capacity] .*= 1e5
    domain[:spe11_fixed_temperature_boundaries, nothing] = true
    CSP11.prepare_spe11_domain!(domain, :b)
    @test all(==(850.0), domain[:rock_heat_capacity])
    CSP11.prepare_spe11_domain!(domain, :b)
    @test all(==(850.0), domain[:rock_heat_capacity])
end
