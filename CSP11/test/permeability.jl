@testset "Signed permeability tensors" begin
    # A signed cross term is valid when the principal permeabilities are positive.
    for (dims, tensor) in (((2, 2), [2.0, -0.5, 1.0]),
            ((2, 2, 2), [2.0, 0.0, 0.0, 2.0, -0.5, 1.0]))
        mesh = CartesianMesh(dims)
        permeability = repeat(1e-13.*tensor, 1, prod(dims))
        domain = reservoir_domain(mesh; permeability)
        @test domain[:permeability] == permeability
        @test minimum(JutulDarcy.permeability_principal_values(permeability, length(dims))) > 0
        invalid = copy(permeability)
        invalid[1, :] .= -1e-13
        @test_throws ArgumentError reservoir_domain(mesh; permeability = invalid)
    end
    @test_throws ArgumentError reservoir_domain(CartesianMesh((2, 2)); permeability = -1e-13)
    # Positive diagonal entries alone do not make a tensor physically valid.
    @test_throws ArgumentError reservoir_domain(CartesianMesh((2, 2)); permeability = [1e-13, -2e-13, 1e-13])
    schedule = (thermal = false, nstep_initialization = 0, nstep_injection1 = 1,
        nstep_injection2 = 1, nstep_migration = 1, use_reporting_steps = false)
    case, = setup_spe11_case((12, 10, 8); case = :c, schedule...)
    domain = reservoir_domain(case.model)
    @test minimum(domain[:permeability][5, :]) < 0
    @test maximum(domain[:permeability][5, :]) > 0
    _, messages = JutulDarcy.CaseValidation.validate(case; info_level = -1)
    @test !any(m -> m[1] != JutulDarcy.CaseValidation.OK && occursin("permeability", m[3]), messages)
    domain[:permeability][6, 1] = -1e-13
    _, messages = JutulDarcy.CaseValidation.validate(case; info_level = -1)
    @test any(m -> m[1] == JutulDarcy.CaseValidation.WARNING && occursin("permeability", m[3]), messages)
end
