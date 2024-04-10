# Exactly 365 days in a year in spec
const spe11_year = 365*si_unit(:day)


function setup_spe11_caseb_from_mrst_grid(basename;
        thermal = true,
        nstep_initialization = thermal*10,
        nstep_injection1 = 50,
        nstep_injection2 = 50,
        nstep_migration = 100,
        composite = false,
        kwarg...
    )
    dirname = joinpath(@__DIR__, "..", "..", "data")

    pth = joinpath(dirname, "$basename.mat")

    domain, wells = reservoir_domain_and_wells_csp11(pth);
    composite = composite || thermal

    if thermal
        othername = "thermal_cv"
    else
        othername = "isothermal"
    end
    name = "spe11b_$(basename)_$othername"

    model, parameters = setup_reservoir_model_csp11(domain;
        wells = wells,
        thermal = thermal,
        backend = :csr,
        general_ad = false,
        dT_max_abs = 30.0,
        composite = composite,
        split_wells = true,
        kwarg...
    );
    ##
    forces, dt = CSP11.setup_reservoir_forces_and_timesteps_csp11b(model,
        nstep_initialization = nstep_initialization,
        # nstep_migration = 0,
        # nstep_injection1 = 0,
        # nstep_injection2 = 0,
        # time_injection1 = 365*si_unit(:day),
        nstep_injection1 = nstep_injection1,
        nstep_injection2 = nstep_injection2,
        nstep_migration = nstep_migration,
    );
    state0 = setup_state0_csp11(model)

    case = JutulCase(model, dt, forces; state0 = state0, parameters = parameters)
    return (case, name)
end

function reservoir_domain_and_wells_csp11(pth::AbstractString, case = :b; kwarg...)
    matdata = MAT.matread(pth)
    raw_G = matdata["G"]
    buffer_cells = Int.(vec(raw_G["bufferCells"]))
    G = UnstructuredMesh(MRSTWrapMesh(raw_G), z_is_depth = true)
    raw_rock = matdata["rock"]
    K = collect(raw_rock["perm"]')
    @. K = max(K, 1e-10*si_unit(:darcy))
    poro = collect(vec(raw_rock["poro"]))
    poro[poro .< 0.05] .= 0.05
    satnum = Int.(vec(raw_rock["regions"]["saturation"]))
    domain = reservoir_domain_csp11(G, case; satnum = satnum, permeability = K, porosity = poro, kwarg...)
    # domain[:volumes] .= raw_G["cells"]["volumes"]
    @. domain[:volumes][buffer_cells] *= raw_G["bufferMult"]
    cc = domain[:cell_centroids]
    z = cc[3, :]
    z_mid = median(z)
    top_cells = Int[]
    bottom_cells = Int[]
    for i in 1:number_of_boundary_faces(G)
        N = domain[:boundary_normals][:, i]
        if abs(N[3]) > abs(N[1]) + abs(N[2])
            c = domain[:boundary_neighbors][i]
            if domain[:boundary_centroids][3, i] > z_mid
                push!(bottom_cells, c)
            else
                push!(top_cells, c)
            end
        end
    end
    domain[:rock_heat_capacity][top_cells] *= 1e5
    domain[:rock_heat_capacity][bottom_cells] *= 1e5

    wc1, wc2 = Int.(vec(raw_G["cells"]["wellCells"]))
    simple_well = true
    I0 = setup_well(domain, [wc1], simple_well = simple_well, name = :INJ0)
    I1 = setup_well(domain, [wc2], simple_well = simple_well, name = :INJ1)
    wells = [I0, I1]

    A = raw_G["cells"]["fractionInA"]
    B = raw_G["cells"]["fractionInB"]
    C = raw_G["cells"]["fractionInC"]

    boundary = Int.(vec(raw_G["bufferCells"]))
    domain[:A, nothing] = A
    domain[:B, nothing] = B
    domain[:C, nothing] = C
    domain[:boundary, nothing] = boundary

    domain[:well_cells, nothing] = [wc1, wc2]
    return domain, wells
end

function reservoir_domain_csp11(G, case = :b; satnum, temperature = 333.15, kwarg...)
    maximum(satnum) == 7 || throw(ArgumentError("Must have 7 as lowest SATNUM region"))
    minimum(satnum) == 1 || throw(ArgumentError("Must have 1 as lowest SATNUM region"))
    nc = number_of_cells(G)
    length(satnum) == nc || throw(ArgumentError("satnum must have number of cells entries ($nc), was $(length(satnum))"))
    domain = reservoir_domain(G; satnum = satnum, kwarg...)
    if case == :b
        rock_thermal_conductivity = fill(0.85, nc)
        diffusion = repeat([1e-9, 2e-8], 1, nc)
        # At approximate 100 bar, 30 deg C. Dependence comes in later.
        fluid_thermal_conductivity = repeat([0.6, 0.088], 1, nc)
        rock_density = fill(2500, nc)
        c_r = [1.9, 1.25, 1.25, 1.25, 0.92, 0.26, 2.0]
        for (i, reg) in enumerate(satnum)
            rock_thermal_conductivity[i] = c_r[reg]
            if reg == 7
                diffusion[:, i] .= 0.0
            end
        end
        domain[:diffusion] = diffusion
        domain[:rock_thermal_conductivity, Cells()] = rock_thermal_conductivity
        domain[:fluid_thermal_conductivity, Cells()] = fluid_thermal_conductivity
        domain[:rock_density, Cells()] = rock_density
        domain[:component_heat_capacity, Cells()] = repeat([4100.0, 950.0], 1, nc)
        domain[:rock_heat_capacity, Cells()] = fill(850.0, nc)
    else
        throw(ArgumentError("Only case b is supported at the moment."))
    end
    domain[:temperature, Cells()] = temperature
    return domain
end

function setup_reservoir_model_csp11(reservoir::DataDomain; include_satfun = true, kwarg...)
    model, parameters = setup_reservoir_model(reservoir, :co2brine; kwarg...)
    if include_satfun
        pth_data = joinpath(@__DIR__, "..", "..", "small_pyopm", "130_62", "CSP11B_DISGAS.DATA")
        case_cpgrid = setup_case_from_data_file(pth_data)
        model_cp = case_cpgrid.model
        # TODO: These have not been added manually yet. Remove once replaced with analytical versions
        kr = deepcopy(model_cp[:Reservoir][:RelativePermeabilities])
        pc = deepcopy(model_cp[:Reservoir][:CapillaryPressure])

        empty!(kr.regions)
        empty!(pc.regions)
        for c in reservoir[:satnum]
            push!(kr.regions, c)
            push!(pc.regions, c)
        end

        if reservoir_model(model).system isa CompositeSystem
            kr = Pair(:flow, kr)
            pc = Pair(:flow, pc)
        end

        set_secondary_variables!(model[:Reservoir],
            RelativePermeabilities = kr,
            CapillaryPressure = pc
        )
    end
    return (model, parameters)
end

function setup_reservoir_forces_and_timesteps_csp11(model; case = :b, kwarg...)
    if case == :b
        out = setup_reservoir_forces_csp11b(model; kwarg...)
    else
        throw(ArgumentError("Only case b supported at the moment."))
    end
    return out
end

function setup_reservoir_forces_and_timesteps_csp11b(model, case = :b;
        well_labels = (:INJ0, :INJ1),
        nstep_initialization = 1,
        nstep_injection1 = 1,
        nstep_injection2 = 1,
        nstep_migration = 1,
        time_initialization = 1000*spe11_year,
        time_injection1 = 25*spe11_year,
        time_injection2 = 25*spe11_year,
        time_migration = 1000*spe11_year - time_injection1 - time_injection2,
        rate_injection1 = (0.035, 0.0).*si_unit(:kilogram)./si_unit(:second),
        rate_injection2 = (0.035, 0.035).*si_unit(:kilogram)./si_unit(:second),
        injection_temperature = convert_to_si(10, :Celsius)
    )
    tables = JutulDarcy.CO2Properties.co2_brine_property_tables()
    H_tab = tables[:enthalpy]
    co2_H_eval(p, T) = H_tab(p, T)[2]
    H_well = co2_H_eval
    # H_well = missing

    rmodel = JutulDarcy.reservoir_model(model)
    rho_brine, rho_co2 = JutulDarcy.reference_densities(rmodel.system)

    dt = Float64[]
    forces = Dict{Symbol, Any}[]

    w1, w2 = well_labels

    new_period!(time_total, nsteps, forces_for_step) = add_timesteps_and_forces!(dt, forces, time_total, nsteps, forces_for_step)
    # Disable all wells for injection and migration
    no_forces = setup_reservoir_forces(model)
    new_period!(time_initialization, nstep_initialization, no_forces)

    # First injection period
    ctrl1 = Dict{Symbol, Any}()
    ctrl1[w1] = rate_to_injection_control(first(rate_injection1), rho_co2, injection_temperature, enthalpy = H_well)
    ctrl1[w2] = rate_to_injection_control(last(rate_injection1), rho_co2, injection_temperature, enthalpy = H_well)

    forces_injection1 = setup_reservoir_forces(model, control = ctrl1)
    new_period!(time_injection1, nstep_injection1, forces_injection1)

    # Second injection period
    ctrl2 = Dict{Symbol, Any}()
    ctrl2[w1] = rate_to_injection_control(first(rate_injection2), rho_co2, injection_temperature, enthalpy = H_well)
    ctrl2[w2] = rate_to_injection_control(last(rate_injection2), rho_co2, injection_temperature, enthalpy = H_well)

    forces_injection2 = setup_reservoir_forces(model, control = ctrl2)
    new_period!(time_injection2, nstep_injection2, forces_injection2)

    # Finally migrate a bit.
    new_period!(time_migration, nstep_migration, no_forces)

    return (forces, dt)
end

function add_timesteps_and_forces!(dt, forces, time_total, nsteps, forces_for_step)
    if time_total > 0
        for i in 1:nsteps
            push!(dt, time_total/nsteps)
            push!(forces, deepcopy(forces_for_step))
        end
    end
end

function rate_to_injection_control(wrate, co2_density, co2_temp; kwarg...)
    if wrate < 0
        error("Only injection supported")
    elseif wrate == 0
        ctrl = DisabledControl()
    else
        t = TotalRateTarget(wrate/co2_density)
        # CO2 injection
        mix = [0.0, 1.0]
        # Water injection
        # mix = [1.0, 0.0]
        ctrl = InjectorControl(t, mix; density = co2_density, temperature = co2_temp, kwarg...)
    end
    return ctrl
end

function setup_state0_csp11(model, case = :b)
    @assert case == :b
    domain = reservoir_domain(model)
    cc = domain[:cell_centroids]
    wc1 = domain[:well_cells][1]
    z = cc[3, :]
    z0 = cc[3, wc1]
    p0 = 3e7
    # TODO: Improve this
    p = @. p0 + (z - z0)*JutulDarcy.gravity_constant*1000.0
    T0 = @. 70 − 0.025.*(1200.0 - z) + 273.15
    state0 = setup_reservoir_state(model, OverallMoleFractions = [1.0, 0.0], Pressure = p, Temperature = T0)
    return state0
end