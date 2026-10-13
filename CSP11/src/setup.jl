# Exactly 365 days in a year in spec
const spe11_year = 365*si_unit(:day)

function get_path_to_matfile(basename)
    dirname = joinpath(@__DIR__, "..", "..", "data")
    pth = joinpath(dirname, "$basename.mat")
end

"""
    setup_spe11_case(domain_or_dims; case = :b, ...)

Set up SPE11B or SPE11C from either an existing reservoir `DataDomain` or a
tuple of Cartesian dimensions. For B, dimensions are `(nx, nz)` (or `(nx, 1,
nz)`); for C they are `(nx, ny, nz)`.

An existing domain is augmented in-place with benchmark wells, reporting
regions, observation points and boundary conditions. It must provide
`:permeability` and `:porosity`; `:satnum` is inferred from the geometry when
missing.

Set `use_wells=false` to inject CO₂ through reservoir source terms instead of
well and facility models. Rates, injection periods and temperature are retained;
SPE11C rates are distributed along each trajectory by its length in each cell.
"""
function setup_spe11_case(domain_or_dims;
        case = :b,
        name = nothing,
        wells = nothing,
        use_wells::Bool = true,
        input_data = Dict{String, Any}("spe11_case" => String(case)),
        domain_kwargs = NamedTuple(),
        well_kwargs = NamedTuple(),
        thermal = true,
        nstep_initialization = thermal*10,
        nstep_injection1 = 50,
        nstep_injection2 = 50,
        nstep_migration = 100,
        use_reporting_steps = true,
        divide_c_wells = false,
        kgrad = :tpfa,
        kwarg...
    )
    case in (:b, :c) || throw(ArgumentError("Only SPE11 cases :b and :c are supported"))
    if domain_or_dims isa Tuple && all(x -> x isa Integer, domain_or_dims)
        dims = Tuple(Int.(domain_or_dims))
        domain = setup_spe11_domain(dims; case = case, domain_kwargs...)
        source_name = "cartesian_$(join(dims, 'x'))"
    elseif domain_or_dims isa DataDomain
        domain = domain_or_dims
        prepare_spe11_domain!(domain, case)
        source_name = "domain"
    else
        throw(ArgumentError("Expected a reservoir DataDomain or a tuple of Cartesian dimensions"))
    end

    if isnothing(wells)
        wells = setup_spe11_wells!(domain, case; use_wells, divide_c_wells = divide_c_wells, well_kwargs...)
    elseif !haskey(domain, :well_cells)
        throw(ArgumentError("A domain used with custom wells must define domain[:well_cells]"))
    end
    if isnothing(name)
        name = source_name
    end
    return _setup_spe11_case(domain, wells;
        case = case,
        name = name,
        input_data = input_data,
        thermal = thermal,
        use_wells = use_wells,
        nstep_initialization = nstep_initialization,
        nstep_injection1 = nstep_injection1,
        nstep_injection2 = nstep_injection2,
        nstep_migration = nstep_migration,
        use_reporting_steps = use_reporting_steps,
        kgrad = kgrad,
        kwarg...
    )
end

function _setup_spe11_case(domain, wells;
        case,
        name,
        input_data,
        thermal,
        use_wells = true,
        nstep_initialization,
        nstep_injection1,
        nstep_injection2,
        nstep_migration,
        use_reporting_steps,
        kgrad,
        kwarg...
    )

    if thermal
        othername = "thermal"
    else
        othername = "isothermal"
    end
    name = "spe11$(case)_$(name)_$(othername)_$kgrad"
    if !use_wells
        name *= "_sources"
        if !isempty(wells)
            # The MAT path supplies separate point wells and prescribed rates.
            # Also permit custom wells by spreading their rate over perforations.
            cells = [copy(physical_representation(w).perforations.reservoir) for w in wells]
            if !haskey(domain, :spe11_source_cells) || domain[:spe11_source_cells] != cells
                domain[:spe11_source_cells, nothing] = cells
                domain[:spe11_source_weights, nothing] = [fill(1/length(c), length(c)) for c in cells]
            end
        elseif !haskey(domain, :spe11_source_cells)
            throw(ArgumentError("Source mode requires SPE11 injection cells"))
        end
    end

    model, parameters = setup_reservoir_model_csp11(domain;
        wells = use_wells ? wells : [],
        thermal = thermal,
        dT_max_abs = 30.0,
        kgrad = kgrad,
        kwarg...
    )
    check_spe11_volumes(model, parameters, case)
    if case == :b
        forces, dt = CSP11.setup_reservoir_forces_and_timesteps_csp11(model,
            case,
            nstep_initialization = nstep_initialization,
            nstep_injection1 = nstep_injection1,
            nstep_injection2 = nstep_injection2,
            nstep_migration = nstep_migration,
            use_reporting_steps = use_reporting_steps,
            use_wells = use_wells
        );
    elseif case == :c
        rate_injection1 = deepcopy(domain[:well_rates])
        rate_injection2 = deepcopy(rate_injection1)
        num_wells = domain[:num_well_cells]
        rate_injection1[num_wells[1]+1:end] .= 0
        well_labels = []
        for i = 0:sum(num_wells)-1
            push!(well_labels, Symbol(:INJ, i))
        end
        forces, dt = CSP11.setup_reservoir_forces_and_timesteps_csp11(model,
            case,
            well_labels = well_labels,
            nstep_initialization = nstep_initialization,
            nstep_injection1 = nstep_injection1,
            nstep_injection2 = nstep_injection2,
            nstep_migration = nstep_migration,
            rate_injection1 = rate_injection1,
            rate_injection2 = rate_injection2,
            use_reporting_steps = use_reporting_steps,
            use_wells = use_wells
        );
    end

    state0 = setup_state0_csp11(model, case) #check for c

    simulation_case = JutulCase(model, dt, forces;
        state0 = state0,
        parameters = parameters,
        input_data = input_data
    )
    return (simulation_case, name)
end

"""
    setup_spe11_case_from_mrst_grid(basename; kwarg...)

Compatibility wrapper for the original MRST/MAT setup path. New code should
call [`setup_spe11_case`](@ref) with a domain or Cartesian dimensions.
"""
function setup_spe11_case_from_mrst_grid(basename;
        case = :b,
        use_wells::Bool = true,
        thermal = true,
        nstep_initialization = thermal*10,
        nstep_injection1 = 50,
        nstep_injection2 = 50,
        nstep_migration = 100,
        use_reporting_steps = true,
        kgrad = :tpfa,
        kwarg...
    )
    pth = get_path_to_matfile(basename)
    domain, wells, matfile = reservoir_domain_and_wells_csp11(pth, case)
    return _setup_spe11_case(domain, wells;
        case = case,
        name = basename,
        input_data = matfile,
        thermal = thermal,
        use_wells = use_wells,
        nstep_initialization = nstep_initialization,
        nstep_injection1 = nstep_injection1,
        nstep_injection2 = nstep_injection2,
        nstep_migration = nstep_migration,
        use_reporting_steps = use_reporting_steps,
        kgrad = kgrad,
        kwarg...
    )
end

"""
    setup_spe11_domain(dims; case = :b, kwarg...)

Construct a complete Cartesian SPE11 reservoir domain directly from the
published geometry. `dims` is `(nx, nz)` for B and `(nx, ny, nz)` for C. The
case-C mesh is warped into physical space according to the benchmark mapping.
"""
function setup_spe11_domain(dims::Tuple; case = :b, kwarg...)
    case in (:b, :c) || throw(ArgumentError("Only SPE11 cases :b and :c are supported"))
    all(x -> x isa Integer, dims) || throw(ArgumentError("Cartesian dimensions must be integers"))
    dims = Tuple(Int.(dims))
    all(>(0), dims) || throw(ArgumentError("All Cartesian dimensions must be positive"))
    if case == :b
        if length(dims) == 2
            nx, nz = dims
            mesh_dims = (nx, 1, nz)
        elseif length(dims) == 3 && dims[2] == 1
            mesh_dims = dims
        else
            throw(ArgumentError("SPE11B dimensions must be (nx, nz) or (nx, 1, nz)"))
        end
        physical_size = (spe11_size_x, 1.0, spe11_size_z)
    else
        length(dims) == 3 || throw(ArgumentError("SPE11C dimensions must be (nx, ny, nz)"))
        mesh_dims = dims
        physical_size = (spe11_size_x, spe11_size_y, spe11_size_z)
    end

    mesh = UnstructuredMesh(CartesianMesh(mesh_dims, physical_size), z_is_depth = true)
    if case == :c
        for (i, p) in enumerate(mesh.node_points)
            mesh.node_points[i] = SVector(p[1], p[2], p[3] - spe11c_elevation_offset(p[2]))
        end
    end

    cc = tpfv_geometry(mesh).cell_centroids
    satnum = Vector{Int}(undef, size(cc, 2))
    for i in eachindex(satnum)
        satnum[i] = spe11_facies(cc[1, i], cc[2, i], cc[3, i]; case = case)
    end
    permeability, porosity = rock_props_from_satnum(satnum, case)
    if case == :c
        permeability = transform_spe11c_permeability(permeability, cc[2, :])
    end
    # `reservoir_domain` currently checks every compact-tensor entry for
    # non-negativity, although physically valid off-diagonal entries may be
    # negative. Initialize with magnitudes and install the signed tensor after
    # the domain has been constructed.
    permeability_for_constructor = case == :c ? abs.(permeability) : permeability
    domain = reservoir_domain_csp11(mesh, case;
        satnum = satnum,
        permeability = permeability_for_constructor,
        porosity = porosity,
        kwarg...
    )
    if case == :c
        domain[:permeability, Cells()] = permeability
    end
    prepare_spe11_domain!(domain, case)
    return domain
end

function transform_spe11c_permeability(permeability, y)
    nc = length(y)
    size(permeability) == (3, nc) || throw(ArgumentError("Expected a 3×$nc diagonal permeability array"))
    out = zeros(eltype(permeability), 6, nc)
    for i in 1:nc
        k_h = permeability[1, i]
        k_v = permeability[3, i]
        # Derivative of the physical elevation offset in Eq. (4.1). Our third
        # coordinate is depth-positive, hence the negative yz cross term.
        slope = -3/25*((y[i] - 2500.0)/2500.0) + 1/500
        out[:, i] .= (k_h, 0.0, 0.0, k_h, -slope*k_h, k_v + slope^2*k_h)
    end
    return out
end

function _add_spe11_thermal_properties!(domain)
    satnum = domain[:satnum]
    nc = number_of_cells(domain)
    conductivity_by_facies = [1.9, 1.25, 1.25, 1.25, 0.92, 0.26, 2.0]
    if !haskey(domain, :rock_thermal_conductivity)
        domain[:rock_thermal_conductivity, Cells()] = conductivity_by_facies[satnum]
    end
    if !haskey(domain, :diffusion)
        diffusion = repeat([1e-9, 2e-8], 1, nc)
        diffusion[:, satnum .== 7] .= 0.0
        domain[:diffusion] = diffusion
    end
    if !haskey(domain, :fluid_thermal_conductivity)
        domain[:fluid_thermal_conductivity, Cells()] = repeat([0.6, 0.088], 1, nc)
    end
    if !haskey(domain, :rock_density)
        domain[:rock_density, Cells()] = fill(2500.0, nc)
    end
    if !haskey(domain, :component_heat_capacity)
        domain[:component_heat_capacity, Cells()] = repeat([4100.0, 950.0], 1, nc)
    end
    if !haskey(domain, :rock_heat_capacity)
        domain[:rock_heat_capacity, Cells()] = fill(850.0, nc)
    end
    if !haskey(domain, :temperature)
        domain[:temperature, Cells()] = fill(333.15, nc)
    end
    return domain
end

function _spe11_reporting_regions!(domain, case)
    haskey(domain, :A) && haskey(domain, :B) && haskey(domain, :C) && return domain
    cc = domain[:cell_centroids]
    nc = size(cc, 2)
    A = zeros(nc)
    B = zeros(nc)
    C = zeros(nc)
    if case == :b
        boxes = ((3300.0, 8300.0, 0.0, 600.0),
                 (100.0, 3300.0, 600.0, 1200.0),
                 (3300.0, 7800.0, 100.0, 400.0))
    else
        boxes = ((3300.0, 8300.0, 0.0, 750.0),
                 (100.0, 3300.0, 750.0, 1350.0),
                 (3300.0, 7800.0, 250.0, 550.0))
    end
    for i in 1:nc
        x = cc[1, i]
        elevation = spe11_size_z - cc[3, i]
        for (weights, box) in zip((A, B, C), boxes)
            xmin, xmax, zmin, zmax = box
            weights[i] = xmin <= x <= xmax && zmin <= elevation <= zmax
        end
    end
    domain[:A, Cells()] = A
    domain[:B, Cells()] = B
    domain[:C, Cells()] = C
    return domain
end

function _spe11_buffer_volumes(domain, case)
    case in (:b, :c) || throw(ArgumentError("Only SPE11 cases :b and :c are supported"))
    satnum = domain[:satnum]
    buffer = zeros(eltype(domain[:volumes]), length(satnum))
    boundary_centroids = domain[:boundary_centroids]
    isempty(domain[:boundary_neighbors]) && return buffer
    xmin, xmax = extrema(boundary_centroids[1, :])
    ymin, ymax = extrema(boundary_centroids[2, :])
    xtol = max(1.0, xmax - xmin)*1e-8
    ytol = max(1.0, ymax - ymin)*1e-8
    facies = case == :b ? (2:5) : (2:6)
    for face in eachindex(domain[:boundary_neighbors])
        cell = domain[:boundary_neighbors][face]
        x, y = boundary_centroids[1, face], boundary_centroids[2, face]
        on_x_side = abs(x - xmin) <= xtol || abs(x - xmax) <= xtol
        on_c_front_or_back = case == :c && (abs(y - ymin) <= ytol || abs(y - ymax) <= ytol)
        if (on_x_side || on_c_front_or_back) && satnum[cell] in facies
            # Add every eligible face, including both faces of corner cells.
            buffer[cell] += spe11_buffer_length*domain[:boundary_areas][face]
        end
    end
    return buffer
end

function _spe11_boundary_conditions!(domain, case)
    if !haskey(domain, :spe11_physical_volumes)
        if haskey(domain, :boundary)
            # An existing domain may already contain legacy buffer volumes.
            physical = tpfv_geometry(physical_representation(domain)).volumes
        else
            physical = copy(domain[:volumes])
        end
        domain[:spe11_physical_volumes, Cells()] = physical
    end
    buffer = _spe11_buffer_volumes(domain, case)
    domain[:spe11_buffer_volumes, Cells()] = buffer
    domain[:volumes, Cells()] = domain[:spe11_physical_volumes] .+ buffer
    is_boundary = buffer .> 0
    domain[:boundary, nothing] = findall(is_boundary)
    domain[:is_boundary, Cells()] = is_boundary

    thermal_faces = Int[]
    for face in eachindex(domain[:boundary_neighbors])
        normal = domain[:boundary_normals][:, face]
        if abs(normal[3]) > abs(normal[1]) + abs(normal[2])
            push!(thermal_faces, face)
        end
    end
    if haskey(domain, :spe11_fixed_temperature_boundaries) && domain[:spe11_fixed_temperature_boundaries]
        # Migrate domains already prepared by the heat-capacity approximation.
        cells = unique(domain[:boundary_neighbors][thermal_faces])
        domain[:rock_heat_capacity][cells] ./= 1e5
        domain[:spe11_fixed_temperature_boundaries, nothing] = false
    end
    domain[:spe11_thermal_boundary_faces, nothing] = thermal_faces
    return domain
end

_spe11_geothermal_temperature(depth) = 70.0 - 0.025*(1200.0 - depth) + 273.15

function _spe11_thermal_boundary_conditions(model)
    rmodel = JutulDarcy.reservoir_model(model)
    JutulDarcy.model_is_thermal(rmodel) || return nothing
    domain = reservoir_domain(model)
    phi = domain[:porosity]
    # The impermeable top/bottom facies contain water. Match the reservoir's
    # parallel rock/fluid conductivity, using the physical boundary geometry.
    fluid = domain[:fluid_thermal_conductivity]
    water = fluid isa AbstractVector ? fluid : vec(fluid[1, :])
    conductivity = (1.0 .- phi).*domain[:rock_thermal_conductivity] .+ phi.*water
    trans = compute_boundary_trans(domain, conductivity)
    # Aggregate by cell so even a one-layer grid has one force per cell. This
    # avoids concurrent writes when forces are evaluated by device kernels.
    conductance = Dict{Int, Float64}()
    weighted_temperature = Dict{Int, Float64}()
    for face in domain[:spe11_thermal_boundary_faces]
        cell = domain[:boundary_neighbors][face]
        t = trans[face]
        temperature = _spe11_geothermal_temperature(domain[:boundary_centroids][3, face])
        conductance[cell] = get(conductance, cell, 0.0) + t
        weighted_temperature[cell] = get(weighted_temperature, cell, 0.0) + t*temperature
    end
    return [FlowBoundaryCondition(cell, 3e7, weighted_temperature[cell]/conductance[cell];
        trans_flow = 0.0, trans_thermal = conductance[cell])
        for cell in sort!(collect(keys(conductance))) if conductance[cell] > 0.0]
end

"""
    check_spe11_volumes(model, parameters, case)

Check physical geometry, lateral buffer storage and initialized model pore
volume. Returns measured and integrated reference volumes in cubic metres.
Facies discretization differences are logged at debug level; inconsistent
geometry or model storage produces warnings with actual and expected values.
"""
function check_spe11_volumes(model, parameters, case)
    domain = reservoir_domain(model)
    reference = spe11_reference_volumes(case)
    physical = sum(domain[:spe11_physical_volumes])
    expected_bulk = physical + sum(_spe11_buffer_volumes(domain, case))
    bulk = sum(domain[:volumes])
    expected_pore = sum(pore_volume(domain))
    pore = sum(pore_volume(model, parameters))
    function check(message, actual, expected; rtol = 1e-8)
        relative_error = abs(actual - expected)/expected
        if !isfinite(relative_error) || relative_error > rtol
            @warn message actual expected relative_error
        end
    end
    check("SPE11 physical bulk volume differs from the specified domain dimensions",
        physical, reference.physical_bulk)
    check("SPE11 effective bulk volume is inconsistent with the mesh and lateral buffer areas",
        bulk, expected_bulk)
    check("SPE11 model pore volume is inconsistent with the reservoir domain",
        pore, expected_pore)
    bulk_error = abs(bulk - reference.bulk)/reference.bulk
    pore_error = abs(pore - reference.pore)/reference.pore
    if max(bulk_error, pore_error) > 1e-3
        # Sampling facies at cell centres changes their volumes on coarse
        # meshes. Report this discretization error without a setup warning.
        @debug "SPE11 discretized storage volumes differ from the integrated facies geometry" bulk pore reference_bulk = reference.bulk reference_pore = reference.pore bulk_relative_error = bulk_error pore_relative_error = pore_error
    end
    return (physical = physical, bulk = bulk, pore = pore, reference = reference,
        bulk_relative_error = bulk_error, pore_relative_error = pore_error)
end

function _spe11_observation_points!(domain, case)
    haskey(domain, :observation_points) && return domain
    if case == :b
        points = ([4500.0, 0.5, 700.0], [5100.0, 0.5, 100.0])
    else
        points = ([4500.0, 2500.0, 545.0], [5100.0, 2500.0, -55.0])
    end
    cc = domain[:cell_centroids]
    cell_points = vec(reinterpret(SVector{3, Float64}, cc))
    observations = zeros(Int, number_of_cells(domain))
    for (i, point) in enumerate(points)
        observations[find_closest_point(cell_points, point)] = i
    end
    domain[:observation_points, Cells()] = observations
    return domain
end

function prepare_spe11_domain!(domain::DataDomain, case)
    case in (:b, :c) || throw(ArgumentError("Only SPE11 cases :b and :c are supported"))
    for key in (:permeability, :porosity)
        haskey(domain, key) || throw(ArgumentError("The supplied domain must define :$key"))
    end
    if !haskey(domain, :satnum)
        cc = domain[:cell_centroids]
        satnum = [spe11_facies(cc[1, i], cc[2, i], cc[3, i]; case = case) for i in axes(cc, 2)]
        domain[:satnum, Cells()] = satnum
    end
    domain[:porosity][domain[:satnum] .== 7] .= spe11_facies7_porosity
    _add_spe11_thermal_properties!(domain)
    _spe11_reporting_regions!(domain, case)
    _spe11_boundary_conditions!(domain, case)
    _spe11_observation_points!(domain, case)
    return domain
end

function _spe11_trajectory_cells(domain, trajectory; n = 501)
    mesh = physical_representation(domain) |> UnstructuredMesh
    cells, extra = Jutul.find_enclosing_cells(mesh, trajectory; n = n, extra_out = true)
    isempty(cells) && error("The SPE11 well trajectory did not intersect the supplied domain")
    lengths = extra[:lengths]
    return cells, lengths
end

function setup_spe11_wells!(domain::DataDomain, case; use_wells = true, divide_c_wells = false, kwarg...)
    default_options = (simple_well = true, radius = 0.15, dir = :y)
    options = merge(default_options, values(kwarg))
    if case == :b
        cc = domain[:cell_centroids]
        points = vec(reinterpret(SVector{3, Float64}, cc))
        well_points = ([2700.0, 0.5, 900.0], [5100.0, 0.5, 500.0])
        well_cells = [find_closest_point(points, p) for p in well_points]
        wells = use_wells ? [setup_well(domain, well_cells[i]; options..., name = Symbol(:INJ, i-1)) for i in 1:2] : []
        domain[:well_cells, nothing] = well_cells
        source_cells = [[c] for c in well_cells]
        source_weights = [[1.0], [1.0]]
    else
        trajectory1 = [2700.0 1000.0 900.0; 2700.0 4000.0 900.0]
        y = collect(range(1000.0, 4000.0, length = 101))
        trajectory2 = hcat(fill(5100.0, length(y)), y, 500.0 .- spe11c_elevation_offset.(y))
        cells1, lengths1 = _spe11_trajectory_cells(domain, trajectory1; n = 1001)
        cells2, lengths2 = _spe11_trajectory_cells(domain, trajectory2; n = 11)
        cells = [cells1; cells2]
        if divide_c_wells
            wells = Any[]
            if use_wells
                for i in eachindex(cells)
                    push!(wells, setup_well(domain, cells[i]; options..., name = Symbol(:INJ, i-1)))
                end
            end
            source_cells = [[c] for c in cells]
            source_weights = [[1.0] for c in cells]
            rates1 = 50.0.*lengths1./sum(lengths1)
            rates2 = 50.0.*lengths2./sum(lengths2)
            n1 = length(cells1)
            n2 = length(cells2)
        else
            wells = use_wells ? [setup_well(domain, cells1; options..., name = :INJ0),
                setup_well(domain, cells2; options..., name = :INJ1)] : []
            source_cells = [cells1, cells2]
            source_weights = [lengths1./sum(lengths1), lengths2./sum(lengths2)]
            rates1 = [50.0]
            rates2 = [50.0]
            n1 = n2 = 1
        end
        domain[:well_cells, nothing] = cells
        domain[:num_well_cells, nothing] = [n1, n2]
        domain[:well_rates, nothing] = [rates1; rates2].*si_unit(:kilogram)./si_unit(:second)
    end
    domain[:spe11_source_cells, nothing] = source_cells
    domain[:spe11_source_weights, nothing] = source_weights
    return wells
end

function reservoir_domain_and_wells_csp11(pth::AbstractString, case = :b; kwarg...)
    matdata = MAT.matread(pth)
    raw_rock = missing
    if haskey(matdata, "G")
        raw_G = matdata["G"]
        if haskey(matdata, "rock")
            raw_rock = matdata["rock"]
        end
    else
        @assert haskey(matdata, "cells")
        raw_G = matdata
    end
    G = UnstructuredMesh(MRSTWrapMesh(raw_G), z_is_depth = true)
    satnum = Int.(vec(raw_G["cells"]["tag"]))
    # TODO: Special perm transform for case C
    if ismissing(raw_rock)
        @assert case == :b
        K, poro = rock_props_from_satnum(satnum, case)
    else
        @info "Perm from rock"
        K = collect(raw_rock["perm"]')
        if case == :c
            @assert size(K, 1) == 6
        end
        @. K = max(K, 1e-10*si_unit(:darcy))
        poro = collect(vec(raw_rock["poro"]))
        poro[satnum .== 7] .= spe11_facies7_porosity
    end
    domain = reservoir_domain_csp11(G, case; satnum = satnum, permeability = K, porosity = poro, kwarg...)
    _spe11_boundary_conditions!(domain, case)

    simple_well = true
    if case == :b
        wc1, wc2 = Int.(vec(raw_G["cells"]["wellCells"]))
        I0 = setup_well(domain, [wc1], simple_well = simple_well, name = :INJ0)
        I1 = setup_well(domain, [wc2], simple_well = simple_well, name = :INJ1)
        wells = [I0, I1]
        domain[:well_cells, nothing] = [wc1, wc2]
        pop_1_coords = [4500, 1, 1200-500]
        pop_2_coords = [5100, 1, 1200-1100]
    elseif case == :c
        w = raw_G["cells"]["wellCells"]
        wc1, wc2 = Int.(vec(w[1])), Int.(vec(w[2]))
        wells1 = Any[]
        wells2 = Any[]
        i=0
        for iw1 in eachindex(wc1)
            push!(wells1, setup_well(domain, wc1[iw1], simple_well = simple_well, name = Symbol(:INJ, i)))
            i += 1
        end
        for iw2 in eachindex(wc2)
            push!(wells2, setup_well(domain, wc2[iw2], simple_well = simple_well, name = Symbol(:INJ, i)))
            i += 1
        end
        wells = [wells1; wells2]
        domain[:well_cells, nothing] = [wc1; wc2]
        domain[:num_well_cells, nothing] = [length(wc1), length(wc2)]
        rates_well_1 = vec(raw_G["cells"]["wellMassRate"][1])
        rates_well_2 = vec(raw_G["cells"]["wellMassRate"][2])
        domain[:well_rates, nothing] = [rates_well_1; rates_well_2].*si_unit(:kilogram)./si_unit(:second)
        pop_1_coords = [4500, 2500, 1200-655]
        pop_2_coords = [5100, 02500, 1200-1255]
    end

    A = raw_G["cells"]["fractionInA"]
    B = raw_G["cells"]["fractionInB"]
    C = raw_G["cells"]["fractionInC"]

    domain[:A, Cells()] = vec(A)
    domain[:B, Cells()] = vec(B)
    domain[:C, Cells()] = vec(C)

    # Obervation points
    nc = number_of_cells(domain)
    pts = domain[:cell_centroids]
    pts = vec(reinterpret(Jutul.StaticArrays.SVector{3, Float64}, pts))
    p1 = find_closest_point(pts, pop_1_coords)
    p2 = find_closest_point(pts, pop_2_coords)

    observation_points = zeros(Int, nc)
    observation_points[p1] = 1
    observation_points[p2] = 2
    domain[:observation_points, Cells()] = observation_points

    # domain[:well_cells, nothing] = [wc1, wc2]
    return domain, wells, matdata
end

function reservoir_domain_csp11(G, case = :b; satnum, temperature = 333.15,
        porosity = last(rock_props_from_satnum(satnum, case)), kwarg...)
    maximum(satnum) == 7 || throw(ArgumentError("Must have 7 as highest SATNUM region"))
    minimum(satnum) == 1 || throw(ArgumentError("Must have 1 as lowest SATNUM region"))
    nc = number_of_cells(G)
    length(satnum) == nc || throw(ArgumentError("satnum must have number of cells entries ($nc), was $(length(satnum))"))
    porosity = porosity isa Real ? fill(porosity, nc) : copy(porosity)
    porosity[satnum .== 7] .= spe11_facies7_porosity
    domain = reservoir_domain(G; satnum = satnum, porosity = porosity, kwarg...)
    # TODO: Move over perm and poro assignment here.
    if case == :b || case == :c
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
        throw(ArgumentError("Only case b and c(?) is supported at the moment."))
    end
    domain[:temperature, Cells()] = temperature
    return domain
end

function rock_props_from_satnum(satnum, case)
    n = length(satnum)
    perm = zeros(3, n)
    poro = zeros(n)
    @assert minimum(satnum) >= 1
    @assert maximum(satnum) <= 7
    MINPERM = 1e-10*si_unit(:darcy)
    MINPORO = spe11_facies7_porosity
    if case == :a
        z_factor = 1.0
        # Note: Last entry was 0, capped to very small perm and small poro
        perm_reg = [4e-11, 5e-10, 1e-9, 2e-9, 4e-9, 1e-8, MINPERM]
        poro_reg = [0.44, 0.43, 0.44, 0.45, 0.43, 0.46, MINPORO]
    else
        z_factor = 0.1
        perm_reg = [1e-16, 1e-13, 2e-13, 5e-13, 1e-12, 2e-12, MINPERM]
        poro_reg = [0.1, 0.2, 0.2, 0.2, 0.25, 0.35, MINPORO]
    end
    @assert length(perm_reg) == length(poro_reg) == 7
    for (i, reg) in enumerate(satnum)
        perm[1, i] = perm[2, i] = perm_reg[reg]
        perm[3, i] = z_factor*perm_reg[reg]
        poro[i] = poro_reg[reg]
    end
    return (perm, poro)
end

const spe11_wetting_immobile_saturation = (0.32, 0.14, 0.12, 0.12, 0.12, 0.10)
const spe11_nonwetting_immobile_saturation = 0.10
const spe11_relative_permeability_exponent = 1.5
const spe11_capillary_pressure_exponent = 1.5
const spe11_max_capillary_pressure = 3.0e7
const spe11_leverett_coefficient = 6.12e-3

function spe11_phase_relative_permeability(immobile, label; n = 65)
    effective_saturation = range(0.0, 1.0; length = n)
    saturation = immobile .+ (1.0 - immobile).*effective_saturation
    kr = effective_saturation.^spe11_relative_permeability_exponent
    return PhaseRelativePermeability(saturation, kr; label = label)
end

function spe11_erf(x)
    ax = abs(x)
    if ax < 0.5
        # The power series is accurate and well-conditioned near zero.
        term = total = ax
        for n in 1:100
            term *= -ax^2/n
            increment = term/(2n + 1)
            total += increment
            abs(increment) <= eps(Float64)*abs(total) && break
        end
        value = 2.0/sqrt(pi)*total
    else
        # Compact approximation with a maximum absolute error around 3e-8.
        t = 1.0/(1.0 + 0.5ax)
        polynomial = 0.17087277
        polynomial = -0.82215223 + t*polynomial
        polynomial = 1.48851587 + t*polynomial
        polynomial = -1.13520398 + t*polynomial
        polynomial = 0.27886807 + t*polynomial
        polynomial = -0.18628806 + t*polynomial
        polynomial = 0.09678418 + t*polynomial
        polynomial = 0.37409196 + t*polynomial
        polynomial = 1.00002368 + t*polynomial
        tail = t*exp(-ax^2 - 1.26551223 + t*polynomial)
        value = 1.0 - tail
    end
    return signbit(x) ? -value : value
end

function spe11_capillary_pressure(sg, immobile, entry_pressure)
    sw = 1.0 - sg
    normalized_sw = max((sw - immobile)/(1.0 - immobile), 0.0)
    if iszero(normalized_sw)
        return spe11_max_capillary_pressure
    end
    pc_unbounded = entry_pressure*normalized_sw^(-1.0/spe11_capillary_pressure_exponent)
    erf_argument = pc_unbounded/spe11_max_capillary_pressure*sqrt(pi)/2.0
    return spe11_max_capillary_pressure*spe11_erf(erf_argument)
end

function spe11_capillary_pressure_table(immobile, entry_pressure;
        relative_tolerance = 1e-3,
        absolute_tolerance = 1.0
    )
    pc(sg) = spe11_capillary_pressure(sg, immobile, entry_pressure)
    mobile_limit = 1.0 - immobile
    effective_saturation = vcat(
        collect(range(0.0, 1.0; length = 17)),
        10.0.^range(-14.0, 0.0; length = 33)
    )
    sort!(effective_saturation)
    unique!(effective_saturation)
    points = 1.0 .- (immobile .+ (1.0 - immobile).*effective_saturation)
    push!(points, 1.0)
    sort!(points)
    unique!(points)

    function refine_interval!(lo, hi, pc_lo, pc_hi, level = 0)
        level == 40 && return
        midpoint = (lo + hi)/2.0
        pc_midpoint = pc(midpoint)
        needs_refinement = any((0.25, 0.5, 0.75)) do fraction
            sg = lo + fraction*(hi - lo)
            pc_exact = pc(sg)
            pc_linear = pc_lo + fraction*(pc_hi - pc_lo)
            tolerance = max(absolute_tolerance, relative_tolerance*abs(pc_exact))
            abs(pc_exact - pc_linear) > tolerance
        end
        if needs_refinement
            push!(points, midpoint)
            refine_interval!(lo, midpoint, pc_lo, pc_midpoint, level + 1)
            refine_interval!(midpoint, hi, pc_midpoint, pc_hi, level + 1)
        end
        return
    end

    initial_points = copy(points)
    for i in 1:(length(initial_points) - 1)
        lo = initial_points[i]
        hi = initial_points[i + 1]
        refine_interval!(lo, hi, pc(lo), pc(hi))
    end
    sort!(points)
    pressure = pc.(points)
    return Jutul.LinearInterpolant(points, pressure; constant_dx = false)
end

"""
    spe11_saturation_functions(satnum)

Construct the SPE11B/C Brooks-Corey relative-permeability and capillary-
pressure functions directly from the benchmark parameters. The original
PYOPMSPE11 deck represented the same functions with 200,000 coupled `SGOF`
rows per permeable facies. Relative permeability is sampled independently on a
small uniform effective-saturation grid, while capillary pressure uses adaptive
sampling around its much sharper transition.
"""
function spe11_saturation_functions(satnum)
    regions = Int.(vec(satnum))
    all(in(1:7), regions) || throw(ArgumentError("SPE11 saturation regions must be in 1:7"))

    gas = spe11_phase_relative_permeability(
        spe11_nonwetting_immobile_saturation, :g)
    impermeable_gas = PhaseRelativePermeability([0.0, 1.0], [0.0, 1.0]; label = :g)
    gas_tables = ntuple(i -> i == 7 ? impermeable_gas : gas, 7)

    impermeable_liquid = PhaseRelativePermeability([0.0, 1.0], [0.0, 1.0]; label = :og)
    liquid_tables = ntuple(7) do i
        if i == 7
            impermeable_liquid
        else
            spe11_phase_relative_permeability(spe11_wetting_immobile_saturation[i], :og)
        end
    end
    relative_permeabilities = JutulDarcy.ReservoirRelativePermeabilities(;
        og = liquid_tables,
        g = gas_tables,
        regions = regions
    )

    porosity = (0.10, 0.20, 0.20, 0.20, 0.25, 0.35)
    horizontal_permeability = (1e-16, 1e-13, 2e-13, 5e-13, 1e-12, 2e-12)
    capillary_tables = ntuple(7) do i
        if i == 7
            # The deck used a two-point fallback for the impermeable facies.
            Jutul.LinearInterpolant(
                [0.0, 1.0], [0.0, spe11_max_capillary_pressure];
                constant_dx = false
            )
        else
            entry_pressure = spe11_leverett_coefficient*sqrt(
                porosity[i]/horizontal_permeability[i])
            spe11_capillary_pressure_table(
                spe11_wetting_immobile_saturation[i], entry_pressure)
        end
    end
    capillary_pressure = JutulDarcy.SimpleCapillaryPressure(
        (capillary_tables, ); regions = regions)
    return relative_permeabilities, capillary_pressure
end

function setup_reservoir_model_csp11(reservoir::DataDomain; include_satfun = true, kwarg...)
    model, parameters = setup_reservoir_model(reservoir, :co2brine;
        co2_source = :csp11,
        extra_out = true,
        kwarg...
    )
    if include_satfun
        kr, pc = spe11_saturation_functions(reservoir[:satnum])
        set_secondary_variables!(model[:Reservoir],
            RelativePermeabilities = kr,
            CapillaryPressure = pc
        )
    end
    return (model, parameters)
end

function setup_reservoir_forces_and_timesteps_csp11(model, case = :b;
        use_wells = true,
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
        use_reporting_steps = false,
        injection_temperature = convert_to_si(10, :Celsius)
    )
    if !(case == :b || case == :c)
        throw(ArgumentError("Only case b and c supported at the moment."))
    end

    if use_reporting_steps
        @info "Using sequence of time-steps required for reporting... Overriding any kwarg set."
        # 5 year intervals, 25 years
        nstep_initialization = 1
        nstep_injection1 = 5
        nstep_injection2 = 5
        # 1000 - 2*25 = 950 / 5 = 190
        nstep_migration = 190
    end

    rmodel = JutulDarcy.reservoir_model(model)
    rho_brine, rho_co2 = JutulDarcy.reference_densities(rmodel.system)

    dt = Float64[]
    forces = Dict{Symbol, Any}[]

    new_period!(time_total, nsteps, forces_for_step) = add_timesteps_and_forces!(dt, forces, time_total, nsteps, forces_for_step)
    tables = JutulDarcy.CO2Properties.co2_brine_property_tables()
    H_tab = tables[:enthalpy]
    co2_H_eval(p, T) = H_tab(p, T)[2]
    H_well = co2_H_eval

    if !use_wells
        no_forces = setup_reservoir_forces(model)
        new_period!(time_initialization, nstep_initialization, no_forces)
        forces_injection1 = spe11_source_forces(model, rate_injection1, injection_temperature)
        new_period!(time_injection1, nstep_injection1, forces_injection1)
        forces_injection2 = spe11_source_forces(model, rate_injection2, injection_temperature)
        new_period!(time_injection2, nstep_injection2, forces_injection2)
        new_period!(time_migration, nstep_migration, no_forces)
    elseif case == :b
        w1, w2 = well_labels

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
    elseif case == :c
        # Disable all wells for injection and migration
        no_forces = setup_reservoir_forces(model)
        new_period!(time_initialization, nstep_initialization, no_forces)

        # First injection period
        ctrl1 = Dict{Symbol, Any}()
        for (i, name) in enumerate(well_labels)
            ctrl1[name] = rate_to_injection_control(rate_injection1[i], rho_co2, injection_temperature)
        end
        forces_injection1 = setup_reservoir_forces(model, control = ctrl1)
        new_period!(time_injection1, nstep_injection1, forces_injection1)

        # Second injection period
        ctrl2 = Dict{Symbol, Any}()
        for (i, name) in enumerate(well_labels)
            ctrl2[name] = rate_to_injection_control(rate_injection2[i], rho_co2, injection_temperature)
        end
        
        forces_injection2 = setup_reservoir_forces(model, control = ctrl2)
        new_period!(time_injection2, nstep_injection2, forces_injection2)

        # Finally migrate a bit.
        new_period!(time_migration, nstep_migration, no_forces)
    end
    if use_reporting_steps
        @assert length(dt) == 201
    end

    bc = _spe11_thermal_boundary_conditions(model)
    if !isnothing(bc)
        for force in forces
            force[:Reservoir] = merge(force[:Reservoir], (bc = bc,))
        end
    end

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
    @assert case == :b || case == :c
    domain = reservoir_domain(model)
    cc = domain[:cell_centroids]
    wc1 = domain[:well_cells][1]
    z = cc[3, :]
    z0 = cc[3, wc1]
    p0 = 3e7
    # TODO: Improve this
    p = @. p0 + (z - z0)*JutulDarcy.gravity_constant*1000.0
    T0 = _spe11_geothermal_temperature.(z)
    state0 = setup_reservoir_state(model, OverallMoleFractions = [1.0, 0.0], Pressure = p, Temperature = T0)
    return state0
end
