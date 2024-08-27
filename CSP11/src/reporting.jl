struct CSP11ReportHelper{T}
    output_path::String
    P1::Union{Int, Nothing}
    P2::Union{Int, Nothing}
    A::Vector{T}
    B::Vector{T}
    C::Vector{T}
    sealing::Vector{Int}
    boundary::Vector{Int}
    satnum::Vector{Int}
end

function CSP11ReportHelper(p1, p2, A, B, C, sealing, boundary, satnum; path)
    A = reformat_weights(A)
    B = reformat_weights(B)
    C = reformat_weights(C)
    return CSP11ReportHelper(path, p1, p2, A, B, C, sealing, boundary, satnum)
end

mutable struct CSP11ReportHandler
    helper::Union{CSP11ReportHelper{Tuple{Int64, Float64}}, Nothing}
    time::Float64
end

function reformat_weights(weights::Vector{Tuple{Int, T}}) where T<:Real
    return weights
end

function reformat_weights(weights::AbstractVecOrMat{Float64})
    out = Tuple{Int, Float64}[]
    for (i, w) in enumerate(weights)
        if w > 0
            push!(out, (i, w))
        end
    end
    return out
end

function CSP11ReportHelper(domain; path, specase)
    A = domain[:A]
    B = domain[:B]
    C = domain[:C]

    boundary = domain[:boundary]
    obs_pts = domain[:observation_points]
    p1 = findfirst(isequal(1), obs_pts)
    p2 = findfirst(isequal(2), obs_pts)

    satnum = domain[:satnum]
    seal = findall(isequal(1), satnum)
    return CSP11ReportHelper(p1, p2, A, B, C, seal, boundary, satnum, path = path)
end

function get_reporting_hook(pth, domain; specase = :b)
    # Create file
    file_pth = joinpath(pth, "spe11$(specase)_time_series.csv")
    f = open(file_pth, "w")
    println(f, "# t [s], p1 [Pa], p2 [Pa], mobA [kg], immA [kg], dissA [kg], sealA [kg], mobB [kg], immB [kg], dissB [kg], sealB [kg], M_C [m], sealTot [kg], boundTot [kg]")
    close(f)
    @info "Writing to $file_pth"
    #  1. Mobile free phase (CO2 at saturations for which the non-wetting
    # relative permeability exceeds 0); 2. immobile free phase (CO2 at
    # saturations for which the non-wetting relative permeability equals 0); 3.
    # dissolved (CO2 in water phase)
    handler = CSP11ReportHandler(nothing, 0.0)

    function hook(done, report, sim, dt, forces, max_iter, cfg)
        if isnothing(handler.helper)
            @info "Setting up helper..."
            rmodel = reservoir_model(sim.model)
            handler.helper = CSP11ReportHelper(rmodel.data_domain, path = file_pth, specase = :b)
        end
        time_offset = 1000.0*spe11_year

        if report[:success]
            handler.time += report[:dt]
            if handler.time - time_offset >= 0.0
                state = sim.storage.state0.Reservoir
                write_reporting_line!(handler.helper, state, handler.time, time_offset)
            end
        end
        return (done, report)
    end
    return (hook, file_pth, f)
end

function write_reporting_line!(helper::CSP11ReportHelper, state, time, time_offset)
    # fmt = x -> round(x, sigdigits = 16)
    fmt = x -> x

    p_at_1 = state.Pressure[helper.P1]
    p_at_2 = state.Pressure[helper.P2]

    satnum = helper.satnum
    satnum::Vector{Int}
    mobA, immA, dissA, sealA = co2_measures(state, helper.A, satnum)
    mobB, immB, dissB, sealB = co2_measures(state, helper.B, satnum)

    sealTot = sealing_co2(state, helper.sealing)
    boundTot = sealing_co2(state, helper.boundary)
    f = open(helper.output_path, "a")
    println(f, "$(fmt(time-time_offset)), $(fmt(p_at_1)), $(fmt(p_at_2)), $(fmt(mobA)), $(fmt(immA)), $(fmt(dissA)), $(fmt(sealA)), $(fmt(mobB)), $(fmt(immB)), $(fmt(dissB)), $(fmt(sealB)), NaN, $(fmt(sealTot)), $(fmt(boundTot))")
    close(f)
end

function find_closest_point(pts, ref_pt)
    val, p = findmin(
        pt -> sqrt(sum((pt - ref_pt).^2)),
        pts
    )
    return p
end

function co2_measures(state, cells, satnum)
    brine = 1
    gas = 2
    kr = state.RelativePermeabilities
    rho = state.PhaseMassDensities
    X = state.LiquidMassFractions
    Y = state.VaporMassFractions
    s = state.Saturations
    vol = state.FluidVolume

    val_mobile = 0.0
    val_immobile = 0.0
    val_diss = 0.0
    val_seal = 0.0

    for (cell, w) in cells
        @assert 0 < w <= 1.0
        mass_co2_in_gas = vol[cell]*s[gas, cell]*rho[gas, cell]*Y[2, cell]
        mass_co2_in_brine = vol[cell]*s[brine, cell]*rho[brine, cell]*X[2, cell]

        tot_co2_mass = state.TotalMasses[2, cell]

        # eq = isapprox(mass_co2_in_gas + mass_co2_in_brine, tot_co2_mass, rtol = 1e-3, atol = 1e-6)
        # @assert eq
        if satnum[cell] == 1
            val_seal += w*tot_co2_mass
        end
        if kr[gas, cell] > 1e-12
            val_mobile += w*mass_co2_in_gas
        else
            val_immobile += w*mass_co2_in_gas
        end
        val_diss += w*mass_co2_in_brine
    end
    return (val_mobile, val_immobile, val_diss, val_seal)
end

function sealing_co2(state, cells)
    tot = 0.0
    for cell in cells
        tot += state.TotalMasses[2, cell]
    end
    return tot
end

function map_to_reporting_grid(state, weights, total_weights, dims)
    # x [m], z [m], pressure [Pa], gas saturation [-], mass fraction of CO2 in liquid [-],
    # mass fraction of H20 in vapor [-], phase mass density gas [kg/m3],
    # phase mass density water [kg/m3], total mass CO2 [kg] T[C]
    # In file spe11b_spatial_map_<Y>y.csv
    # Give data in reference configuration.
    # Do this in two parts, this part only does the mapping
    # The origin of the coordinate system should be located in the lower-left corner with the x-axis positively
    # oriented towards the right and the z-axis positively oriented towards the top. (The reported x and y
    # values refer to the lower-left corners of each cell in the uniform report grid.) Moreover, note that
    # intensive variables (pressure, saturation, and mass fractions) should be reported as cell-center values,
    # while extensive variables (total mass) should be reported as integral/average values for the cell
    nx, ny, nz = dims
    # Case B is [840, 120]
    # Case C is [168, 100, 120]
    is_case_b = nx == 840 && ny == 120 && nz == 1
    # total_weights tell us how much weight a cell has been assigned in total
    # relative weight of a total mass for a given cell is w_i / total_weights to
    # give the fraction of the total mass assigned to a cell. Sum it up directly
    # from there.
    if is_case_b
        # x, z
        nz = ny
        ny = 1
        # 8400 x 1200 m
    else
        # x, y, z
        # 8400 x 5000 x 1200 m
        @assert dims == [168, 100, 120]
    end
    
    n = prod(dims)
    p = zeros(n)
    sg = zeros(n)
    x_co2 = zeros(n)
    y_h2o = zeros(n)
    rho_w = zeros(n)
    rho_g = zeros(n)
    T = zeros(n)
    co2_mass = zeros(n)

    P = state[:Pressure]
    Temperature = state[:Temperature]
    Sg = state[:Saturations][2, :]
    X_co2 = state[:LiquidMassFractions][2, :]
    Y_h2o = state[:VaporMassFractions][1, :]
    Rho_w = state[:PhaseMassDensities][1, :]
    Rho_g = state[:PhaseMassDensities][2, :]
    CO2_mass = state[:TotalMasses][2, :]

    for (i, j, w) in weights
        p[i] += P[j]*w
        sg[i] += Sg[j]*w
        x_co2[i] += X_co2[j]*w
        y_h2o[i] += Y_h2o[j]*w
        rho_w[i] += Rho_w[j]*w
        rho_g[i] += Rho_g[j]*w
        T[i] += Temperature[j]*w

        co2_mass[i] += CO2_mass[j]*w/total_weights[j]
    end
    return Dict(
        :p => p,
        :sg => sg,
        :x_co2 => x_co2,
        :y_h2o => y_h2o,
        :rho_w => rho_w,
        :rho_g => rho_g,
        :T => T,
        :co2_mass => co2_mass
    )
end

function map_to_reporting_grid(case, states::AbstractVector)
    yr = spe11_year
    # Double check that the case matches the reporting intervals
    @assert case.dt[1] ≈ 1000yr
    @assert length(case.dt) == 201
    for i in 1:200
        @assert case.dt[i+1] ≈ 5yr
    end
    nc = number_of_cells(case.model[:Reservoir].domain)
    rg = case.input_data["G"]["reportingGrid"]
    w = case.input_data["G"]["reportingGrid"]["map"][:, 3]
    I = Int.(case.input_data["G"]["reportingGrid"]["map"][:, 1])
    J = Int.(case.input_data["G"]["reportingGrid"]["map"][:, 2])

    dims = Int.(rg["dims"])
    weights = [(I[i], J[i], w[i]) for i in eachindex(w, I, J)]
    @assert maximum(I) == prod(dims)
    @assert maximum(J) == nc

    total_weights = zeros(nc)
    for (j, w_i) in zip(J, w)
        total_weights[j] += w_i
    end

    return map(x -> map_to_reporting_grid(x, weights, total_weights, dims), states)
end

function write_reporting_grid(case, states, pth, specase::Symbol)
    @assert specase in (:c, :b)
    @assert isdir(pth)
    states = map_to_reporting_grid(case, states)
    @assert length(states) == 201
    # Case B is [840, 120]
    # Case C is [168, 100, 120]

    if specase == :b
        pdims = [8400.0, 1.0, 1200.0]
        nx = 840
        ny = 1
        nz = 120
        # 8400 x 1200 m
        is_3d = false
    else
        # x, y, z
        pdims = [8400.0, 5000.0, 1200.]
        # 8400 x 5000 x 1200 m
        nx, ny, nz = 168, 100, 120
        is_3d = true
    end
    dx = pdims[1]/nx
    dy = pdims[2]/ny
    dz = pdims[3]/nz

    @assert length(states[1][:p]) == nx*ny*nz
    # spe11b_performance_spatial_map_1000y

    for (tno, state) in enumerate(states)
        get_data(k) = reshape(state[k], nx, ny, nz)
        p = get_data(:p)
        sg = get_data(:sg)
        x_co2 = get_data(:x_co2)
        y_h2o = get_data(:y_h2o)
        rho_w = get_data(:rho_w)
        rho_g = get_data(:rho_g)
        T = get_data(:T)
        co2_mass = get_data(:co2_mass)

        nyr = 5*(tno-1)
        file_pth = joinpath(pth, "spe11$(specase)_spatial_map_$(nyr)y.csv")
        @info "Writing step $tno to $file_pth"
        @time write_timestep_dense(file_pth, nx, ny, nz, dx, dy, dz, p, sg, x_co2, y_h2o, rho_g, rho_w, co2_mass, T, is_3d)
    end
end

function write_timestep_dense(file_pth, nx, ny, nz, dx, dy, dz, p, sg, x_co2, y_h2o, rho_g, rho_w, co2_mass, T, is_3d)
    f = open(file_pth, "w")
    if is_3d
        header = "# x [m], y [m], z [m], pressure [Pa], gas saturation [-], mass fraction of CO2 in liquid [-], mass fraction of H20 in vapor [-], phase mass density gas [kg/m3], phase mass density water [kg/m3], total mass CO2 [kg], temperature [C]"
    else
        header = "# x [m], z [m], pressure [Pa], gas saturation [-], mass fraction of CO2 in liquid [-], mass fraction of H20 in vapor [-], phase mass density gas [kg/m3], phase mass density water [kg/m3], total mass CO2 [kg], temperature [C]"
    end
    println(f, header)
    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                x = (i-0.5)*dx
                y = (j-0.5)*dy
                z = (k-0.5)*dz
                # println("$tno: $i $j $k ($nx $ny $nz)")
                if true
                write_dense_line!(
                    f, x, y, z,
                    p[i, j, k],
                    sg[i, j, k], x_co2[i, j, k],
                    y_h2o[i, j, k],
                    rho_g[i, j, k],
                    rho_w[i, j, k],
                    co2_mass[i, j, k],
                    T[i, j, k],
                    is_3d
                )
                end
            end
        end
    end
    close(f)
end

function write_dense_line!(f, x, y, z, p, sg, x_co2, y_h2o, rhog, rhow, mass_co2, T, is_3d)
    # "# x [m], y [m], z [m], pressure [Pa], gas saturation [-], mass fraction of CO2 in liquid [-], mass fraction of H20 in vapor [-], phase mass density gas [kg/m3], phase mass density water [kg/m3], total mass CO2 [kg], temperature [C]"
    function print_entry(x, is_last = false)
        print(f, "$x")
        if !is_last
            print(f, ", ")
        end
    end
    print_entry(x)
    if is_3d
        print_entry(y)
    end
    print_entry(z)
    print_entry(p)
    print_entry(sg)
    print_entry(x_co2)
    print_entry(y_h2o)
    print_entry(rhog)
    print_entry(rhow)
    print_entry(mass_co2)
    print_entry(T)
    print(f, "\n")
end
