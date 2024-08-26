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
    if specase == :b
        pop_1_coords = [4500, 1, 1200-500]
        pop_2_coords = [5100, 1, 1200-1100]
        # POP 1: (4500, 500)
        # POP 2: (5100, 1100)
    elseif specase == :c
        pop_1_coords = [4500, 2500, 1200-655]
        pop_2_coords = [5100, 02500, 1200-1255]
    else
        @error "Only case b and c supported."
    end
    # TODO: Put POPs into data_domain

    # pts = domain.representation.node_points
    pts = domain[:cell_centroids]
    pts = vec(reinterpret(Jutul.StaticArrays.SVector{3, Float64}, pts))
    p1 = find_closest_point(pts, pop_1_coords)
    p2 = find_closest_point(pts, pop_2_coords)

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
