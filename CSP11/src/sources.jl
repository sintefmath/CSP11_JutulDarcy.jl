"""Prescribed pure-CO₂ mass injection, with its stream temperature in kelvin."""
struct SPE11Source <: Jutul.JutulForce
    cell::Int
    value::Float64
    temperature::Float64
end

function spe11_source_forces(model, rates, temperature)
    domain = reservoir_domain(model)
    cells = domain[:spe11_source_cells]
    weights = domain[:spe11_source_weights]
    length(cells) == length(weights) == length(rates) ||
        throw(ArgumentError("SPE11 source cells, weights and injection rates must match"))
    # Trajectories can share a cell on coarse meshes. Aggregate before launching
    # cell kernels so each residual entry has just one writer.
    cell_rates = Dict{Int, Float64}()
    for (well_cells, well_weights, rate) in zip(cells, weights, rates)
        rate >= 0 || throw(ArgumentError("SPE11 sources only support injection"))
        length(well_cells) == length(well_weights) || throw(DimensionMismatch())
        isapprox(sum(well_weights), 1.0) || throw(ArgumentError("SPE11 source weights must sum to one"))
        for (cell, weight) in zip(well_cells, well_weights)
            cell_rates[cell] = get(cell_rates, cell, 0.0) + rate*weight
        end
    end
    sources = [SPE11Source(cell, cell_rates[cell], temperature)
        for cell in sort!(collect(keys(cell_rates))) if cell_rates[cell] > 0]
    return setup_reservoir_forces(model; sources)
end

function Jutul.apply_forces_to_equation!(acc, storage, model::SimulationModel,
        eq::ConservationLaw{:TotalMasses}, eq_s, sources::AbstractVector{SPE11Source}, time)
    Jutul.threaded_loop(length(sources), model.context) do i
        @inbounds source = sources[i]
        # SPE11 component ordering is H₂O, CO₂, with rates prescribed in kg/s.
        @inbounds acc[2, source.cell] -= source.value
    end
    return acc
end

spe11_source_heat_capacity(variable, state, cell, pressure, temperature) =
    state.ComponentHeatCapacity[2, cell]

function spe11_source_heat_capacity(variable::JutulDarcy.PressureTemperatureDependentVariable,
        state, cell, pressure, temperature)
    region = JutulDarcy.region(variable.regions, cell)
    table = JutulDarcy.table_by_region(variable.tab, region)
    return table(pressure, temperature)[2]
end

function Jutul.apply_forces_to_equation!(acc, storage, model::SimulationModel,
        eq::ConservationLaw{:TotalThermalEnergy}, eq_s, sources::AbstractVector{SPE11Source}, time)
    state = Jutul.evaluation_state(storage)
    density = model.secondary_variables[:PhaseMassDensities]
    capacity = get(model.secondary_variables, :ComponentHeatCapacity, nothing)
    map = Jutul.global_map(model.domain)
    Jutul.threaded_loop(length(sources), model.context) do i
        @inbounds source = sources[i]
        cell = Jutul.full_cell(source.cell, map)
        @inbounds pressure = state.Pressure[cell]
        rho_co2 = density.tab(pressure, source.temperature)[2]
        cv = spe11_source_heat_capacity(capacity, state, cell, pressure, source.temperature)
        # Match the well path's constant-volume thermal model for pure CO₂:
        # h = u + p/rho, evaluated at injection temperature and local pressure.
        enthalpy = cv*source.temperature + pressure/rho_co2
        @inbounds acc[1, source.cell] -= source.value*enthalpy
    end
    return acc
end

function Jutul.subforce(sources::AbstractVector{SPE11Source}, model)
    map = Jutul.global_map(model.domain)
    local_sources = SPE11Source[]
    for source in sources
        Jutul.global_cell_inside_domain(source.cell, map) || continue
        cell = Jutul.interior_cell(Jutul.local_cell(source.cell, map), map)
        isnothing(cell) && continue
        push!(local_sources, SPE11Source(cell, source.value, source.temperature))
    end
    return local_sources
end
