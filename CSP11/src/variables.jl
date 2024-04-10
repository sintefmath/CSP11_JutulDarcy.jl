struct CSP11Enthalpy{T} <: VectorVariables
    tab_h::T
    tab_rho_pure::T
    function CSP11Enthalpy(tab::T, tab_rho::T) where T
        new{T}(tab, tab_rho)
    end
end

function Jutul.values_per_entity(model, ::CSP11Enthalpy)
    return 2
end

@jutul_secondary function update_internal_energy!(H_phases, var::CSP11Enthalpy, model, Pressure, Temperature, LiquidMassFractions, VaporMassFractions, PhaseMassDensities, ix)
    fsys = JutulDarcy.flow_system(model.system)
    @assert !JutulDarcy.has_other_phase(fsys)
    @assert JutulDarcy.number_of_components(fsys) == 2
    l, v = JutulDarcy.phase_indices(fsys)

    X, Y = LiquidMassFractions, VaporMassFractions
    rho = PhaseMassDensities
    for c in ix
        p, T = Pressure[c], Temperature[c]
        H_h2o, H_co2 = var.tab_h(p, T)
        rho_h2o_pure, rho_co2_pure = var.tab_rho_pure(p, T)

        # This is heat capacity
        H_h2o -= p/rho_h2o_pure
        H_co2 -= p/rho_co2_pure
        # Should not be needed?
        H_l = X[1, c]*H_h2o + X[2, c]*H_co2
        H_v = Y[1, c]*H_h2o + Y[2, c]*H_co2

        H_phases[l, c] = H_l + p/rho[l, c]
        H_phases[v, c] = H_v + p/rho[v, c]
    end
    return H_phases
end

struct CSP11InternalEnergy{T} <: VectorVariables
    tab_h::T
    tab_rho_pure::T
    function CSP11InternalEnergy(tab::T, tab_rho::T) where T
        new{T}(tab, tab_rho)
    end
end

function Jutul.values_per_entity(model, ::CSP11InternalEnergy)
    return 2
end

@jutul_secondary function update_internal_energy!(C_phases, var::CSP11InternalEnergy, model, Pressure, Temperature, LiquidMassFractions, VaporMassFractions, PhaseMassDensities, ix)
    fsys = JutulDarcy.flow_system(model.system)
    @assert !JutulDarcy.has_other_phase(fsys)
    @assert JutulDarcy.number_of_components(fsys) == 2
    l, v = JutulDarcy.phase_indices(fsys)

    X, Y = LiquidMassFractions, VaporMassFractions
    rho = PhaseMassDensities
    for c in ix
        p, T = Pressure[c], Temperature[c]
        H_h2o, H_co2 = var.tab_h(p, T)
        rho_h2o_pure, rho_co2_pure = var.tab_rho_pure(p, T)

        rho_h2o = PhaseMassDensities[l, c]
        rho_co2 = PhaseMassDensities[v, c]

        # This is heat capacity*T
        C_h2o = H_h2o - p/rho_h2o_pure
        C_co2 = H_co2 - p/rho_co2_pure

        # C_h2o *= rho_h2o_pure/rho_h2o
        # C_co2 *= rho_co2_pure/rho_co2


        C_phases[l, c] = C_h2o
        C_phases[v, c] = C_co2

        # C_phases[l, c] = C_h2o*X[1, c] + C_co2*X[2, c]
        # C_phases[v, c] = C_h2o*Y[1, c] + C_co2*Y[2, c]
    end
    return C_phases
end

struct CSP11InternalEnergy2{T} <: VectorVariables
    tab_h::T
    function CSP11InternalEnergy2(tab::T) where T
        new{T}(tab)
    end
end

function Jutul.values_per_entity(model, ::CSP11InternalEnergy2)
    return 2
end

@jutul_secondary function update_internal_energy!(C_phases, var::CSP11InternalEnergy2, model, Pressure, Temperature, LiquidMassFractions, VaporMassFractions, PhaseMassDensities, ix)
    fsys = JutulDarcy.flow_system(model.system)
    @assert !JutulDarcy.has_other_phase(fsys)
    @assert JutulDarcy.number_of_components(fsys) == 2
    l, v = JutulDarcy.phase_indices(fsys)
    # error()
    X, Y = LiquidMassFractions, VaporMassFractions
    rho = PhaseMassDensities
    for c in ix
        p, T = Pressure[c], Temperature[c]
        E_h2o, E_co2 = var.tab_h(p, T)

        C_phases[l, c] = E_h2o*X[1, c] + E_co2*X[2, c]
        C_phases[v, c] = E_h2o*Y[1, c] + E_co2*Y[2, c]
    end
    # @info "Internal energy" value.(extrema(C_phases))

    return C_phases
end

struct CSP11InternalEnergy3{T} <: VectorVariables
    tab_h::T
    function CSP11InternalEnergy3(tab::T) where T
        new{T}(tab)
    end
end

function Jutul.values_per_entity(model, ::CSP11InternalEnergy3)
    return 2
end

@jutul_secondary function update_internal_energy!(C_phases, var::CSP11InternalEnergy3, model, Pressure, Temperature, LiquidMassFractions, VaporMassFractions, PhaseMassDensities, ix)
    fsys = JutulDarcy.flow_system(model.system)
    @assert !JutulDarcy.has_other_phase(fsys)
    @assert JutulDarcy.number_of_components(fsys) == 2
    l, v = JutulDarcy.phase_indices(fsys)
    # error()
    X, Y = LiquidMassFractions, VaporMassFractions
    rho = PhaseMassDensities
    for c in ix
        p, T = Pressure[c], Temperature[c]
        E_h2o, E_co2 = var.tab_h(p, T)

        C_phases[l, c] = E_h2o
        C_phases[v, c] = E_co2

    end
    # @info "Internal energy" value.(extrema(C_phases))

    return C_phases
end
