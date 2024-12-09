using CSV, DataFrames
using Jutul

specase = :b
cases, = get_RSC_basenames(grids=[:C, :HC], resolutions=["10k"])
name1 = cases[1]
name2 = cases[2]


# specase = :c
# name1 = "horz_ndg_cut_PG_50x50x50_thermal_cv"
# name2 = "struct50x50x50_thermal_cv"

output_folder = jutul_output_path("spe11$(specase)_$(name1)_thermal_cv_tpfa", subfolder = "csp11_rsc")
output_folder2 = jutul_output_path("spe11$(specase)_$(name2)_thermal_cv_tpfa", subfolder = "csp11_rsc")

file_name = "spe11$(specase)_time_series.csv"

function resample_table(t, vals)
    t_rep = 3.1536e6 # 0.1 year
    start = 0.0
    stop = t[end]
    t_i = collect(range(start, stop, step = t_rep))

    new_vals = []
    for v in vals
        I = get_1d_interpolator(t, v)
        push!(new_vals, I.(t_i))
    end

    return (t_i, new_vals)
end

function read_file(pth; resample = false)
    df = CSV.read(pth, DataFrame)

    time = df[:, 1]

    p1 = df[:, 2]
    p2 = df[:, 3]

    mobA = df[:, 4]
    immA = df[:, 5]
    dissA = df[:, 6]
    sealA = df[:, 7]

    mobB = df[:, 8]
    immB = df[:, 9]
    dissB = df[:, 10]
    sealB = df[:, 11]
    M = df[:, 12]
    sealTot = df[:, 13]
    boundTot = df[:, 14]

    if resample
        time, new_tab = resample_table(time, (p1, p2, mobA, immA, dissA, sealA, mobB, immB, dissB, sealB, M, sealTot, boundTot))
        p1, p2, mobA, immA, dissA, sealA, mobB, immB, dissB, sealB, M, sealTot, boundTot = new_tab
    end
    uyear = 365.0*3600*24.0
    return (t = time/uyear,
            p1 = p1, p2 = p2,
                A = (mob = mobA, imm=immA, diss=dissA, seal=sealA),
                B = (mob = mobB, imm=immB, diss=dissB, seal=sealB),
            M = M, sealTot = sealTot, boundTot = boundTot)
end


tab = read_file(joinpath(output_folder, file_name))
tab_name = "C"

tab2 = read_file(joinpath(output_folder2, file_name))
tab2_name = "HC"
##
# using CairoMakie
using GLMakie

SCATTER_ARG = (markersize = 4,)

function plot_co2!(ax, tab, reg = :A; use_scatter = false, name = "", kwarg...)
    if isnothing(tab)
        return tab
    end
    function genlabel(x)
        if name == ""
            return x
        else
            return "$x ($name)"
        end
    end
    if use_scatter
        f! = scatter!
        ARG = SCATTER_ARG
    else
        f! = lines!
        ARG = NamedTuple()
    end
    t = tab.t
    X = tab[reg]
    f!(ax, t, X.mob, label = genlabel("mobile"); ARG...)
    f!(ax, t, X.imm, label = genlabel("immobile"); ARG...)
    f!(ax, t, X.diss, label = genlabel("dissolved"); ARG...)
    f!(ax, t, X.seal, label = genlabel("in sealing"); ARG...)
    f!(ax, t, X.mob + X.imm + X.diss, label = genlabel("total"); ARG...)
end

function plot_pressure!(ax, tab; use_scatter = false, name = "", kwarg...)
    if isnothing(tab)
        return tab
    end
    if use_scatter
        f! = scatter!
        ARG = SCATTER_ARG
    else
        f! = lines!
        ARG = NamedTuple()
    end
    function genlabel(x)
        if name == ""
            return x
        else
            return "$x ($name)"
        end
    end
    t = tab.t
    f!(t, tab.p1/1e6, label = genlabel("P1"); ARG...)
    f!(t, tab.p2/1e6, label = genlabel("P2"); ARG...)
end

function plot_other!(ax, tab; use_scatter = false, name = "", kwarg...)
    if isnothing(tab)
        return tab
    end
    if use_scatter
        f! = scatter!
        ARG = SCATTER_ARG
    else
        f! = lines!
        ARG = NamedTuple()
    end
    function genlabel(x)
        if name == ""
            return x
        else
            return "$x ($name)"
        end
    end
    t = tab.t
    f!(t, tab.sealTot, label = genlabel("seal"); ARG...)
    f!(t, tab.boundTot, label = genlabel("bound"); ARG...)
end

fig = Figure(size = (1200, 800))

ax = Axis(fig[1, 1], title = "CO₂ in region A", xlabel = "Time / year", ylabel = "Mass / kg")
plot_co2!(ax, tab, :A, name = tab_name)
plot_co2!(ax, tab2, :A, name = tab2_name, use_scatter = true)

axislegend(position = :rc)

ax = Axis(fig[1, 2], title = "CO₂ in region B", xlabel = "Time / year", ylabel = "Mass / kg")
plot_co2!(ax, tab, :B, name = tab_name)
plot_co2!(ax, tab2, :B, name = tab2_name, use_scatter = true)

axislegend(position = :rc)

ax = Axis(fig[2, 1], title = "Pressure at observation points", xlabel = "Time / year", ylabel = "Pressure / MPa")
plot_pressure!(ax, tab, name = tab_name)
plot_pressure!(ax, tab2, name = tab2_name, use_scatter = true)

axislegend(position = :rc)

ax = Axis(fig[2, 2], title = "CO₂ in sealing/bound", xlabel = "Time / year", ylabel = "Mass / kg")
plot_other!(ax, tab, name = tab_name)
plot_other!(ax, tab2, name = tab2_name, use_scatter = true)

axislegend(position = :lt)

# display(GLMakie.Screen(), fig)
fig

