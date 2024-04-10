using CSV, DataFrames
using Jutul

name1 = "horizon-cut_100x50_thermal_cv"
name2 = "horizon-cut_100x50_thermal_cv"

output_folder = jutul_output_path("spe11b_$name1", subfolder = "csp11")
output_folder2 = jutul_output_path("spe11b_$name2", subfolder = "csp11")

file_name = "spe11b_time_series.csv"

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
tab_name = "Case1"

tab2 = read_file(joinpath(output_folder2, file_name))
tab2_name = "Case2"
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
    else
        f! = lines!
    end
    t = tab.t
    X = tab[reg]
    f!(ax, t, X.mob, label = genlabel("mobile"); SCATTER_ARG...)
    f!(ax, t, X.imm, label = genlabel("immobile"); SCATTER_ARG...)
    f!(ax, t, X.diss, label = genlabel("dissolved"); SCATTER_ARG...)
    f!(ax, t, X.seal, label = genlabel("in sealing"); SCATTER_ARG...)
    f!(ax, t, X.mob + X.imm + X.diss, label = genlabel("total"); SCATTER_ARG...)
end

function plot_pressure!(ax, tab; use_scatter = false, name = "", kwarg...)
    if isnothing(tab)
        return tab
    end
    if use_scatter
        f! = scatter!
    else
        f! = lines!
    end
    function genlabel(x)
        if name == ""
            return x
        else
            return "$x ($name)"
        end
    end
    t = tab.t
    f!(t, tab.p1/1e6, label = genlabel("P1"); SCATTER_ARG...)
    f!(t, tab.p2/1e6, label = genlabel("P2"); SCATTER_ARG...)
end

function plot_other!(ax, tab; use_scatter = false, name = "", kwarg...)
    if isnothing(tab)
        return tab
    end
    if use_scatter
        f! = scatter!
    else
        f! = lines!
    end
    function genlabel(x)
        if name == ""
            return x
        else
            return "$x ($name)"
        end
    end
    t = tab.t
    f!(t, tab.sealTot, label = genlabel("seal"); SCATTER_ARG...)
    f!(t, tab.boundTot, label = genlabel("bound"); SCATTER_ARG...)
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

