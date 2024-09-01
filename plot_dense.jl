using CSV, DataFrames, GLMakie


paths = String[]
pnames = String[]

function addcase!(pth, name = "Case $(length(paths)+1)")
    push!(paths, pth)
    push!(pnames, name)
end

specase = :b
addcase!(joinpath(@__DIR__, "data/compare/b/spe11b_spatial_map_1000y.csv"))
addcase!(joinpath(@__DIR__, "data/compare/generated/spe11b_spatial_map_1000y.csv"))

# specase = :c
# pth = "/media/moyner/Data/jutul_output/csp11_delivery/spe11c_c_10x10x10_thermal_cv_tpfa/dense/spe11b_spatial_map_1000y.csv"


# specase = :c
# pth = joinpath(@__DIR__, "data/compare/generated/spe11c_spatial_map_1000y.csv")


# specase = :c
# pth = joinpath(@__DIR__, "data/compare/c/spe11c_spatial_map_1000y.csv")

# paths = [pth]

if specase == :b
    dims = [840, 120]
else
    @assert specase == :c
    dims = [168, 100, 120]
end

function get_data(pth)
    df = CSV.read(pth, DataFrame)
    ##
    f = " total mass CO2 [kg]"
    f = " total mass CO2 [kg]"
    f = " phase mass density water [kg/m3]"
    # f = "# x [m]"
    f = " gas saturation [-]"
    # f = 
    val = df[:, f]
    val = reshape(val, dims...)
    @info "Sum" sum(val)
    return val
end

data = map(get_data, paths);
##
ndata = length(data)
fig = Figure(size = (1200, 1000))

begin
    cmin = Inf
    cmax = -Inf
    for (i, val) in enumerate(data)
        if specase == :b
            f! = heatmap!
        else
            f! = volume!
        end
        cmin = min(minimum(val), cmin)
        cmax = max(maximum(val), cmax)
        ax = Axis(fig[1, i], title = pnames[i])
        f!(ax, val)
    end
end
Colorbar(fig[1, ndata+1], colorrange = (cmin, cmax))
fig
