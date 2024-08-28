using CSV, DataFrames, GLMakie



specase = :b
pth = joinpath(@__DIR__, "data/compare/b/spe11b_spatial_map_1000y.csv")

# specase = :b
# pth = joinpath(@__DIR__, "data/compare/generated/spe11b_spatial_map_1000y.csv")


# specase = :c
# pth = joinpath(@__DIR__, "data/compare/c/spe11c_spatial_map_1000y.csv")


if specase == :b
    dims = [840, 120]
else
    @assert specase == :c
    dims = [168, 100, 120]
end
df = CSV.read(pth, DataFrame)
##
f = " total mass CO2 [kg]"
# f = 
val = df[:, f]
val = reshape(val, dims...)
@info "Sum" sum(val)

if specase == :b
    fig, ax, plt = heatmap(val)
else
    fig, ax, plt = volume(val)
end
Colorbar(fig[1, 2], plt)
fig