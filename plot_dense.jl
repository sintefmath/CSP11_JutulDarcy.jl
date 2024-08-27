using CSV, DataFrames, GLMakie

specase = :b
pth = joinpath(@__DIR__, "data/compare/b/spe11b_spatial_map_1000y.csv")

specase = :c
pth = joinpath(@__DIR__, "data/compare/c/spe11c_spatial_map_1000y.csv")


if specase == :b
    dims = [840, 120]
else
    @assert specase == :c
    dims = [168, 100, 120]
end
df = CSV.read(pth, DataFrame)
##
rhow = df[:, " phase mass density water [kg/m3]"]
rhow = reshape(rhow, dims...)

if specase == :b
    heatmap(rhow)
else
    volume(rhow)
end
