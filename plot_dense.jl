using CSV, DataFrames, GLMakie


paths = String[]
pnames = String[]

function addcase!(pth, name = "Case $(length(paths)+1)")
    push!(paths, pth)
    push!(pnames, name)
end

if false
    specase = :b
    # addcase!(joinpath(@__DIR__, "data/compare/b/spe11b_spatial_map_1000y.csv"))
    # addcase!(joinpath(@__DIR__, "data/compare/generated/spe11b_spatial_map_1000y.csv"))

    addcase!(
        raw"D:\jobb\spe11_delivery\aug31\delivery_sintef\delivery\spe11b\result1\spe11b_spatial_map_1000y.csv",
        "TPFA 819x117"
    )

    addcase!(
        raw"D:\jobb\spe11_delivery\aug31\delivery_sintef\delivery\spe11b\result2\spe11b_spatial_map_1000y.csv",
        "MPFA 819x117"
    )

    
    addcase!(
        raw"D:\jobb\spe11_delivery\aug31\delivery_sintef\delivery\spe11b\result3\spe11b_spatial_map_1000y.csv",
        "TPFA 400x60"
    )

    addcase!(
        raw"D:\jobb\spe11_delivery\aug31\delivery_sintef\delivery\spe11b\result4\spe11b_spatial_map_1000y.csv",
        "MPFA 400x60"
    )


else
    specase = :c
    # addcase!(joinpath(@__DIR__, "data/compare/generated/spe11c_spatial_map_1000y.csv"))

    addcase!(
        raw"D:\jobb\spe11_delivery\aug31\delivery_sintef\delivery\spe11c\result1\spe11c_spatial_map_1000y.csv",
        "TPFA 170x100x50"
    )

    addcase!(
        raw"D:\jobb\spe11_delivery\aug31\delivery_sintef\delivery\spe11c\result2\spe11c_spatial_map_1000y.csv",
        "MPFA 170x100x50"
    )

    
    addcase!(
        raw"D:\jobb\spe11_delivery\aug31\delivery_sintef\delivery\spe11c\result3\spe11c_spatial_map_1000y.csv",
        "TPFA 170x100100"
    )

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

function get_data(en)
    i, pth = en
    df = CSV.read(pth, DataFrame)
    ##
    f = " total mass CO2 [kg]"
    f = " total mass CO2 [kg]"
    # f = " phase mass density water [kg/m3]"
    # f = "# x [m]"
    # f = " # y [m]"
    # if i < 4
    #     f = " y [m]"
    # else
    #     f = " z [m]"
    # end
    # f = " gas saturation [-]"
    # f = 
    x = df[:, "# x [m]"]
    y = df[:, " y [m]"]
    z  = df[:, " z [m]"]

    Lx = maximum(x)
    Ly = maximum(y)
    Lz = maximum(z)
    
    val = df[:, f]
    function sortfunction(i)
        xi = x[i]
        yi = y[i]
        zi = z[i]

        return zi*(Lx*Ly) + yi*Lx + xi
    end
    ix = sort(eachindex(val), by = sortfunction)

    val = reshape(val[ix], dims...)
    @info "Sum" sum(val)
    return val
end

data = map(get_data, enumerate(paths));
##
ndata = length(data)
fig = Figure(size = (1200, 1000))

begin
    cmin = Inf
    cmax = -Inf
    for (i, val) in enumerate(data)
        if specase == :b
            f! = heatmap!
            A = Axis
        else
            f! = volume!
            A = Axis3
        end
        cmin = min(minimum(val), cmin)
        cmax = max(maximum(val), cmax)
        ax = A(fig[1, i], title = pnames[i])
        f!(ax, val)
    end
end
Colorbar(fig[1, ndata+1], colorrange = (cmin, cmax))
fig
