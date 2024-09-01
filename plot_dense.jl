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

    addcase!(joinpath(@__DIR__, "data/compare/b/spe11b_spatial_map_1000y.csv"), "OPM")

else
    specase = :c
    # addcase!(joinpath(@__DIR__, "data/compare/generated/spe11c_spatial_map_1000y.csv"))

    addcase!(raw"D:/jobb/jutul_output/csp11_delivery_sep1_test/spe11c_nudge_c_10x10x10_thermal_cv_tpfa/dense/spe11c_spatial_map_1000y.csv")
    if false

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
    end
    addcase!(joinpath(@__DIR__, "data/compare/c/spe11c_spatial_map_1000y.csv"), "OPM")
end
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

function get_corrected_coord_case_c()
    dx = 50.0
    dy = 50.0
    dz = 10.0

    ni = 168
    nj = 120
    nk = 100

    x = Float64[]
    y = Float64[]
    z = Float64[]
    # Original: k, i, j
    for k in 1:nk
        for j in 1:nj
            for i in 1:ni
                xi = (i-0.5)*dx
                yi = (k-0.5)*dy
                zi = (j-0.5)*dz

                push!(x, xi)
                push!(y, yi)
                push!(z, zi)
                # push!(z, yi)
                # push!(y, zi)

            end
        end
    end

    return (x, y, z)
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

    if i < 2
        x = df[:, "# x [m]"]
        y = df[:, " y [m]"]
        z  = df[:, " z [m]"]
    
        x, y, z = get_corrected_coord_case_c()
    else

        x = df[:, "# x [m]"]
        y = df[:, " y [m]"]
        z  = df[:, " z [m]"]
    
        # z = reverse(z)
        # y, z = z, y
        # , z = z, x
    end
    Lx = maximum(x)
    Ly = maximum(y)
    Lz = maximum(z)
    @info "??" Lx Ly Lz
    
    val = df[:, f]
    function sortfunction(i)
        xi = x[i]
        yi = y[i]
        zi = z[i]
        # return i
        return zi*(Lx*Ly) + yi*Lx + xi
    end
    ix = sort(eachindex(val), by = sortfunction)
    # ix = eachindex(val)
    # if i == 1
    #     val = reshape(val[ix], dims[1], dims[3], dims[2])
    # else
        val = reshape(val[ix], dims...)
    # end
    @info "Sum" sum(val)
    return val
end

data = map(get_data, enumerate(paths));

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
