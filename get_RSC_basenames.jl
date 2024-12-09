function get_RSC_basenames(;grids = [:C, :HC, :CC, :PEBI, :QT, :T],
    resolutions=["10k"], specase = :b)
    if specase == :b
        # Dictionary mapping resolution symbols to their specs [cartesian, QT_ref, T_ref]
        resolutionDict = Dict(
            "10k"  => ["140x75", "1_2", "1_25"],
            "50k"  => ["500x100", "0_54", "0_52"],
            "100k" => ["840x120", "0_38", "0_36"],
            "200k" => ["1180x170", "0_28", "0_25"],
            "500k" => ["1870x270", "0_175", "0_16"],
            "1M"   => ["2640x380", "0_125", "0_112"]
        )
    else
        error("specase $specase not yet supported")
    end
    
    # Get specs for requested resolutions
    resolutionSpecs = [resolutionDict[res] for res in resolutions]
    
    cartGrids = []
    QTGrids = []
    TGrids = []
    
    # Sort grids into appropriate categories
    for grid in grids
        if grid in [:C, :HC, :CC, :PEBI]
            push!(cartGrids, "b_" * string(grid))
        elseif grid == :QT
            push!(QTGrids, "b_QT")
        elseif grid == :T
            push!(TGrids, "b_T")
        end
    end

    allcases = []
    gridcases = []  # For MATLAB-style grid names
    
    for specs in resolutionSpecs
        res = specs[1]
        refQT = specs[2]
        refT = specs[3]
        
        # Parse nx and nz from resolution string (e.g., "140x75")
        nx, nz = parse.(Int, split(res, "x"))
        
        # Parse QT and T values from their reference strings
        qt = parse(Float64, replace(refQT, "_" => "."))
        t = parse(Float64, replace(refT, "_" => "."))
        
        for grid in cartGrids
            push!(allcases, grid*"_"*res)
            
            # Add corresponding MATLAB grid name based on grid type
            if occursin("b_C_", grid)
                push!(gridcases, "struct$(nx)x$(nz)")
            elseif occursin("b_HC_", grid)
                push!(gridcases, "horz_ndg_cut_PG_$(nx)x$(nz)")
            elseif occursin("b_CC_", grid)
                push!(gridcases, "cart_ndg_cut_PG_$(nx)x$(nz)")
            elseif occursin("b_PEBI_", grid)
                push!(gridcases, "cPEBI_$(nx)x$(nz)")
            end
        end
        
        for grid in QTGrids
            push!(allcases, grid*refQT)
            push!(gridcases, "gq_pb$qt")
        end
        
        for grid in TGrids
            push!(allcases, grid*refT)
            push!(gridcases, "5tetRef$t")
        end
    end
    
    return allcases, gridcases
end