function get_RSC_basenames(;grids = [:C, :HC, :CC, :PEBI, :QT, :T],
    resolutions=["10k"], specase = :b)
    if specase == :b
        # Dictionary mapping resolution symbols to their specs [cartesian, QT_ref, T_ref]
        resolutionDict = Dict(
            "10k"  => ["140x75", "1_2", "1_3"],
            "50k"  => ["500x100", "0_38", "0_54"],
            "100k" => ["840x120", "0_25", "0_37"],
            "200k" => ["1180x170", "0_17", "0_26"],
            "500k" => ["1870x270", "0_103", "0_160"],
            "1M"   => ["2640x380", "0_071", "0_112"]
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
    for specs in resolutionSpecs
        res = specs[1]
        refQT = specs[2]
        refT = specs[3]
        for grid in cartGrids
            push!(allcases, grid*"_"*res)
        end
        for grid in QTGrids
            push!(allcases, grid*refQT)
        end
        for grid in TGrids
            push!(allcases, grid*refT)
        end
    end
    return allcases
end
