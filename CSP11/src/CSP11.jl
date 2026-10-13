module CSP11
    using Jutul, JutulDarcy, StaticArrays

    export setup_state0_csp11
    export setup_spe11_case
    export setup_spe11_domain
    export spe11_facies, spe11c_elevation_offset
    export reservoir_domain_csp11
    export setup_reservoir_model_csp11

    include("specification.jl")
    include("geometry.jl")
    include("setup.jl")
    include("reporting.jl")

end # module CSP11
