module CSP11
    using Jutul, JutulDarcy, StaticArrays, MultiComponentFlash
    using Statistics, MAT

    export csp_phase_property_table
    export read_solubility_table
    export read_component_table
    export setup_state0_csp11
    export setup_spe11_case
    export setup_spe11_domain
    export setup_spe11_case_from_mrst_grid
    export spe11_facies, spe11c_elevation_offset
    export reservoir_domain_csp11
    export reservoir_domain_and_wells_csp11
    export setup_reservoir_model_csp11

    include("reporting.jl")
    include("variables.jl")
    include("geometry.jl")
    # include("kvalues.jl")
    # include("props.jl")
    # include("reading.jl")
    include("setup.jl")

end # module CSP11
