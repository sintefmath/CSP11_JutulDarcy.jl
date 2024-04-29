module CSP11
    using Jutul, JutulDarcy, StaticArrays, MultiComponentFlash
    using Statistics, MAT

    export csp_phase_property_table
    export read_solubility_table
    export read_component_table
    export setup_state0_csp11
    export setup_spe11_case_from_mrst_grid
    export reservoir_domain_csp11
    export reservoir_domain_and_wells_csp11
    export setup_reservoir_model_csp11

    include("reporting.jl")
    include("variables.jl")
    # include("kvalues.jl")
    # include("props.jl")
    # include("reading.jl")
    include("setup.jl")

end # module CSP11
