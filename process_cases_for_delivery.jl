cases_to_run = [
    "c_10x10x10"
]
kgrad_to_run = [:tpfa, :avgmpfa]
for kg in kgrad_to_run
    for basename in cases_to_run
        specase = Symbol(basename[1])
        case, name = setup_spe11_case_from_mrst_grid(basename,
            case = specase,
            thermal = true,
            use_reporting_steps = true,
            kgrad = kg
        )
        pth = jutul_output_path(name, subfolder = "csp11_delivery")
        dense_pth = joinpath(pth, "dense")
        states, reports = read_results(pth)
        rstates = map(x -> x[:Reservoir], states)

        mkpath(dense_pth)
        CSP11.write_reporting_grid(case, rstates, dense_pth, specase)

        sparse_pth = joinpath(pth, "sparse")
        CSP11.print_performance(sparse_pth, rstates, reports, specase)
    end
end
