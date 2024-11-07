using Jutul, JutulDarcy, HYPRE, CSP11

cases_b = [
    "b_400x60",
    "b_819x117",
    "b_2640x380"
]

cases_c = [
    "c_10x10x10",
    "c_85x50x40",
    "c_170x100x50"
]

allcases = [cases_b, cases_c]
allcases = [["b_400x60"]]

kgrad_to_run = [:tpfa, :avgmpfa]
for cases_to_run in allcases
    for kg in kgrad_to_run
        for basename in cases_to_run
            specase = Symbol(basename[1])
            case, name = setup_spe11_case_from_mrst_grid(basename,
                case = specase,
                thermal = true,
                use_reporting_steps = true,
                kgrad = kg
            );
            domain = reservoir_domain(case.model)
            nc = number_of_cells(domain)
            println("Running case $name with $nc cells")
            pth = jutul_output_path(name, subfolder = "csp11_delivery_aug30")
            hook, csv_pth, f = CSP11.get_reporting_hook(pth, domain, specase = specase)
            simulate_reservoir(case,
                output_path = pth,
                max_timestep = 1*si_unit(:year),
                report_level = 0,
                tol_cnv = 0.01,
                post_ministep_hook = hook,
                relaxation = SimpleRelaxation(),
                max_nonlinear_iterations = 20,
                max_timestep_cuts = 100,
                output_states = false,
                output_reports = false,
                info_level = 1
            );
        end
    end
end
