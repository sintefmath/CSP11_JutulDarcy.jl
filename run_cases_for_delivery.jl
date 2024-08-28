using Jutul, JutulDarcy, HYPRE, CSP11
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
        );
        pth = jutul_output_path(name, subfolder = "csp11_delivery")
        hook, csv_pth, f = CSP11.get_reporting_hook(pth, domain);
        simulate_reservoir(case,
            output_path = pth,
            max_timestep = 1*si_unit(:year),
            report_level = 1,
            post_ministep_hook = hook,
            relaxation = SimpleRelaxation(),
            max_nonlinear_iterations = 20,
            max_timestep_cuts = 100,
            info_level = 1
        );
    end
end
