using Jutul, JutulDarcy, HYPRE, CSP11
include("get_RSC_basenames.jl")
allcases = get_RSC_basenames(grids=[:C, :CC], resolutions=["10k"])
allcases = [allcases]

# kgrad_to_run = [:tpfa, :avgmpfa]
kgrad_to_run = [:tpfa]

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
            pth = jutul_output_path(name, subfolder = "csp11_rsc")
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
                output_states = true,
                output_reports = true,
                info_level = 1
            );
        end
    end
end
