using Jutul, JutulDarcy, HYPRE, CSP11
include("get_RSC_basenames.jl")
allcases,  = get_RSC_basenames(grids=[:C, :HC, :CC, :PEBI, :QT, :T], resolutions=["100k"])
# allcases,  = get_RSC_basenames(grids=[:C], resolutions=["10k"])

allcases = [allcases]

restart = true
save_mrst_output = true
# kgrad_to_run = [:tpfa, :avgmpfa]
# kgrad_to_run = [:tpfa, :avgmpfa, :ntpfa]
kgrad_to_run = [:ntpfa]


if save_mrst_output
    extra_outputs = [:RelativePermeabilities,
                     :LiquidMassFractions,
                     :VaporMassFractions,
                     :PhaseMassDensities]
else
    extra_outputs = false
end
for cases_to_run in allcases
    for kg in kgrad_to_run
        for basename in cases_to_run
            specase = Symbol(basename[1])
            case, name = setup_spe11_case_from_mrst_grid(basename,
                case = specase,
                thermal = false,
                use_reporting_steps = true,
                kgrad = kg,
                extra_outputs = extra_outputs
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
                output_states = false,
                output_reports = false,
                info_level = 0,
                restart = restart
            );
            if save_mrst_output
                mrst_output_path = pth*"_mrst"
                res = simulate_reservoir(case,
                    output_path = pth,
                    restart = true
                );
                ws, states = res;
                JutulDarcy.write_reservoir_simulator_output_to_mrst(case.model[:Reservoir], states, res.result.reports, case.forces, mrst_output_path)
            end
        end
    end
end
