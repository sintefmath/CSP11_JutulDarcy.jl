using Jutul, JutulDarcy, HYPRE, CSP11, GLMakie

# basename = "cPEBI_2640x380"
# basename = "cPEBI_819x117"
basename = "horizon-cut_100x50"; specase = :b
basename = "struct20x20x20"; specase = :c
# basename = "struct50x50x50"; specase = :c

# basename = "struct100x100x100"; specase = :c

case, name = setup_spe11_case_from_mrst_grid(basename, case = specase, thermal = true);
pth = jutul_output_path(name, subfolder = "csp11")
domain = reservoir_domain(case.model)

hook, csv_pth, f = CSP11.get_reporting_hook(pth, domain);
##
sim, config = setup_reservoir_simulator(case,
    output_path = pth,
    max_timestep = 1*si_unit(:year),
    report_level = 1,
    tol_cnv = 0.001,
    post_ministep_hook = hook,
    relaxation = SimpleRelaxation(),
    max_nonlinear_iterations = 20,
    max_timestep_cuts = 100,
);
##
restart = false #'true' will continue sim where it left off
res = simulate_reservoir(case,
    simulator = sim,
    config = config,
    restart = restart,
);
ws, states = res;
##write to mrst format
mrst_output_path = pth*"_mrst"
JutulDarcy.write_reservoir_simulator_output_to_mrst(case.model[:Reservoir], states, res.result.reports, case.forces, mrst_output_path)

## Plotting
model = case.model
using GLMakie
plot_facies = false
fig = plot_reservoir(model, states, title = name, transparency = false, alpha = 0.5 + 0.5*!plot_facies)
if plot_facies
    ax = fig.current_axis[]
    plot_cell_data!(ax, domain.representation, Float64.(domain[:satnum]), transparency = false, alpha = 0.5, colormap = :gray1)
end
fig
