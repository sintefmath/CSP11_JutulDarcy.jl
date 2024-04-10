using Jutul, JutulDarcy, HYPRE, CSP11, GLMakie

basename = "cPEBI_2640x380"
basename = "cPEBI_819x117"
basename = "horizon-cut_100x50"

case, name = setup_spe11_caseb_from_mrst_grid(basename, thermal = true);
pth = jutul_output_path(name, subfolder = "csp11")
domain = reservoir_domain(case.model)

hook, csv_pth, f = CSP11.get_reporting_hook(pth, domain)

res = simulate_reservoir(case,
    output_path = pth,
    max_timestep = 1*si_unit(:year),
    report_level = 0,
    tol_cnv = 0.001,
    post_ministep_hook = hook,
    relaxation = SimpleRelaxation(),
    max_nonlinear_iterations = 20,
    max_timestep_cuts = 100,
    info_level = 1
)
ws, states = res;

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
