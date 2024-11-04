using Jutul, JutulDarcy, HYPRE, CSP11, GLMakie

# basename = "cPEBI_2640x380"
# basename = "cPEBI_819x117"; specase = :b
basename = "horizon-cut_100x50"; specase = :b

basename = "b_400x60"; specase = :b

basename = "c_85x50x40"; specase = :c

# basename = "test_cart_100_50_50"; specase = :c
# basename = "cart_a_30_20"
# basename = "struct20x20x20"; specase = :c

# basename = "c_10x10x10"; specase = :c

# basename = "cart_a_30_20"; specase = :a

# basename = "struct50x50x50"; specase = :c

# basename = "struct100x100x100"; specase = :c

basename = "nudge_c_170x100x100"; specase = :c
kg = :avgmpfa
# kg = :ntpfa
# kg = :tpfa_test
# kg = :tpfa

case, name = setup_spe11_case_from_mrst_grid(basename,
    case = specase,
    thermal = true,
    nstep_initialization = 0,
    use_reporting_steps = false,
    # dz_max = 0.1,
    kgrad = kg
);
pth = jutul_output_path(name, subfolder = "csp11")
domain = reservoir_domain(case.model)

##
hook, csv_pth, f = CSP11.get_reporting_hook(pth, domain);
sim, config = setup_reservoir_simulator(case,
    output_path = pth,
    max_timestep = 1*si_unit(:year),
    report_level = 0,
    tol_cnv = 0.01,
    nonlinear_tolerance = Inf,
    post_ministep_hook = hook,
    relaxation = SimpleRelaxation(),
    max_nonlinear_iterations = 20,
    max_timestep_cuts = 100,
    target_ds = 0.1,
    info_level = 3
);
##
res = simulate_reservoir(case,
    simulator = sim,
    config = config,
);
ws, states = res;
##write to mrst format
# mrst_output_path = pth*"_mrst"
# JutulDarcy.write_reservoir_simulator_output_to_mrst(case.model[:Reservoir], states, res.result.reports, case.forces, mrst_output_path)

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
##
if false
    states_ref = CSP11.map_to_reporting_grid(case, states);
    ##
    if specase == :b
        x = states_ref[end][:sg]
        x = states_ref[end][:co2_mass]
        x = reshape(x, 840, 120)
        heatmap(x)
    end
end
##
pth = joinpath(@__DIR__, "data", "compare", "generated")
mkpath(pth)
CSP11.write_reporting_grid(case, states, pth, specase)