using Jutul, JutulDarcy, HYPRE, CSP11

basename = "cPEBI_2640x380"; specase = :b
basename = "cPEBI_819x117"; specase = :b
basename = "horizon-cut_100x50"; specase = :b

kg = :tpfa

case, name = setup_spe11_case_from_mrst_grid(basename, case = specase, thermal = true, nstep_initialization = 0, use_reporting_steps = false, kgrad = kg);
pth = jutul_output_path(name, subfolder = "csp11")
domain = reservoir_domain(case.model)

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

  
