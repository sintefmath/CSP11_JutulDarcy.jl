using Jutul, JutulDarcy, CSP11

specase = :b
dims = (400, 60)
# specase = :a; dims = (280, 120)
# specase = :c; dims = (85, 50, 40)

case, name = setup_spe11_case(dims; case = specase)
pth = jutul_output_path(name; subfolder = "csp11")
hook, = CSP11.get_reporting_hook(pth, reservoir_domain(case.model))
res = simulate_reservoir(case;
    output_path = pth,
    max_timestep = specase == :a ? 10*si_unit(:minute) : CSP11.spe11_year,
    post_ministep_hook = hook,
    tol_cnv = 0.01,
    relaxation = SimpleRelaxation(),
    max_nonlinear_iterations = 20,
    max_timestep_cuts = 100,
    info_level = 3
)
