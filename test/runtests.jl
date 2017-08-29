using OPFRecourse, Gurobi, Distributions, Base.Test

data_file = string(Pkg.dir(),"/AlternatingOPF/test/data/nesta_case30_ieee_prob.m")
@time ref = PowerModels.build_ref(PowerModels.parse_file(data_file));
for (i,u) in ref[:uncertainty]; u["std"] /= 100.0 end # scale down volatility
ref[:branch][2]["rate_a"] = 0.83 # line tightening
ref = OPFRecourse.NetworkReference(ref, bus_prob = 0.95, line_prob = 0.95)
@time ccopf = OPFRecourse.ChanceConstrainedOPF(ref, Gurobi.GurobiSolver());
@time JuMP.solve(ccopf.model, method=:Reformulate)
@time fullccopf = OPFRecourse.FullChanceConstrainedOPF(ref, Gurobi.GurobiSolver());
@time JuMP.solve(fullccopf.model, method=:Reformulate)
@time m = OPFRecourse.SingleScenarioOPF(ref, Gurobi.GurobiSolver(OutputFlag=0));
@time JuMP.solve(m.model)

scenarios = OPFRecourse.OPFScenarios(ref, m, nsamples = 1000)

@test unique(scenarios.whichbasis,1) == [1 1]
@test length(scenarios.cbases) == 1
@test length(scenarios.rbases) == 1
br = OPFRecourse.BasisRecourse(ref, m, scenarios.cbases[1], scenarios.rbases[1]);
for i in 1:size(scenarios.scenarios,1)
    @test isapprox(
        OPFRecourse.get_opf_solution(br, scenarios.scenarios[i,:]),
        scenarios.solutions[i,:]
    )
end
@test isapprox(-JuMP.getvalue(fullccopf.α[1:2,:]), br.linearterms, atol=1e-6)
