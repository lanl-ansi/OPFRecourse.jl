using OPFRecourse, Gurobi, Distributions, Base.Test

data_file = string(Pkg.dir(),"/AlternatingOPF/test/data/nesta_case30_ieee_prob.m")

@testset "Lower Volatility (Only 1 Optimal Basis)" begin
    @time ref = PowerModels.build_ref(PowerModels.parse_file(data_file));
    for (i,u) in ref[:uncertainty]; u["std"] /= 100.0 end
    ref[:branch][2]["rate_a"] = 0.83 # line tightening
    ref = OPFRecourse.NetworkReference(ref, bus_prob = 0.95, line_prob = 0.95)
    @time ccopf = OPFRecourse.ChanceConstrainedOPF(ref, Gurobi.GurobiSolver(OutputFlag=0));
    @time JuMP.solve(ccopf.model, method=:Reformulate)
    @time fullccopf = OPFRecourse.FullChanceConstrainedOPF(ref, Gurobi.GurobiSolver(OutputFlag=0));
    @time JuMP.solve(fullccopf.model, method=:Reformulate)
    @time m = OPFRecourse.SingleScenarioOPF(ref, Gurobi.GurobiSolver(OutputFlag=0));
    @time JuMP.solve(m.model)

    @time scenarios = OPFRecourse.OPFScenarios(ref, m, nsamples = 1000);
    br = OPFRecourse.BasisRecourse(ref, m, scenarios.cbases[1], scenarios.rbases[1]);

    @test unique(scenarios.whichbasis,1) == [1 1]
    @test length(scenarios.cbases) == 1
    @test length(scenarios.rbases) == 1
    for i in 1:size(scenarios.scenarios,1)
        @test isapprox(
            OPFRecourse.reproject(
                ref,
                OPFRecourse.get_opf_solution(br, scenarios.scenarios[i,:])
            ),
            scenarios.solutions[i,:],
            atol = 1e-6
        )
    end
    @test isapprox(-JuMP.getvalue(fullccopf.α[1:2,:]), br.linearterms, atol=1e-6)
end

@testset "Higher Volatility (With 2 Optimal Basis)" begin
    @time ref = PowerModels.build_ref(PowerModels.parse_file(data_file));
    for (i,u) in ref[:uncertainty]; u["std"] /= 50.0 end
    ref[:branch][2]["rate_a"] = 0.83 # line tightening
    ref = OPFRecourse.NetworkReference(ref, bus_prob = 0.95, line_prob = 0.95)
    @time m = OPFRecourse.SingleScenarioOPF(ref, Gurobi.GurobiSolver(OutputFlag=0));
    @time JuMP.solve(m.model)
    @time scenarios = OPFRecourse.OPFScenarios(ref, m, nsamples = 1000);

    @test unique(scenarios.whichbasis,1) == [1 2; 1 1]
    @test length(scenarios.cbases) == 1
    @test length(scenarios.rbases) == 2
    br1 = OPFRecourse.BasisRecourse(ref, m, scenarios.cbases[1], scenarios.rbases[1]);
    br2 = OPFRecourse.BasisRecourse(ref, m, scenarios.cbases[1], scenarios.rbases[2]);

    @test OPFRecourse.ngenerationviolations(ref, OPFRecourse.get_opf_solution(br1, scenarios.scenarios[1,:])) == 0
    @test OPFRecourse.ngenerationviolations(ref, OPFRecourse.get_opf_solution(br2, scenarios.scenarios[1,:])) == 0
    @test OPFRecourse.ntransmissionviolations(ref, OPFRecourse.get_opf_solution(br1, scenarios.scenarios[1,:]), scenarios.scenarios[1,:]) == 1
    @test OPFRecourse.ntransmissionviolations(ref, OPFRecourse.get_opf_solution(br2, scenarios.scenarios[1,:]), scenarios.scenarios[1,:]) == 0
end
