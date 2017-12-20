using OPFRecourse, Gurobi, Base.Test

data_file = string(Pkg.dir(),"/OPFRecourse/test/data/nesta_case30_ieee.m")

@testset "Lower Volatility (Only 1 Optimal Basis)" begin
    @time ref = OPFRecourse.NetworkReference(data_file, bus_prob = 0.95, line_prob = 0.95, σscaling = 0.05);
    
    @time ccopf = OPFRecourse.ChanceConstrainedOPF(ref, Gurobi.GurobiSolver(OutputFlag=0));
    @time JuMP.solve(ccopf.model, method=:Reformulate)

    @time fullccopf = OPFRecourse.FullChanceConstrainedOPF(ref, Gurobi.GurobiSolver(OutputFlag=0));
    @time JuMP.solve(fullccopf.model, method=:Reformulate)

    @time m = OPFRecourse.SingleScenarioOPF(ref, Gurobi.GurobiSolver(OutputFlag=0));
    @time JuMP.solve(m.model)

    @time scenarios = OPFRecourse.OPFScenarios(ref, m, nsamples = 1000);
    @test length(OPFRecourse.get_opf_solution(ccopf, scenarios.scenarios[1,:])) == ref.ngen
    @test length(OPFRecourse.get_opf_solution(fullccopf, scenarios.scenarios[1,:])) == ref.ngen
    @test length(OPFRecourse.get_opf_solution(m, scenarios.scenarios[1,:])) == ref.ngen

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
end

@testset "Higher Volatility (With 2 Optimal Basis)" begin
    @time ref = OPFRecourse.NetworkReference(data_file, bus_prob = 0.6, line_prob = 0.6, σscaling = 0.5);
    
    @time fullccopf = OPFRecourse.FullChanceConstrainedOPF(ref, Gurobi.GurobiSolver(OutputFlag=0));
    @time JuMP.solve(fullccopf.model, method=:Reformulate)

    @time m = OPFRecourse.SingleScenarioOPF(ref, Gurobi.GurobiSolver(OutputFlag=0));
    @time JuMP.solve(m.model)
    
    srand(1234)
    @time scenarios = OPFRecourse.OPFScenarios(ref, m, nsamples = 1000);
    @test length(scenarios.cbases) == 2
    @test length(scenarios.rbases) == 2

    uniquebases = unique(scenarios.whichbasis,1)

    basisrecourses = [
        OPFRecourse.BasisRecourse(ref, m,
            scenarios.cbases[uniquebases[i,1]],
            scenarios.rbases[uniquebases[i,2]]
        )
        for i in 1:size(uniquebases,1)
    ]

    ensemble = OPFRecourse.EnsembleRecourse(ref, fullccopf, basisrecourses)
    for i in 1:size(scenarios.scenarios,1)
        @test isapprox(
            OPFRecourse.get_opf_solution(ensemble, scenarios.scenarios[i,:]),
            scenarios.solutions[i,:],
            atol = 1e-4
        )
    end
end
