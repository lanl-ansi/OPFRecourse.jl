using OPFRecourse, Gurobi, Distributions, Base.Test

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
    @test isapprox(-JuMP.getvalue(fullccopf.α[1:2,:]), br.linearterms, atol=1e-6)
end

@testset "Higher Volatility (With 2 Optimal Basis)" begin
    @time ref = OPFRecourse.NetworkReference(data_file, bus_prob = 0.6, line_prob = 0.6, σscaling = 0.5);
    
    @time fullccopf = OPFRecourse.FullChanceConstrainedOPF(ref, Gurobi.GurobiSolver(OutputFlag=0));
    @time JuMP.solve(fullccopf.model, method=:Reformulate)

    @time m = OPFRecourse.SingleScenarioOPF(ref, Gurobi.GurobiSolver(OutputFlag=0));
    @time JuMP.solve(m.model)
    
    srand(1234)
    @time scenarios = OPFRecourse.OPFScenarios(ref, m, nsamples = 1000);
    @test unique(scenarios.whichbasis,1) == [2 1; 1 2]
    @test length(scenarios.cbases) == 2
    @test length(scenarios.rbases) == 2

    br1 = OPFRecourse.BasisRecourse(ref, m, scenarios.cbases[2], scenarios.rbases[1]);
    br2 = OPFRecourse.BasisRecourse(ref, m, scenarios.cbases[1], scenarios.rbases[2]);

    @test OPFRecourse.ngenerationviolations(ref, OPFRecourse.get_opf_solution(br1, scenarios.scenarios[1,:])) == 0
    @test OPFRecourse.ngenerationviolations(ref, OPFRecourse.get_opf_solution(br2, scenarios.scenarios[1,:])) == 0
    @test OPFRecourse.ntransmissionviolations(ref, OPFRecourse.get_opf_solution(br1, scenarios.scenarios[1,:]), scenarios.scenarios[1,:]) == 0
    @test OPFRecourse.ntransmissionviolations(ref, OPFRecourse.get_opf_solution(br2, scenarios.scenarios[1,:]), scenarios.scenarios[1,:]) == 1

    ensemble = OPFRecourse.EnsembleRecourse(ref, fullccopf, [br1, br2]);
    for i in 1:size(scenarios.scenarios,1)
        @test isapprox(
            OPFRecourse.reproject(
                ref,
                OPFRecourse.get_opf_solution(ensemble, scenarios.scenarios[i,:]),
            ),
            scenarios.solutions[i,:],
            atol = 1e-4
        )
    end
end

@testset "Construction of NetworkReference object" begin
    for f in [
        "pglib_opf_case3_lmbd.m",
        "pglib_opf_case5_pjm.m",
        "pglib_opf_case14_ieee.m",
        "pglib_opf_case24_ieee_rts.m",
        "pglib_opf_case30_as.m",
        "pglib_opf_case30_fsr.m",
        "pglib_opf_case30_ieee.m",
        "pglib_opf_case39_epri.m",
        "pglib_opf_case57_ieee.m",
        "pglib_opf_case73_ieee_rts.m",
        "pglib_opf_case89_pegase.m",
        "pglib_opf_case118_ieee.m",
        "pglib_opf_case162_ieee_dtc.m",
        "pglib_opf_case200_pserc.m",
        "pglib_opf_case240_pserc.m",
        "pglib_opf_case300_ieee.m",
        # "pglib_opf_case1354_pegase.m",
        # "pglib_opf_case1888_rte.m",
        # "pglib_opf_case1951_rte.m",
        # "pglib_opf_case2383wp_k.m",
        # "pglib_opf_case2736sp_k.m",
        # "pglib_opf_case2737sop_k.m",
        # "pglib_opf_case2746wop_k.m",
        # "pglib_opf_case2746wp_k.m",
        # "pglib_opf_case2848_rte.m",
        # "pglib_opf_case2868_rte.m",
        # "pglib_opf_case2869_pegase.m",
        # "pglib_opf_case3012wp_k.m",
        # "pglib_opf_case3120sp_k.m"
    ]
        @testset "$f" begin
            @time ref = OPFRecourse.NetworkReference(string(Pkg.dir(),"/OPFRecourse/test/data/pglib-opf/", f))
            @time m = OPFRecourse.SingleScenarioOPF(ref, Gurobi.GurobiSolver(OutputFlag=0));
            @time JuMP.solve(m.model)
        end
    end

    for f in [
        "pglib_opf_case3_lmbd.m",
        "pglib_opf_case5_pjm.m",
        "pglib_opf_case14_ieee.m",
        "pglib_opf_case24_ieee_rts.m",
        "pglib_opf_case30_as.m",
        "pglib_opf_case30_fsr.m",
        "pglib_opf_case30_ieee.m",
        "pglib_opf_case39_epri.m",
        "pglib_opf_case57_ieee.m",
        "pglib_opf_case73_ieee_rts.m",
        "pglib_opf_case89_pegase.m"
    ]
        @testset "$f" begin
            @time ref = OPFRecourse.NetworkReference(string(Pkg.dir(),"/OPFRecourse/test/data/pglib-opf/", f))
            @time ccopf = OPFRecourse.ChanceConstrainedOPF(ref, Gurobi.GurobiSolver(OutputFlag=0, TimeLimit=30));
            @time JuMP.solve(ccopf.model, method=:Cuts, silent=true)
        end
    end
end
