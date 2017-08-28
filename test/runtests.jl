using OPFRecourse, Gurobi, Base.Test

data_file = string(Pkg.dir(),"/AlternatingOPF/test/data/nesta_case30_ieee_prob.m")
@time ref = PowerModels.build_ref(PowerModels.parse_file(data_file));
for (i,u) in ref[:uncertainty]; u["std"] /= 100.0 end # scale down volatility
ref[:branch][2]["rate_a"] = 0.83 # line tightening
ref = OPFRecourse.NetworkReference(ref, bus_prob = 0.95, line_prob = 0.95)
@time ccopf = OPFRecourse.ChanceConstrainedOPF(ref, Gurobi.GurobiSolver());
@time JuMP.solve(ccopf.model, method=:Reformulate)

@time m = OPFRecourse.SingleScenarioOPF(ref, Gurobi.GurobiSolver(OutputFlag=0));
@time JuMP.solve(m.model)

ω = Distributions.MvNormal([ref.refω[i]["mean"] for i in 1:ref.nuncertain], Array(ref.sqrtΣ^2))
nsamples = 20_000
ωsamples = rand(ω, nsamples)
status = Array{Symbol}(nsamples);
soln_p = zeros(nsamples, ref.ngen);
JuMP.build(m.model);
nv = Gurobi.num_vars(m.model.internalModel.inner)
nc = Gurobi.num_constrs(m.model.internalModel.inner)
cbases = Dict{NTuple{nv, Symbol},Vector{Int}}()
rbases = Dict{NTuple{nc, Symbol},Vector{Int}}()
bases = Dict{NTuple{nv+nc, Symbol},Vector{Int}}()
noptimal = 0

for i in 1:nsamples
    status[i] = OPFRecourse.get_opf_solution(m, ωsamples[:,i])

    if status[i] == :Optimal
        # if the scenario was feasible, keep track of basis
        soln_p[i,:] = JuMP.getvalue(m.p)
        noptimal += 1
        cbasis, rbasis = MathProgBase.getbasis(m.model.internalModel)
        cbasis, rbasis = tuple(cbasis...), tuple(rbasis...)
        basis = tuple(cbasis..., rbasis...)
        cbases[cbasis] = get(cbases, cbasis, Int[])
        rbases[rbasis] = get(rbases, rbasis, Int[])
        bases[basis] = get(bases, basis, Int[])
        push!(cbases[cbasis], noptimal)
        push!(rbases[rbasis], noptimal)
        push!(bases[basis], noptimal)
    end
end
@assert noptimal == sum(status .== :Optimal)

sample_p = Float32.(soln_p[status .== :Optimal, :])
sample_ω = Float32.(ωsamples[:,status .== :Optimal])'

OPFRecourse.get_opf_solution(m, ωsamples[:,1])
cbasis, rbasis = MathProgBase.getbasis(m.model.internalModel)
br = OPFRecourse.BasisRecourse(ref, m, cbasis, rbasis);
@show OPFRecourse.get_opf_solution(br, ωsamples[:,1])
# versus optimal solution (for comparison)
@show sample_p[1,:]
