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

mutable struct BasisPolicy
    fixed::Dict{Int,Float64}
    varying::Dict{Int,Vector{Float64}}
    constant::Vector{Float64}
end

bp = BasisPolicy(Dict(),Dict(),[])

OPFRecourse.get_opf_solution(m, ωsamples[:,1])
cbasis, rbasis = MathProgBase.getbasis(m.model.internalModel)

basic_indices = Int[]
for i in 1:ref.ngen
    if cbasis[i] == :Basic
        push!(basic_indices, i)
    elseif cbasis[i] == :NonbasicAtLower
        # @assert(m.model.colNames[i] == "p[$i]", m.model.colNames[i])
        bp.fixed[i] = m.model.colLower[i]
    elseif cbasis[i] == :NonbasicAtUpper
        # @assert(m.model.colNames[i] == "p[$i]", m.model.colNames[i])
        bp.fixed[i] = m.model.colUpper[i]
    else
        error("Unrecognised basis status: $(cbasis[i]) at index $i")
    end
end
numbasic = sum(rbasis .!== :Basic)
basiscol = Dict(zip(basic_indices,1:numbasic))
@assert length(basic_indices) == numbasic
@assert issorted(basic_indices)
@assert length(m.model.linconstr) == 2*ref.nline + 1
basis = zeros(Float64, numbasic, numbasic)
ωmatrix = zeros(Float64, numbasic, ref.nuncertain)
c = 0
# basis[end,:] = ones(numbasic)
for i in 1:(2*ref.nline + 1)
    terms = m.model.linconstr[i].terms
    if rbasis[i] == :Basic
        # ignore
    else
        c += 1
        rhs = if rbasis[i] == :NonbasicAtLower
            @assert m.model.internalModel.lb[i] == m.model.linconstr[i].lb
            m.model.internalModel.lb[i]
        elseif rbasis[i] == :NonbasicAtUpper
            @assert m.model.internalModel.ub[i] == m.model.linconstr[i].ub
            m.model.internalModel.ub[i]
        else
            error("Unrecognised basis status: $(rbasis[i]) at constraint $i")
        end
        rhs -= terms.constant
        for (v,coeff) in zip(terms.vars, terms.coeffs)
            if startswith(v.m.colNames[v.col], "p")
                @assert v.m.colNames[v.col] == "p[$(v.col)]" "$i"
                if cbasis[v.col] == :Basic
                    basis[c,basiscol[v.col]] += coeff
                else
                    rhs -= coeff*bp.fixed[v.col]
                end
            else
                @assert v.m.colNames[v.col] == "ω[$(v.col-ref.ngen)]" "$i"
                @assert cbasis[v.col] !== :Basic
                ωmatrix[c,v.col-ref.ngen] -= coeff
            end
        end
        push!(bp.constant, rhs)
    end
end
@assert c == numbasic
@assert length(bp.constant) == length(bp.fixed) + length(bp.varying) == numbasic
sample_ω*(inv(basis)*ωmatrix)' .+ (inv(basis)*bp.constant)' # predictions based on optimal basis
sample_p # optimal solutions (for comparison)