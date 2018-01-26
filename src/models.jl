mutable struct SingleScenarioOPF
    model::JuMP.Model
    p::Vector{JuMP.Variable}
    ω::Vector{JuMP.Variable}
end

function SingleScenarioOPF(
        ref::NetworkReference,
        solver::MathProgBase.AbstractMathProgSolver
    )
    model = JuMP.Model(solver=solver)
    JuMP.@variable(model, ref.gen[i].pmin <= p[i in 1:ref.ngen] <= ref.gen[i].pmax, start=ref.gen[i].pstart)
    JuMP.@variable(model,                    ω[i in 1:ref.nbus])
    JuMP.@expression(model, busvalue[i in 1:ref.nbus],
        sum(p[g] for g in ref.bus[i].gens) + ω[i] - ref.bus[i].pd - ref.bus[i].gs 
    )
    lineflow(l) = ref.line[l].β*(
        θ(ref,busvalue,ref.line[l].frombus) - θ(ref,busvalue,ref.line[l].tobus)
    )
    JuMP.@constraints model begin
        [l in 1:ref.nline], lineflow(l) <= ref.line[l].rate
        [l in 1:ref.nline], lineflow(l) >= -ref.line[l].rate
        0 == sum(sum(p[g] for g in ref.bus[i].gens) + ω[i] - ref.bus[i].pd - ref.bus[i].gs
                 for i in 1:ref.nbus)
    end
    JuMP.@objective(model, Min, cost(ref, p))
    SingleScenarioOPF(model,p,ω)
end

SingleScenarioOPF(filename::String; kwargs...) =
    SingleScenarioOPF(PM.build_ref(PM.parse_file(filename)); kwargs...)

SingleScenarioOPF(ref::Dict{Symbol,Any}; kwargs...) =
    SingleScenarioOPF(NetworkReference(ref); kwargs...)

function get_opf_solution(opf::SingleScenarioOPF, ω)
    for i in eachindex(ω); JuMP.fix(opf.ω[i], ω[i]) end
    @assert JuMP.solve(opf.model) == :Optimal
    JuMP.getvalue(opf.p)
end
