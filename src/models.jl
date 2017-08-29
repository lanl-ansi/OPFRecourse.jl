"computes B̃(p + μ + ω - d)"
function busfromRHS(ref, p, omega)
    function busvalue(i)
        result = - ref.bus[i]["pd"] - ref.bus[i]["gs"]
        for g in ref.busgens[i]; result += p[g] end
        for j in 1:ref.nuncertain
            (ref.refω[j]["bus"] == i) && (result += omega[j])
        end
        result 
    end
    [busvalue(i) for i in 1:ref.nbus]
end

"computes B̃(p - α⋅ω + μ + ω - d)"
function busfromRHS(ref, p, α, omega)
    function busvalue(i)
        result = - ref.bus[i]["pd"] - ref.bus[i]["gs"]
        for g in ref.busgens[i]; result += p[g] - α[g]*sum(omega) end
        for j in 1:ref.nuncertain
            (ref.refω[j]["bus"] == i) && (result += omega[j])
        end
        result 
    end
    [busvalue(i) for i in 1:ref.nbus]
end

function fullbusfromRHS(ref, p, α, omega)
    function busvalue(i)
        result = - ref.bus[i]["pd"] - ref.bus[i]["gs"]
        for g in ref.busgens[i]; result += p[g] - sum(α[g,j]*omega[j] for j in 1:ref.nuncertain) end
        for j in 1:ref.nuncertain
            (ref.refω[j]["bus"] == i) && (result += omega[j])
        end
        result 
    end
    [busvalue(i) for i in 1:ref.nbus]
end

"computes B_f*θ where B_f is the matrix of power transfer factors"
Bftheta(ref, θ) =
    [β(ref,l)*(θ[frombus(ref,l)] - θ[tobus(ref,l)]) for l in 1:ref.nline]

lineflow(ref, p, ω)    = Bftheta(ref, ref.π*busfromRHS(ref, p, ω))
lineflow(ref, p, α, ω) = Bftheta(ref, ref.π*busfromRHS(ref, p, α, ω))
fulllineflow(ref, p, α, ω) = Bftheta(ref, ref.π*fullbusfromRHS(ref, p, α, ω))

mutable struct ChanceConstrainedOPF
    model::JuMP.Model
    p
    α
end

function ChanceConstrainedOPF(
        ref::NetworkReference,
        solver::MathProgBase.AbstractMathProgSolver
    )
    model = JuMPChance.ChanceModel(solver=solver)
    JuMP.@variable(model, pmin(ref,i) <= p[i in 1:ref.ngen] <= pmax(ref,i), start=pstart(ref,i))
    JuMP.@variable(model,                α[i in 1:ref.ngen] >= 0)
    JuMPChance.@indepnormal(model,       ω[j in 1:ref.nuncertain],
        mean=ref.refω[j]["mean"], var=ref.refω[j]["std"]^2
    )
    f = lineflow(ref, p, α, ω)
    JuMP.@constraints model begin
        sum(α[i] for i in 1:ref.ngen) == 1
        sum(α[g] for g in ref.busgens[ref.r]) == 0
    end
    for i in 1:ref.ngen
        JuMP.@constraint(model, p[i] - sum(ω)*α[i] <= pmax(ref,i), with_probability=ref.bus_prob)
        JuMP.@constraint(model, p[i] - sum(ω)*α[i] >= pmin(ref,i), with_probability=ref.bus_prob)
    end
    for l in 1:ref.nline
        JuMP.@constraint(model, f[l] <= rate(ref, l), with_probability=ref.line_prob)
        JuMP.@constraint(model, f[l] >= -rate(ref, l), with_probability=ref.line_prob)
    end
    JuMP.@objective(model, Min, sum(cost(ref,i,1)*p[i] + cost(ref,i,2)*p[i] + cost(ref,i,3) for i in 1:ref.ngen))
    ChanceConstrainedOPF(model,p,α)
end

ChanceConstrainedOPF(filename::String; kwargs...) =
    ChanceConstrainedOPF(PM.build_ref(PM.parse_file(filename)); kwargs...)

ChanceConstrainedOPF(ref::Dict{Symbol,Any}; kwargs...) =
    ChanceConstrainedOPF(NetworkReference(ref); kwargs...)

mutable struct FullChanceConstrainedOPF
    model::JuMP.Model
    p
    α
end

function FullChanceConstrainedOPF(
        ref::NetworkReference,
        solver::MathProgBase.AbstractMathProgSolver
    )
    model = JuMPChance.ChanceModel(solver=solver)
    JuMP.@variable(model, pmin(ref,i) <= p[i in 1:ref.ngen] <= pmax(ref,i), start=pstart(ref,i))
    JuMP.@variable(model,                α[i in 1:ref.ngen, j in 1:ref.nuncertain] >= 0)
    JuMPChance.@indepnormal(model,       ω[j in 1:ref.nuncertain],
        mean=ref.refω[j]["mean"], var=ref.refω[j]["std"]^2
    )
    f = fulllineflow(ref, p, α, ω)
    JuMP.@constraints model begin
        [j in 1:ref.nuncertain], sum(α[i,j] for i in 1:ref.ngen) == 1
    end
    for i in 1:ref.ngen
        JuMP.@constraint(model, p[i] - sum(α[i,j]*ω[j] for j in 1:ref.nuncertain) <= pmax(ref,i), with_probability=ref.bus_prob)
        JuMP.@constraint(model, p[i] - sum(α[i,j]*ω[j] for j in 1:ref.nuncertain) >= pmin(ref,i), with_probability=ref.bus_prob)
    end
    for l in 1:ref.nline
        JuMP.@constraint(model, f[l] <= rate(ref, l), with_probability=ref.line_prob)
        JuMP.@constraint(model, f[l] >= -rate(ref, l), with_probability=ref.line_prob)
    end
    JuMP.@objective(model, Min, sum(cost(ref,i,1)*p[i] + cost(ref,i,2)*p[i] + cost(ref,i,3) for i in 1:ref.ngen))
    FullChanceConstrainedOPF(model,p,α)
end

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
    JuMP.@variable(model,  pmin(ref,i) <= p[i in 1:ref.ngen] <= pmax(ref,i),  start = pstart(ref,i))
    JuMP.@variable(model,                 ω[i in 1:ref.nuncertain])
    f = lineflow(ref, p, ω)
    JuMP.@constraints model begin
        f .<= [rate(ref,l) for l in 1:ref.nline]
        f .>= [-rate(ref,l) for l in 1:ref.nline]
        0 == sum(sum(p[g] for g in ref.busgens[i]) +
                 sum(ω[j] for j in 1:ref.nuncertain if ref.refω[j]["bus"] == i) -
                 ref.bus[i]["pd"] - ref.bus[i]["gs"]
                 for i in 1:ref.nbus)
    end
    JuMP.@objective(model, Min, sum(cost(ref,i,1)*p[i] + cost(ref,i,2)*p[i] + cost(ref,i,3) for i in 1:ref.ngen))
    SingleScenarioOPF(model,p,ω)
end

SingleScenarioOPF(filename::String; kwargs...) =
    SingleScenarioOPF(PM.build_ref(PM.parse_file(filename)); kwargs...)

SingleScenarioOPF(ref::Dict{Symbol,Any}; kwargs...) =
    SingleScenarioOPF(NetworkReference(ref); kwargs...)

function get_opf_solution(opf::SingleScenarioOPF, ω)
    for i in eachindex(ω); JuMP.fix(opf.ω[i], ω[i]) end
    JuMP.solve(opf.model)
end
