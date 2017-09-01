struct Bus
    pd::Float64
    gs::Float64
    gens::Vector{Int}
end

struct Gen
    bus::Int
    pmin::Float64
    pmax::Float64
    pstart::Float64
    cost::Vector{Float64}
end

struct Line
    rate::Float64
    frombus::Int
    tobus::Int
    β::Float64
end

mutable struct NetworkReference
    ref::Dict{Symbol,Any}
    
    nbus::Int
    ngen::Int
    nline::Int

    r::Int # reference bus

    bus::Vector{Bus}
    gen::Vector{Gen}
    line::Vector{Line}

    originalindices::Dict{Symbol,Vector{Int}}

    B::Matrix{Float64}
    π::Matrix{Float64}

    stdω::Vector{Float64}
    line_prob::Float64
    bus_prob::Float64
end

function NetworkReference(ref::Dict{Symbol,Any};
        line_prob::Float64 = 0.9, bus_prob::Float64 = 0.9,
        σscaling::Float64 = 0.05
    )
    function generateindices(d::Dict)
        originalindices = sort(collect(keys(d))); nindices = length(originalindices)
        reverseindices = Dict(zip(originalindices,1:nindices))
        nindices, originalindices, reverseindices
    end
    ngen, genindices, gen_index = generateindices(ref[:gen])
    gen = [Gen(
        ref[:gen][genindices[i]]["gen_bus"],
        ref[:gen][genindices[i]]["pmin"],
        ref[:gen][genindices[i]]["pmax"],
        PowerModels.getstart(ref[:gen],genindices[i],"pg_start"),
        ref[:gen][genindices[i]]["cost"]
    ) for i in 1:ngen]
    nbus, busindices, bus_index = generateindices(ref[:bus])
    bus = [Bus(
        ref[:bus][busindices[i]]["pd"],
        ref[:bus][busindices[i]]["gs"],
        [gen_index[g] for g in ref[:bus_gens][busindices[i]]]
    ) for i in 1:nbus]
    nline, lineindices, line_index = generateindices(ref[:branch])
    line = [Line(
        ref[:branch][lineindices[l]]["rate_a"],
        bus_index[ref[:branch][lineindices[l]]["f_bus"]],
        bus_index[ref[:branch][lineindices[l]]["t_bus"]],
        1 / ref[:branch][lineindices[l]]["br_x"]
    ) for l in 1:nline]
    originalindices = Dict(:bus => busindices, :gen => genindices, :line => lineindices)

    length(ref[:ref_buses]) > 1 && warn("Using only the first reference bus")
    r = bus_index[ref[:ref_buses][first(keys(ref[:ref_buses]))]["bus_i"]]

    nonref_indices = [b for b in 1:nbus if b != r]
    B = admittancematrix(ref, bus_index)
    π = zeros(nbus,nbus)
    π[nonref_indices, nonref_indices] = inv(B[nonref_indices, nonref_indices])

    stdω = [σscaling*(0.01 + ref[:bus][busindices[i]]["pd"]) for i in 1:nbus]

    NetworkReference(ref, nbus, ngen, nline, r, bus, gen, line,
        originalindices, B, π, stdω, line_prob, bus_prob)
end

NetworkReference(filename::String; kwargs...) = NetworkReference(
    PowerModels.build_ref(PowerModels.parse_file(filename)); kwargs...
)

"The total cost of power generation plan `p`."
cost(ref::NetworkReference, p::Vector) = sum(
    ref.gen[i].cost[1]*p[i] + ref.gen[i].cost[2]*p[i] + ref.gen[i].cost[3]
    for i in 1:ref.ngen
)

"Computes the bus admittancematrix B, with indices given by `bus_index::Dict`"
function admittancematrix(ref::Dict{Symbol,Any}, bus_index::Dict{Int,Int})
    nbus = length(ref[:bus])
    B = zeros(nbus,nbus)
    for (i,branch) in ref[:branch]
        f_bus = bus_index[ref[:bus][branch["f_bus"]]["index"]]
        t_bus = bus_index[ref[:bus][branch["t_bus"]]["index"]]
        B[f_bus, t_bus] += (-branch["br_x"]/(branch["br_x"]^2+branch["br_r"]^2))
        B[t_bus, f_bus] += (-branch["br_x"]/(branch["br_x"]^2+branch["br_r"]^2))
        B[f_bus, f_bus] += (branch["br_x"]/(branch["br_x"]^2+branch["br_r"]^2))
        B[t_bus, t_bus] += (branch["br_x"]/(branch["br_x"]^2+branch["br_r"]^2))
    end
    B
end

"""
Computes the `i`-th entry of the phase angle `inv(B̃)*busvalue`

where `busvalue(i)` should correspond to the `i`-th entry of `p + μ + ω - d`
"""
θ(ref, busvalue::Function, i) = sum(ref.π[i,j]*busvalue(j) for j in 1:ref.nbus)

"""
Computes the `i`-th entry of the phase angle `inv(B̃)*busvalue`

where `busvalue[i]` should correspond to the `i`-th entry of `p + μ + ω - d`
"""
θ(ref, busvalue::Vector, i) = sum(ref.π[i,j]*busvalue[j] for j in 1:ref.nbus)

"computes B_f*inv(B̃)(p + μ + ω - d) corresponding to line `l`"
function lineflow(ref, p, ω, l)
    function busvalue(i)
        result = ω[i] - ref.bus[i].pd - ref.bus[i].gs
        if !isempty(ref.bus[i].gens)
            result += sum(p[g] for g in ref.bus[i].gens)
        end
        result
    end
    ref.line[l].β*(θ(ref,busvalue,ref.line[l].frombus) - θ(ref,busvalue,ref.line[l].tobus))
end

"computes B_f*inv(B̃)(p - α⋅ω + μ + ω - d) with aggregated affine recourse corresponding to line `l`"
function lineflow(ref, p, α, ω, l)
    function busvalue(i)
        result = ω[i] - ref.bus[i].pd - ref.bus[i].gs
        if !isempty(ref.bus[i].gens)
            result += sum(p[g] - α[g]*sum(ω) for g in ref.bus[i].gens)
        end
        result
    end
    ref.line[l].β*(θ(ref,busvalue,ref.line[l].frombus) - θ(ref,busvalue,ref.line[l].tobus))    
end

"computes B_f*inv(B̃)(p - α⋅ω + μ + ω - d) with full affine recourse corresponding to line `l`"
function fulllineflow(ref, p, α, ω, l)
    function busvalue(i)
        result = ω[i] - ref.bus[i].pd - ref.bus[i].gs
        if !isempty(ref.bus[i].gens)
            result += sum(
                p[g] - sum(α[g,j]*ω[j] for j in 1:ref.nbus)
                for g in ref.bus[i].gens
            )
        end
        result
    end
    ref.line[l].β*(θ(ref,busvalue,ref.line[l].frombus) - θ(ref,busvalue,ref.line[l].tobus))
end
