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
    nbus, busindices, bus_index = generateindices(ref[:bus])
    bus = [Bus(
        ref[:bus][busindices[i]]["pd"],
        ref[:bus][busindices[i]]["gs"],
        Int.(ref[:bus_gens][busindices[i]])
    ) for i in 1:nbus]
    ngen, genindices, gen_index = generateindices(ref[:gen])
    gen = [Gen(
        ref[:gen][genindices[i]]["gen_bus"],
        ref[:gen][genindices[i]]["pmin"],
        ref[:gen][genindices[i]]["pmax"],
        PowerModels.getstart(ref[:gen],genindices[i],"pg_start"),
        ref[:gen][genindices[i]]["cost"]
    ) for i in 1:ngen]
    nline, lineindices, line_index = generateindices(ref[:branch])
    line = [Line(
        ref[:branch][lineindices[l]]["rate_a"],
        ref[:branch][lineindices[l]]["f_bus"],
        ref[:branch][lineindices[l]]["t_bus"],
        1 / ref[:branch][lineindices[l]]["br_x"]
    ) for l in 1:nline]
    originalindices = Dict(:bus => busindices, :gen => genindices, :line => lineindices)

    length(ref[:ref_buses]) > 1 && warn("Using only the first reference bus")
    r = bus_index[ref[:ref_buses][1]["bus_i"]]

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

cost(ref::NetworkReference, p::Vector) = sum(
    ref.gen[i].cost[1]*p[i] + ref.gen[i].cost[2]*p[i] + ref.gen[i].cost[3]
    for i in 1:ref.ngen
)

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
