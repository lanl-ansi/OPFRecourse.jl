mutable struct NetworkReference
    ref::Dict{Symbol,Any}
    nbus::Int
    ngen::Int
    nline::Int
    r::Int # reference bus
    stdω::Vector{Float64}
    bus::Dict{Int,Dict{String,Any}}
    busgens::Dict{Int,Vector{Any}}
    B::Matrix{Float64}
    π::Matrix{Float64}
    line_prob::Float64
    bus_prob::Float64
end

function NetworkReference(ref::Dict{Symbol,Any};
        line_prob::Float64 = 0.9, bus_prob::Float64 = 0.9,
        σscaling::Float64 = 0.05
    )
    @assert validate_ref(ref)
    nbus, r = length(ref[:bus]), ref[:ref_buses][1]["bus_i"]
    nonref_indices = [b for b in 1:nbus if b != r]
    B = admittancematrix(ref)
    stdω = [σscaling*(0.01 + ref[:bus][i]["pd"]) for i in 1:nbus]
    π = zeros(nbus,nbus)
    π[nonref_indices, nonref_indices] = inv(B[nonref_indices, nonref_indices])
    NetworkReference(ref, nbus, length(ref[:gen]), length(ref[:branch]),
        r, stdω, ref[:bus], ref[:bus_gens], B, π, line_prob, bus_prob)
end

NetworkReference(filename::String; kwargs...) = NetworkReference(
    PowerModels.build_ref(PowerModels.parse_file(filename)); kwargs...
)

rate(ref::NetworkReference, l) = ref.ref[:branch][l]["rate_a"]

frombus(ref::NetworkReference, l) = ref.ref[:branch][l]["f_bus"]

tobus(ref::NetworkReference, l) = ref.ref[:branch][l]["t_bus"]

β(ref::NetworkReference, l) = 1 / ref.ref[:branch][l]["br_x"]

pmin(ref::NetworkReference, i) = ref.ref[:gen][i]["pmin"]

pmax(ref::NetworkReference, i) = ref.ref[:gen][i]["pmax"]

pstart(ref::NetworkReference, i) = PowerModels.getstart(ref.ref[:gen],i,"pg_start")

cost(ref::NetworkReference, i, c) = ref.ref[:gen][i]["cost"][c]

cost(ref::NetworkReference, p::Vector) = sum(
    cost(ref,i,1)*p[i] + cost(ref,i,2)*p[i] + cost(ref,i,3) for i in 1:ref.ngen
)

function admittancematrix(ref::Dict{Symbol,Any})
    nbus = length(ref[:bus])
    B = zeros(nbus,nbus)
    for (i,branch) in ref[:branch]
        f_bus = ref[:bus][branch["f_bus"]]["index"]
        t_bus = ref[:bus][branch["t_bus"]]["index"]
        B[f_bus, t_bus] += (-branch["br_x"]/(branch["br_x"]^2+branch["br_r"]^2))
        B[t_bus, f_bus] += (-branch["br_x"]/(branch["br_x"]^2+branch["br_r"]^2))
        B[f_bus, f_bus] += (branch["br_x"]/(branch["br_x"]^2+branch["br_r"]^2))
        B[t_bus, t_bus] += (branch["br_x"]/(branch["br_x"]^2+branch["br_r"]^2))
    end
    B
end

function validate_ref(ref::Dict{Symbol,Any})
    nbus, ngen, nline = length(ref[:bus]), length(ref[:gen]), length(ref[:branch])
    sort(collect(keys(ref[:branch]))) == collect(1:nline) || return false
    sort(collect(keys(ref[:gen]))) == collect(1:ngen) || return false
    sort(collect(keys(ref[:bus]))) == collect(1:nbus) || return false
    sort([l for (l,i,j) in ref[:arcs_from]]) == collect(1:nline) || return false
    length(ref[:ref_buses]) == 1 || return false
    for i in keys(ref[:branch])
        i == ref[:branch][i]["index"] || return false
    end
    for i in keys(ref[:bus])
        i == ref[:bus][i]["index"] || return false
    end
    return true
end
