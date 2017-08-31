const PM = PowerModels

mutable struct NetworkReference
    ref::Dict{Symbol,Any}
    nbus::Int
    ngen::Int
    nline::Int
    nuncertain::Int
    r::Int # reference bus
    refω::Dict{Int,Dict{String,Any}}
    bus::Dict{Int,Dict{String,Any}}
    busgens::Dict{Int,Vector{Any}}
    busarcs::Dict{Int,Vector{Any}}
    arcsfrom::Vector{Tuple{Int,Int,Int}}
    nonref_indices::Vector{Int}
    sqrtΣ::AbstractMatrix{Float64}
    B::Matrix{Float64}
    π::Matrix{Float64}
    varΩ::Float64
    line_prob::Float64
    bus_prob::Float64
end

function NetworkReference(ref::Dict{Symbol,Any};
        line_prob::Float64 = 0.9, bus_prob::Float64 = 0.9
    )
    @assert validate_ref(ref)
    nbus, r, refω = length(ref[:bus]), ref[:ref_buses][1]["bus_i"], ref[:uncertainty]
    nonref_indices = [b for b in 1:nbus if b != r]
    B = admittancematrix(ref); nuncertain = length(refω)
    sqrtΣ = [ref[:corr][i]["col_$j"]*refω[i]["std"]*refω[j]["std"]
            for i in 1:nuncertain, j in 1:nuncertain] ^ 0.5
    varΩ = sum(ref[:corr][i]["col_$j"]*refω[i]["std"]*refω[j]["std"]
               for i in 1:nuncertain, j in 1:nuncertain)
    π = zeros(nbus,nbus)
    π[nonref_indices, nonref_indices] = inv(B[nonref_indices, nonref_indices])
    NetworkReference(ref, nbus, length(ref[:gen]), length(ref[:branch]), nuncertain,
        r, refω, ref[:bus], ref[:bus_gens], ref[:bus_arcs], ref[:arcs_from],
        nonref_indices, sqrtΣ, B, π, varΩ, line_prob, bus_prob)
end

rate(ref::NetworkReference, l) = ref.ref[:branch][l]["rate_a"]

frombus(ref::NetworkReference, l) = ref.ref[:branch][l]["f_bus"]

tobus(ref::NetworkReference, l) = ref.ref[:branch][l]["t_bus"]

β(ref::NetworkReference, l) = 1 / ref.ref[:branch][l]["br_x"]

pmin(ref::NetworkReference, i) = ref.ref[:gen][i]["pmin"]

pmax(ref::NetworkReference, i) = ref.ref[:gen][i]["pmax"]

pstart(ref::NetworkReference, i) = PM.getstart(ref.ref[:gen],i,"pg_start")

θstart(ref::NetworkReference, i) = PM.getstart(ref.ref[:bus],i,"t_start")

fstart(ref::NetworkReference, l) = PM.getstart(ref.ref[:branch],l,"p_start")

cost(ref::NetworkReference, i, c) = ref.ref[:gen][i]["cost"][c]

cost(ref::NetworkReference, p::Vector) = sum(
    cost(ref,i,1)*p[i] + cost(ref,i,2)*p[i] + cost(ref,i,3) for i in 1:ref.ngen
)

πmatrix(ref::NetworkReference, m) = [ref.π[m, ref.refω[j]["bus"]] for j in 1:ref.nuncertain]

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
    nbus, ngen = length(ref[:bus]), length(ref[:gen])
    # nline, nuncertain = length(ref[:branch]), length(ref[:uncertainty])
    nline = length(ref[:branch])
    # sort(collect(keys(ref[:uncertainty]))) == collect(1:nuncertain) || return false
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
