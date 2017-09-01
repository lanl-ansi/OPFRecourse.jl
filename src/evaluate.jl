
"""
Projects the values in p to within the generation limits, with a slack generator
to account for total power balance.

Assumes that the last entry of p is the slack generator for now.
"""
function reproject(ref::NetworkReference, p)
    @assert length(p) > 1
    @assert ref.ngen == length(p)
    totalpower = sum(p[i] for i in 1:(ref.ngen-1))
    newp = zeros(Float64, ref.ngen)
    for i in 1:(ref.ngen-1)
        newp[i] = max(ref.gen[i].pmin, min(ref.gen[i].pmax, p[i]))
    end
    newp[ref.ngen] = p[ref.ngen] + totalpower - sum(newp)
    newp
end

"""Returns the number of transmission lines with flows exceeding their rates

Evaluated under generation plan `p` and forecast deviations `ω`.
"""
function ntransmissionviolations(ref::NetworkReference, p, ω::Vector; atol::Float64 = 1e-5)
    sum(abs(lineflow(ref, p, ω, l)) > ref.line[l].rate + atol for l in 1:ref.nline)
end

"""Returns the total number of transmission lines with flows exceeding their rates

Evaluated under generation plan `p` evaluated over all forecast deviations `ω[i,:]`.
"""
function ntransmissionviolations(ref::NetworkReference, p, ω::Matrix; atol::Float64 = 1e-5)
    sum(ntransmissionviolations(ref, p, ω[i,:], atol=atol) for i in 1:size(ω,1))
end

"Returns the number of generation limits that were violated under generation plan `p`"
function ngenerationviolations(ref::NetworkReference, p; atol::Float64 = 1e-5)
    sum(ref.gen[i].pmin - atol > p[i] for i in 1:ref.ngen) +
    sum(ref.gen[i].pmax + atol < p[i] for i in 1:ref.ngen)
end

"""Returns the total number of constraints that were violated

Evaluated under generation plan `p` and forecast deviations `ω`.
"""
function nviolations(ref::NetworkReference, p, ω; atol::Float64 = 1e-5)
    ngenerationviolations(ref, p) + ntransmissionviolations(ref, p, ω)
end
