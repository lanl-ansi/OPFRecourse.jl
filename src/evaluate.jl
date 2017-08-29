
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
        newp[i] = max(pmin(ref,i), min(pmax(ref,i), p[i]))
    end
    newp[ref.ngen] = p[ref.ngen] + totalpower - sum(newp)
    newp
end
