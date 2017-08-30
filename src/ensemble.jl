mutable struct EnsembleRecourse
    ref::NetworkReference
    baseline
    recoursef::Vector
end

function get_opf_solution(er::EnsembleRecourse, ω)
    incumbent_p = get_opf_solution(er.baseline, ω)
    incumbent_cost = cost(er.ref, incumbent_p)
    feasible_p = nviolations(er.ref, incumbent_p, ω) == 0
    for rf in er.recoursef
        p = get_opf_solution(rf, ω)
        if nviolations(er.ref, p, ω) == 0 # feasible solution
            curr_cost = cost(er.ref, p)
            if !feasible_p
                # prioritize feasibility
                incumbent_cost = curr_cost
                incumbent_p = p
                feasible_p = true
            elseif curr_cost < incumbent_cost
                # use the solution with lower cost
                incumbent_cost = curr_cost
                incumbent_p = p
            end
        end
    end
    incumbent_p
end
