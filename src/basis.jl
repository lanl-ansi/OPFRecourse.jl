mutable struct BasisRecourse
    ngen::Int
    fixed::Dict{Int,Float64}
    varying::Vector{Int}
    basiscol::Dict{Int,Int}
    linearterms::Matrix{Float64}
    constant::Vector{Float64}
end

function BasisRecourse(ref::NetworkReference, m::SingleScenarioOPF, cbasis, rbasis)
    br = BasisRecourse(ref.ngen,Dict(),Float64[],Dict(),Matrix{Float64}(0,0),[])

    basic_indices = Int[]
    for i in 1:ref.ngen
        if cbasis[i] == :Basic
            push!(basic_indices, i)
        elseif cbasis[i] == :NonbasicAtLower
            # @assert(m.model.colNames[i] == "p[$i]", m.model.colNames[i])
            br.fixed[i] = m.model.colLower[i]
        elseif cbasis[i] == :NonbasicAtUpper
            # @assert(m.model.colNames[i] == "p[$i]", m.model.colNames[i])
            br.fixed[i] = m.model.colUpper[i]
        else
            error("Unrecognised basis status: $(cbasis[i]) at index $i")
        end
    end
    br.varying = basic_indices
    numbasic = sum(rbasis .!== :Basic)
    basiscol = br.basiscol = Dict(zip(basic_indices,1:numbasic))
    @assert length(basic_indices) == numbasic
    @assert issorted(basic_indices)
    @assert length(m.model.linconstr) == 2*ref.nline + 1
    basis = zeros(Float64, numbasic, numbasic)
    ωmatrix = zeros(Float64, numbasic, ref.nuncertain)
    c = 0
    for i in 1:(2*ref.nline + 1)
        terms = m.model.linconstr[i].terms
        if rbasis[i] !== :Basic
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
                if v.col <= ref.ngen # power generator
                    if cbasis[v.col] == :Basic
                        basis[c,basiscol[v.col]] += coeff
                    else
                        rhs -= coeff*br.fixed[v.col]
                    end
                else
                    @assert v.col <= ref.ngen + ref.nuncertain # uncertainty
                    @assert cbasis[v.col] !== :Basic # they are fixed for each scenario
                    ωmatrix[c,v.col-ref.ngen] -= coeff
                end
            end
            push!(br.constant, rhs)
        end
    end
    @assert c == numbasic == length(br.varying) == length(br.constant)
    br.linearterms = inv(basis)*ωmatrix
    br.constant = inv(basis)*br.constant
    br
end

function get_opf_solution(br::BasisRecourse, ω)
    basic_values = br.linearterms*ω + br.constant
    soln = zeros(Float64, br.ngen)
    for i in keys(br.fixed); soln[i] = br.fixed[i] end # nonbasic components
    for i in br.varying; soln[i] = basic_values[br.basiscol[i]] end # basic components
    soln
end
