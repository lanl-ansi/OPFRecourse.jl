mutable struct OPFScenarios
    ref::NetworkReference
    opf::SingleScenarioOPF
    ω::Distributions.Distribution
    scenarios::Matrix{Float64} # nscenarios x nuncertain
    solutions::Matrix{Float64} # nscenarios x ngens
    cbases::Vector{Vector{Symbol}} # nbasis x num_vars
    rbases::Vector{Vector{Symbol}} # nbasis x num_constrs
    whichbasis::Matrix{Int} # nscenarios x 2 (indicating index of cbasis and rbasis)
    whichscenario::Dict{Tuple{Int,Int},Vector{Int}}
end

function OPFScenarios(ref::NetworkReference, m::SingleScenarioOPF; nsamples::Int=1000)
    nonzeroindices = (1:length(ref.stdω))[ref.stdω .> 1e-5]
    ω = Distributions.MvNormal(
        zeros(length(nonzeroindices)),
        diagm(ref.stdω[nonzeroindices])^2
    )
    ωsamples = zeros(ref.nbus, nsamples)
    ωsamples[nonzeroindices, :] = rand(ω, nsamples)
    status = Array{Symbol}(nsamples);
    soln_p = zeros(nsamples, ref.ngen);
    JuMP.build(m.model);
    nv = Gurobi.num_vars(m.model.internalModel.inner)
    nc = Gurobi.num_constrs(m.model.internalModel.inner)
    cbases = Dict{Vector{Symbol},Vector{Int}}()
    rbases = Dict{Vector{Symbol},Vector{Int}}()
    noptimal = 0

    for i in 1:nsamples
        for j in eachindex(m.ω); JuMP.fix(m.ω[j], ωsamples[j,i]) end
        status[i] = JuMP.solve(m.model)

        if status[i] == :Optimal
            # if the scenario was feasible, keep track of basis
            soln_p[i,:] = JuMP.getvalue(m.p)
            noptimal += 1
            cbasis, rbasis = MathProgBase.getbasis(m.model.internalModel)
            cbases[cbasis] = get(cbases, cbasis, Int[])
            rbases[rbasis] = get(rbases, rbasis, Int[])
            push!(cbases[cbasis], noptimal)
            push!(rbases[rbasis], noptimal)
        end
    end
    @assert noptimal == sum(status .== :Optimal)

    sample_p = soln_p[status .== :Optimal, :]
    sample_ω = ωsamples[:,status .== :Optimal]'

    colbases = map(collect,keys(cbases))
    rowbases = map(collect,keys(rbases))
    whichcol = Dict(zip(colbases,1:length(colbases)))
    whichrow = Dict(zip(rowbases,1:length(rowbases)))
    whichbasis = zeros(Int, noptimal, 2)
    for k in keys(cbases); whichbasis[cbases[k],1] = whichcol[collect(k)] end
    for k in keys(rbases); whichbasis[rbases[k],2] = whichrow[collect(k)] end

    whichscenario = Dict{Tuple{Int,Int},Vector{Int}}()
    for i in 1:noptimal
        basiskey = (whichbasis[i,1], whichbasis[i,2])
        whichscenario[basiskey] = get(whichscenario, basiskey, Int[])
        push!(whichscenario[basiskey], i)
    end

    OPFScenarios(ref, m, ω, sample_ω, sample_p, colbases, rowbases, whichbasis, whichscenario)
end
