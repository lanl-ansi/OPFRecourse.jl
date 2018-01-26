# OPFRecourse

[![Build Status](https://travis-ci.org/yeesian/OPFRecourse.jl.svg?branch=master)](https://travis-ci.org/yeesian/OPFRecourse.jl) [![Coverage Status](https://coveralls.io/repos/yeesian/OPFRecourse.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/yeesian/OPFRecourse.jl?branch=master) [![codecov.io](http://codecov.io/github/yeesian/OPFRecourse.jl/coverage.svg?branch=master)](http://codecov.io/github/yeesian/OPFRecourse.jl?branch=master)

Package containing supplementary code for the paper "Statistical Learning For DC Optimal Power Flow" by Ng, Misra, A. Roald, and Backhaus.

## Installation Instructions
**OPFRecourse.jl** is not a registered julia package. To install it, you should first install the following dependencies by running `Pkg.add(...)` in a Julia REPL:

- Distributions
- JuMP
- MathProgBase
- PowerModels
- ProgressMeter
- Clp

Secondly, you can install the package by running `Pkg.clone("https://github.com/lanl-ansi/OPFRecourse.jl.git")`, or manually downloading this github repository into your julia package directory (you can run `Pkg.dir()` in a Julia REPL to locate it).

Finally, you should clone a copy of the [PGLIB OPF](https://github.com/power-grid-lib/pglib-opf) repository into your `OPFRecourse/test/data` folder. The current tests are based on [release v17.08](https://github.com/power-grid-lib/pglib-opf/releases/tag/v17.08).

To test if it is successfully installed, run `Pkg.test("OPFRecourse")` in a Julia REPL.

## Quick Example

We begin by importing the relevant packages, and getting the location of a matpower file, e.g.
```julia
using OPFRecourse, JuMP, Clp

data_file = string(Pkg.dir(),"/OPFRecourse/test/data/nesta_case30_ieee.m")
```

### NetworkReference
This package provides a `NetworkReference` object that uses [PowerModels](https://github.com/lanl-ansi/PowerModels.jl) to parse the matpower file, before pulling out the relevant information for the DC-OPF problem. To construct it, run

```julia
ref = OPFRecourse.NetworkReference(data_file, bus_prob = 0.95, line_prob = 0.95, σscaling = 0.05);
```
where optional parameters (i.e. `bus_prob`, `line_prob` and `σscaling`) can be passed in for the Chance Constraint OPF.

### Model Formulations
To construct and solve the DC OPF model, run the following code:

```julia
m = OPFRecourse.SingleScenarioOPF(ref, Clp.ClpSolver());
@time JuMP.solve(m.model)
```

### OPFScenarios
The `ref::NetworkReference` also contains the distributional information for generating scenarios of forecast deviations. In additional, we identify the optimal bases for for each scenario by solving single-scenario OPFs, and store them all within an `OPFScenarios` object:

```
@time scenarios = OPFRecourse.OPFScenarios(ref, m, nsamples = 1000);
```

The set of all unique column bases can be accessed using `scenarios.cbases`, and the set of all unique column bases can be accessed using `scenarios.rbases` respectively. These correspond to the power generation values and line flow constraints respectively, and contain information about system configurations (provided by `scenarios.whichbasis`) in which a subset of them are set to be active. Read the paper for more details.

### Recourse Policies
`SingleScenarioOPF` implements a `get_opf_solution(model, scenario)` method that computes the power generation plan for the given `scenario`. So the following code

```julia
get_opf_solution(m, scenarios.scenarios[1,:])
```

would get the power generation plans corresponding to the first scenario.

### BasisRecourse and EnsembleRecourse
This package also provides a `BasisRecourse` object that implements the recourse policy corresponding to a basis of the OPF problem. You can construct it, e.g. by

```julia
br1 = BasisRecourse(ref, m, scenarios.cbases[2], scenarios.rbases[1])
br2 = BasisRecourse(ref, m, scenarios.cbases[1], scenarios.rbases[2])
```
and be able to call `get_opf_solution(br, scenarios.scenarios[i,:])` just like the other models.

You can also combine them into a `EnsembleRecourse` object that implements an ensemble of policies. For example,

```julia
ensemble = OPFRecourse.EnsembleRecourse(ref, br0, [br1, br2])
```
would construct an ensemble of the policies by `br0`, `br1` and `br2`, falling back to the `br0` policy (because it provides a probabilistic guarantee) if none of them provided a feasible solution.

## Development

Community-driven development and enhancement of OPFRecourse.jl are welcome and encouraged. Please fork this repository and share your contributions to the master with pull requests.

## Citing OPFRecourse
If you find OPFRecourse useful in your work, we kindly request that you cite the following [technical report](http://arxiv.org/abs/1801.07809):
```
@misc{1801.07809,
  author = {Yeesian Ng and Sidhant Misra and Line A. Roald and Scott Backhaus},
  title = {Statistical Learning For DC Optimal Power Flow},
  year = {2018},
  eprint = {arXiv:1801.07809},
  url = {http://arxiv.org/abs/1801.07809}
}
```
We also encourage citing [PowerModels.jl](https://github.com/lanl-ansi/PowerModels.jl) and [PGLIB OPF](https://github.com/power-grid-lib/pglib-opf) if you use the packages in your work.

## License

This code has been developed as part of the Advanced Network Science Initiative at Los Alamos National Laboratory, and is provided under a BSD license as part of the Multi-Infrastructure Control and Optimization Toolkit (MICOT) project, LA-CC-13-108.
