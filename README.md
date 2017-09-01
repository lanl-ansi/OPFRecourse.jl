# OPFRecourse

[![Build Status](https://travis-ci.org/yeesian/OPFRecourse.jl.svg?branch=master)](https://travis-ci.org/yeesian/OPFRecourse.jl)

[![Coverage Status](https://coveralls.io/repos/yeesian/OPFRecourse.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/yeesian/OPFRecourse.jl?branch=master)

[![codecov.io](http://codecov.io/github/yeesian/OPFRecourse.jl/coverage.svg?branch=master)](http://codecov.io/github/yeesian/OPFRecourse.jl?branch=master)

Package containing supplementary data and code for the paper "Optimal Recourse Design for Stochastic Optimal Power Flow".

## Installation Instructions
The CCOPF models were implemented using the [JuMPChance](https://github.com/mlubin/JuMPChance.jl) extension to [JuMP](https://github.com/JuliaOpt/JuMP.jl) in the [Julia programming language](https://julialang.org/downloads/). Additionally, we used [Gurobi](https://www.gurobi.com) 7.5 in our numerical experiments. Gurobi is a commercial solver which must be installed and licensed separately (one may easily use a different solver if Gurobi is not available, see the JuMP documentation).

You should first install the following packages by running `Pkg.add(...)` in a Julia REPL:

- Distributions
- JuMP
- MathProgBase
- Gurobi
- PowerModels
- JuMPChance

Finally, you can install the package by running `Pkg.clone("https://github.com/yeesian/OPFRecourse.jl.git")`, or manually downloading this github repository into your julia package directory (you can run `Pkg.dir()` in a Julia REPL to locate it).

To test if it is successfully installed, run `Pkg.test("OPFRecourse")` in a Julia REPL.
