# SamplingHS3
[![Build Status](https://travis-ci.org/lithom/SamplingHS3.jl.svg?branch=master)](https://travis-ci.org/lithom/SamplingHS3.jl)  [![Build status](https://ci.appveyor.com/api/projects/status/le0qla5a7i34adal/branch/master?svg=true)](https://ci.appveyor.com/project/lithom/samplinghs3-jl/branch/master)
  [![Coverage Status](https://coveralls.io/repos/github/lithom/SamplingHS3.jl/badge.svg?branch=master)](https://coveralls.io/github/lithom/SamplingHS3.jl?branch=master)

A collection of simple MCMC samplers, with tools for running nested sampling and methods for parameter space exploration.


# MCMC samplers
The provided MCMC samplers sample according to a specific density function f. Additionally, it is possible to restrict the sampling to "likelihood shells", given by an upper bound on a second function f_cost

Let's define some density function over a 2d polytope (a square), and then use the mcmc samplers:
```julia
# density function
f_ec = (x::Vector{Float64}) -> ( sin(x[1])^2 + sin(x[2])^2 ) / ( 0.1 + x[1]*x[1] + x[2]*x[2] )^(1/2.5)
f_density = f_ec
# polytope:
xlim = [-20.,20.]; ylim = [-20.,20.]
bounding_box_G =[eye(2);-eye(2)]
bounding_box_h =[-xlim[1];xlim[2];-ylim[1];ylim[2]]
```

Sample normally:
```julia
# settings for mcmc sampling:                    steps  cov. matrix   %hr   %mh
conf_only_mh = SamplingHS3.DecorrelationConfig( 100000 , 1.0*eye(2) , 0.0 , 1.0 )
# starting point:
x0 = [5.;5.]
# sample
#                                              these two parameter can be used to restrict the sampled space
#                                                                         ___|___
#                                                                        |       | 
mcmc_result = SamplingHS3.run_mcmc( x0 , f_density(x0) , f_density , (x) -> 0 , 1.0 , bounding_box_G  , bounding_box_h , conf_only_mh )
```

Plot result (binned):
```julia
using Plots, GR
gr()
Plots.plot(layout=(1,2))
z_sampled = SamplingHS3.bin2d( mcmc_result.dd.xx , 40 , 40 )
Plots.plot!( linspace(xlim[1],xlim[2],40), linspace(ylim[1],ylim[2],40) , f_likelihood_2p , subplot=1 , st = [:contourf])
Plots.plot!( linspace(xlim[1],xlim[2],40), linspace(ylim[1],ylim[2],40) , z_sampled[:]    , subplot=2 , st = [:contourf])
```

![Inverted egg crate function](https://github.com/lithom/SamplingHS3.jl/blob/master/res/double_ec_densities.png)







Thomas Liphardt, June 2018
