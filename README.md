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
Plots.plot!( linspace(xlim[1],xlim[2],40), linspace(ylim[1],ylim[2],40) , (x,y) -> f_density([x;y]) , subplot=1 , st = [:contourf])
Plots.plot!( z_sampled[2] , z_sampled[3] , z_sampled[1][:]    , subplot=2 , st = [:contourf])
```

![Inverted egg crate function](https://github.com/lithom/SamplingHS3.jl/blob/master/resources/two_ec_densities.png
)


## MCMC sampling with restriction
Let's restrict the sampling to torus-shaped "shell"s of the function f_cost_torus:
```julia
f_cost_torus = (x) -> (sqrt( sum(x.^2) ) - 5)^2
mcmc_result_2 = SamplingHS3.run_mcmc( x0 , f_density(x0) , f_density , f_cost_torus , 1.0 , bounding_box_G  , bounding_box_h , conf_only_mh )
#Plot
using Plots, GR
gr()
Plots.plot(layout=(1,2))
z_sampled_tor = SamplingHS3.bin2d( mcmc_result_torus.dd.xx , 40 , 40 )
Plots.plot!( linspace(xlim[1],xlim[2],40), linspace(ylim[1],ylim[2],40) , (x,y) -> f_density([x;y]) , subplot=1 , st = [:contourf])
Plots.plot!( linspace(xlim[1],xlim[2],40), linspace(ylim[1],ylim[2],40) , (x,y) -> 0. , subplot=2 , st = [:contourf]) # just black backround
Plots.plot!( z_sampled_tor[2], z_sampled_tor[3] , z_sampled_tor[1][:]    , subplot=2 , st = [:contourf] , xlim=[-20;20] , ylim=[-20;20])
```
![Inverted egg crate function, restricted to torus at distance 5 from center](https://github.com/lithom/SamplingHS3.jl/blob/master/resources/two_ec_densities_torus.png
)

## Adaptive MCMC sampling
Let's use adaptive MCMC sampling for the same function:

```julia
config_mcmc = SamplingHS3.MCMCConfig() # default config
mcmc_result_adaptive_01 = SamplingHS3.run_adaptive_mcmc( 10000 , x0 , f_density(x0) , f_density , (x) -> 0 , 1.0 , bounding_box_G  , bounding_box_h , config_mcmc )
```

Then we can plot and compare the sampling results (full example and plotting code is in "adaptive_mcmc_example_01.jl"):
Center: non-adaptive, Right: adaptive.
![Adaptive vs. non adaptive MCMC sampling](https://github.com/lithom/SamplingHS3.jl/blob/master/resources/non_adaptive_vs_adaptive_01.png
)







Thomas Liphardt, June 2018
