

using SamplingHS3

# inverted egg crate function
f_ec = (x::Vector{Float64}) -> ( sin(x[1])^2 + sin(x[2])^2 ) / ( 0.1 + x[1]*x[1] + x[2]*x[2] )^(1/2.5)

# polytope
xlim = [-20.,20.]; ylim = [-20.,20.]
bounding_box_G =[eye(2);-eye(2)]
bounding_box_h =[-xlim[1];xlim[2];-ylim[1];ylim[2]]

x0 = [5.0;5.0]

# non-adaptive mcmc sampling using hr and mh:
conf_only_hr = SamplingHS3.DecorrelationConfig( 10000 , 1.0*eye(2) , 0.5 , 0.5 )
f_density    = f_ec
# this is quite slow.. most time is spent for the generation of the multivar. normal random numbers..
tic()
mcmc_result_01 = SamplingHS3.run_mcmc( x0 , f_density(x0) , f_density , (x) -> 0 , 1.0 , bounding_box_G  , bounding_box_h , conf_only_hr )
toc()

# run with adaptive mcmc
tic()
config_mcmc = SamplingHS3.MCMCConfig() # default config
mcmc_result_adaptive_01 = SamplingHS3.run_adaptive_mcmc( 10000 , x0 , f_density(x0) , f_density , (x) -> 0 , 1.0 , bounding_box_G  , bounding_box_h , config_mcmc )
toc()

# Plot results

using Plots, GR
Plots.gr()
Plots.plot(layout=(2,3))
z_sampled_1 = SamplingHS3.bin2d( mcmc_result_01.dd.xx , 40 , 40 )
z_sampled_2 = SamplingHS3.bin2d( mcmc_result_adaptive_01.dd.xx , 40 , 40 )
layout_fancy = @layout [a{0.33w} b{0.33w} c{0.33w}; d{0.5h}]
Plots.plot(layout=layout_fancy)
Plots.plot!( linspace(xlim[1],xlim[2],40), linspace(ylim[1],ylim[2],40) , f_likelihood_2p , subplot=1 , st = [:contourf])
Plots.plot!( z_sampled_1[2], [z_sampled_1[3]] , z_sampled_1[1][:]  , subplot=2 , st = [:contourf] , xlim=[-20;20] , ylim=[-20;20])
Plots.plot!( z_sampled_2[2], [z_sampled_2[3]] , z_sampled_2[1][:]  , subplot=3 , st = [:contourf] , xlim=[-20;20] , ylim=[-20;20])
# evolution of estimated covariance matrix:
ecov_steps = sort(collect(keys(mcmc_result_adaptive_01.ecov_regularized)))
get_ecov_entries = (i::Int) -> map( x -> mcmc_result_adaptive_01.ecov_regularized[x][i] , sort(ecov_steps) )
Plots.plot!( ecov_steps , get_ecov_entries(1) ,  subplot=4 , lab="a")
Plots.plot!( ecov_steps , get_ecov_entries(2) ,  subplot=4 , lab="b/c")
#Plots.plot!( ecov_steps , get_ecov_entries(3) ,  subplot=4 , lab="b/c")
Plots.plot!( ecov_steps , get_ecov_entries(4) ,  subplot=4 , lab="d")
Plots.title!( "Estimated Covariance Matrix" , subplot=4 )

#Plots.plot!( linspace(xlim[1],xlim[2],40), linspace(ylim[1],ylim[2],40) , z_sampled[:]    , subplot=2 , st = [:contourf])
abc = Plots.savefig("C:\\Temp\\non_adaptive_vs_adaptive_01.png" )
