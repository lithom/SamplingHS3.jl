


using SamplingHS3

# inverted egg crate function
f_ec = (x::Vector{Float64}) -> ( sin(x[1])^2 + sin(x[2])^2 ) / ( 0.1 + x[1]*x[1] + x[2]*x[2] )^(1/2.5)


# non-adaptive mcmc sampling using hr and mh:
conf_only_hr = SamplingHS3.DecorrelationConfig( 10000 , 1.0*eye(2) , 0.5 , 0.5 )
f_density    = f_ec
# this is quite slow.. most time is spent for the generation of the multivar. normal random numbers..
tic()
mcmc_result_01 = SamplingHS3.run_mcmc( x0 , f_density(x0) , f_density , (x) -> 0 , 1.0 , bounding_box_G  , bounding_box_h , conf_only_hr )
toc()

# run with adaptive mcmc
tic()
config_mcmc = SamplingHS3.MCMCConfig()
mcmc_result_adaptive_01 = SamplingHS3.run_adaptive_mcmc( 10000 , x0 , f_density(x0) , f_density , (x) -> 0 , 1.0 , bounding_box_G  , bounding_box_h , config_mcmc )
toc()

# Plot results
