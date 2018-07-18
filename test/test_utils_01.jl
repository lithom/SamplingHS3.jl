# First run some mcmc sampling, then test the util functions..

f_ec = (x::Vector{Float64}) -> ( sin(x[1])^2 + sin(x[2])^2 ) / ( 0.1 + x[1]*x[1] + x[2]*x[2] )^(1/2.5)
f_double_ec = (x::Vector{Float64}) -> max( f_ec( x+[pi/4;pi/4] ) , f_ec( x-[pi/4;pi/4] ) )


xlim = [-20.,20.]
ylim = [-20.,20.]
bounding_box_G =[eye(2);-eye(2)]
bounding_box_h =[-xlim[1];xlim[2];-ylim[1];ylim[2]]

# and create likelihood from "inverted egg crate"
f_likelihood = (x::Vector{Float64}) -> f_ec(x)
f_likelihood_2p = (x,y) -> f_ec([x;y])
#f_likelihood = (x::Vector{Float64}) -> f_double_ec(x)
#f_likelihood_2p = (x,y) -> f_double_ec([x;y])

# starting point
x0 = [5.;5.]
# pure mcmc sampling using mh:
conf_only_mh = SamplingHS3.DecorrelationConfig( 20000 , 1.0*eye(2) , 0.5 , 0.5 )
f_density    = f_likelihood
# this is quite slow.. most time is spent for the generation of the multivar. normal random numbers..
tic()
mcmc_result = SamplingHS3.run_mcmc( x0 , f_density(x0) , f_density , (x) -> 0 , 1.0 , bounding_box_G  , bounding_box_h , conf_only_mh )
toc()

# Test utilts:
z_sampled      = SamplingHS3.bin2d( mcmc_result.dd.xx , 40 , 40 )
z_sampled_mean = SamplingHS3.bin2d( mcmc_result.dd.xx , mcmc_result.dd.xx_fdens , 40 , 40 , mean)

# And test interp1_lin:
xx_01  = collect(linspace(0.0,10.0,2000))
yy_01  = sin.(xx_01)
x_test = [1.0;2.0;5.0;5.23]
y_test = SamplingHS3.interp1_lin(xx_01,yy_01, x_test )
@test maximum( abs.( y_test - sin.(x_test) ) ) < 0.001
