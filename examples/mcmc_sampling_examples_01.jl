
using SamplingHS3


#f_const = (x::Vector{Float64}) -> 1.0

# "inverted" egg crate function for defining likelihood:
#f_ec = (x::Vector{Float64}) -> ( cos(x[1])^2 + cos(x[2])^2 ) / ( x[1]*x[1] + x[2]*x[2] )^(1/3)
f_ec = (x::Vector{Float64}) -> ( sin(x[1])^2 + sin(x[2])^2 ) / ( 0.1 + x[1]*x[1] + x[2]*x[2] )^(1/2.5)
f_double_ec = (x::Vector{Float64}) -> max( f_ec( x+[pi/4;pi/4] ) , f_ec( x-[pi/4;pi/4] ) )

#ndgrid(xlim,ylim,xres,yres) = ( [ xlim[1]+ (0.5+i)*((xlim[2]-xlim[1])/(xres+1)) for i=1:xres, j=1:yres ] , [ ylim[1]+ (0.5+j)*((ylim[2]-ylim[1])/(yres+1)) for i=1:xres, j=1:yres ] )
xres = 100; yres = 100;
xygrids = SamplingHS3.xygrid([-10,10,],[-10,10],xres,yres)
zgrid   = NaN*ones(xres,yres)
for zi=1:length(zgrid); zgrid[zi] = f_ec([xygrids[1][zi];xygrids[2][zi]] ); end

using Plots
plotly()
Plots.heatmap(zgrid)

#scatter( xygrids[1][:] , xygrids[2][:] , ma = 0.75 , zcolor=zgrid[:] , m=2, msw=0 , c=:lightrainbow) )

xlim = [-20.,20.]; ylim = [-20.,20.]
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
conf_only_mh = SamplingHS3.DecorrelationConfig( 100000 , 1.0*eye(2) , 0.0 , 1.0 )
f_density    = f_likelihood
# this is quite slow.. most time is spent for the generation of the multivar. normal random numbers..
tic()
mcmc_result = SamplingHS3.run_mcmc( x0 , f_density(x0) , f_density , (x) -> 0 , 1.0 , bounding_box_G  , bounding_box_h , conf_only_mh )
toc()


# pure mcmc sampling using hr:
conf_only_hr = SamplingHS3.DecorrelationConfig( 10000 , 1.0*eye(2) , 1.0 , 0.0 )
f_density    = f_likelihood
# this is quite slow.. most time is spent for the generation of the multivar. normal random numbers..
tic()
mcmc_result_hr = SamplingHS3.run_mcmc( x0 , f_density(x0) , f_density , (x) -> 0 , 1.0 , bounding_box_G  , bounding_box_h , conf_only_hr )
toc()


using Plots
#plotly()
if(false)
    gr()
    Plots.plot(linspace(xlim[1],xlim[2],40),linspace(ylim[1],ylim[2],40),f_likelihood_2p, st = [:contourf])
    Plots.scatter!( mcmc_result.dd.xx[1,:] , mcmc_result.dd.xx[2,:] , ma = 0.75 , m=2, msw=0 )
    Plots.plot!( mcmc_result.dd.xx[1,:] , mcmc_result.dd.xx[2,:] , la = 0.1)
end
using Plots, GR
Plots.gr()
Plots.plot(layout=(1,2))
z_sampled = SamplingHS3.bin2d( mcmc_result_hr.dd.xx , 40 , 40 )
Plots.plot!( linspace(xlim[1],xlim[2],40), linspace(ylim[1],ylim[2],40) , f_likelihood_2p , subplot=1 , st = [:contourf])
Plots.plot!( z_sampled[2], [z_sampled[3]] , z_sampled[1][:]  , subplot=2 , st = [:contourf] , xlim=[-20;20] , ylim=[-20;20])
#Plots.plot!( linspace(xlim[1],xlim[2],40), linspace(ylim[1],ylim[2],40) , z_sampled[:]    , subplot=2 , st = [:contourf])
abc = Plots.savefig("C:\\Temp\\two_ec_densities.png")

# Lets define a restriction to a torus around the center
f_torus = (x) -> (sqrt( sum(x.^2) ) - 5)^2
f_torus_threshold = 1.0
x0 = [5.0;0.0]
mcmc_result_torus = SamplingHS3.run_mcmc( x0 , f_density(x0) , f_density , f_torus , f_torus_threshold , bounding_box_G  , bounding_box_h , conf_only_mh )

# Plot result binned:
using Plots, GR
gr()
Plots.plot(layout=(1,2))
z_sampled_tor = SamplingHS3.bin2d( mcmc_result_torus.dd.xx , 40 , 40 )
Plots.plot!( linspace(xlim[1],xlim[2],40), linspace(ylim[1],ylim[2],40) , f_likelihood_2p , subplot=1 , st = [:contourf])
Plots.plot!( linspace(xlim[1],xlim[2],40), linspace(ylim[1],ylim[2],40) , (x,y) -> 0. , subplot=2 , st = [:contourf])
Plots.plot!( z_sampled_tor[2], z_sampled_tor[3] , z_sampled_tor[1][:]    , subplot=2 , st = [:contourf] , xlim=[-20;20] , ylim=[-20;20])
abc = Plots.savefig("C:\\Temp\\two_ec_densities_torus.png")
