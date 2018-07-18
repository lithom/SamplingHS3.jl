
using SamplingHS3
using Plots
using GR
# 3d gaussian mixture test..

function plot_3d_view( x_3d::Array{Float64,2} ; xf = [] )
    if(isempty(xf)); xf = ones(size(x_3d)); end
    gr()
    p1 = scatter( x_3d[1,:] , x_3d[2 ,:] , x_3d[3,:] , zcolor = xf , ma = 0.2 , m=3, msw=0 , c=:lightrainbow , leg=:none)
    p2 = scatter( x_3d[1,:] , x_3d[2 ,:]  , zcolor = xf , ma = 0.2 , m=3, msw=0 , c=:lightrainbow , leg=:none)
    p3 = scatter( x_3d[1,:] , x_3d[3 ,:]  , zcolor = xf , ma = 0.2 , m=3, msw=0 , c=:lightrainbow , leg=:none )
    p4 = scatter( x_3d[2,:] , x_3d[3 ,:]  , zcolor = xf , ma = 0.2 , m=3, msw=0 , c=:lightrainbow , leg=:none )

    display(plot(p1,p2,p3,p4,layout=(2,2)))
end

## Create Gausisan mixture
gm3_mix = []
for zi=1:10;
    found_es = false
    vd_pre   = NaN * ones(3,3)
    while(!found_es)
        vd = qr(rand(3,3))[1]
        vd_pre = vd*diagm(rand(3))*vd'
        found_es = issymmetric(vd_pre)
    end
    gpi = SamplingHS3.GaussParam( 8*rand(3) , vd_pre )
    gm3_mix = [ gm3_mix ;  gpi ]
end
gm3 = SamplingHS3.GaussMixture(gm3_mix,0.5+rand(10))
bb3_G = [eye(3);-eye(3)]
bb3_h = [10.;10;10;10;10;10]

## Create cost function (i.e. -log(f(x)) for a pdf given by f(x) )
f_cost_3d = (x) -> -log(SamplingHS3.pdf_gaussmixture_nonnormalized(gm3,x))

## Test Plot
xx3  = SamplingHS3.draw_gaussmixture( gm3 , 100_000 )
plot_3d_view(xx3[:,1:10_000],xf=1:10000)
#scatter(xx3[1,1:2000] , xx3[2 ,1:2000] , xx3[3,1:2000] , ma = 0.1 , m=2, msw=0 , c=:lightrainbow , lab="grad")



## Sample with SamplingHS3.jl
X_init = SamplingHS3.sample_unif_hr( bb3_G , bb3_h , zeros(3) , 400 , 2000 )
ns_result = SamplingHS3.nested_sampling( X_init , f_cost_3d , bb3_G , bb3_h , SamplingHS3.NSConfig(8,0.95,SamplingHS3.StopCriterion(max_time=10.0),SamplingHS3.MCMCConfig()) )


## To generate an actual sample, we have to resample the ns run:
rs_result = SamplingHS3.resample_ns_data( ns_result.ns_iter_data )

## The function sim_shell_compression simulates the shrinking process, and eval_integral
ssc_data = SamplingHS3.sim_shell_compression(ns_result.ns_iter_data)
volume_bounding_box = 20. * 20. * 20.
ei_result = SamplingHS3.eval_integral(volume_bounding_box,ns_result.ns_iter_data)

# plot the shell compression process:
p_comp = plot( log.( ssc_data[2]' ) )
#plot the integration process:
p_integ = plot( ei_result[3]'  , la=0.1 , leg=:none )
p_cum_integ = plot( cumsum( ei_result[3] , 2 )' , lc=:grey , la=0.1 , leg=:none )

## Plot the animcated sampling process (just because it's cool)
x_live_1 = SamplingHS3.reconstruct_live_samples(ns_result.ns_iter_data)
x_live  = x_live_1[1]
xf_live = x_live_1[2]
for zi in 1:length(xf_live)
    #xf_zi = 1:length(xf_live) )
    xf_zi = xf_live[zi]
    display( plot_3d_view( x_live[zi] , xf = xf_zi ) )
end



#binned_a = SamplingHS3.bin2d( ns_x_sequence[[1,2],:] , ones(size(ns_x_sequence,2)) , 40,40 , x -> sum(x) )
binned_a = SamplingHS3.bin2d( ns_x_sequence[[1,2],1:3000] , ones(3000) , 40,40 , x -> sum(x) )
heatmap(binned_a)
binned_b = SamplingHS3.bin2d( ns_x_sequence[[1,3],1:3000] , ones(3000) , 40,40 , x -> sum(x) )
heatmap(binned_b)
scatter!(rs_result[1][1,:],rs_result[1][2,:],rs_result[1][3,:])
