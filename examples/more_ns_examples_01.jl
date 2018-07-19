
using SamplingHS3



gm2_mix = []
for zi=1:10;
    found_es = false
    vd_pre   = NaN * ones(2,2)
    while(!found_es)
        vd = qr(rand(2,2))[1]
        vd_pre = vd*diagm(rand(2))*vd'
        found_es = issymmetric(vd_pre)
    end
    gpi = SamplingHS3.GaussParam( 8*rand(2) , vd_pre )
    gm2_mix = [ gm2_mix ;  gpi ]
end

gm2 = SamplingHS3.GaussMixture(gm2_mix,0.5+rand(10))
xx  = SamplingHS3.draw_gaussmixture( gm2 , 100000 )


using Plots
#SamplingHS3.gr_hist_2d( xx , nbins=40 )


bb_G = [eye(2);-eye(2)]
bb_h = [10.;10;10;10]

X_init = SamplingHS3.sample_unif_hr( bb_G , bb_h , zeros(2) , 1000 , 2000 )

# now sample dat with nested sampling
f_cost = (x) -> -log(SamplingHS3.pdf_gaussmixture_nonnormalized(gm2,x))

x_unif = SamplingHS3.unifsamp_box([-10. 10;-10 10],2000)
x_unif_cost = [  f_cost(x_unif[:,zi]) for zi in (1:2000) ]

x_unif_cost = [  SamplingHS3.mvnpdf([0.;0],1.0*eye(2),x_unif[:,zi:zi]) for zi in (1:2000) ]

using Plots
gr()
gui()
Plots.scatter( x_unif[1,:] , x_unif[2,:] ; zcolor=x_unif_cost[:],c=:lightrainbow,ms=4.0 )

# sample..
#ns_result = SamplingHS3.nested_sampling( X_init , f_cost , bb_G , bb_h , SamplingHS3.NSConfig(8,0.8,SamplingHS3.StopCriterion(max_time=3.),SamplingHS3.MCMCConfig()) )
ns_result = SamplingHS3.nested_sampling( X_init[:,:] , f_cost , bb_G , bb_h , SamplingHS3.NSConfig(8,0.95,SamplingHS3.StopCriterion(max_time=10.0),SamplingHS3.MCMCConfig()) )

x_live_samples_1 = SamplingHS3.reconstruct_live_samples(ns_result.ns_iter_data)
x_live  = x_live_samples_1[1]
xf_live = x_live_samples_1[2]
# plot live samples..
# P_subplots = plot(layout=(4,4) , size=(2000,2000))
for zi in 1:length(x_live)
    display( Plots.scatter(x_live[zi][1,:] , x_live[zi][2,:]  , ma = 0.75 , zcolor=(1:size(x_live[zi],2)) , m=2, msw=0 , c=:lightrainbow) )
    #scatter!( x_live[zi][1,:] , x_live[zi][2,:] , subplot = zi , ma = 0.75 , zcolor=(1:size(x_live[zi],2)) , m=2, msw=0 , c=:lightrainbow , lab="grad" )
end
gui()
gui(plot(P_subplots , size=(4000,4000)))

using Plots
gr()
for n in 1:100
  display(histogram(randn(10000)))
  sleep(0.05)
end


ns_x_sequence = reduce( (x,y) -> [x y] , map( x -> x.X_sorted_out , ns_result.ns_iter_data ) )
scatter( ns_x_sequence[1,:] , ns_x_sequence[2,:] , ma = 0.75 , zcolor=(1:size(ns_x_sequence,2)) , m=2, msw=0 , c=:lightrainbow , lab="grad" )

ssc_data = SamplingHS3.sim_shell_compression(ns_result.ns_iter_data)
ei_result = SamplingHS3.eval_integral(400.,ns_result.ns_iter_data)

rs_result = SamplingHS3.resample_ns_data( ns_result.ns_iter_data )

scatter!(rs_result[1][1,:],rs_result[1][2,:])

scatter(ns_result.ns_iter_data[16].X_sorted_out[1,:],ns_result.ns_iter_data[16].X_sorted_out[2,:])

scatter!( ns_result.X_sequence[1,1:500] , ns_result.X_sequence[2,1:500] )
scatter( ns_result.X_sequence[1,1:500] , ns_result.X_sequence[2,1:500] )
scatter!( ns_result.X_sequence[1,1000:1500] , ns_result.X_sequence[2,1000:1500] )
scatter( ns_result.X_sequence[1,1300:1500] , ns_result.X_sequence[2,1300:1500] )
plot!()


# create plot of some deocrrelation trajectories:
heatmap()


## 3 d test:

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
xx3  = SamplingHS3.draw_gaussmixture( gm3 , 100000 )
bb3_G = [eye(3);-eye(3)]
bb3_h = [10.;10;10;10;10;10]


f_cost_3d = (x) -> -log(SamplingHS3.pdf_gaussmixture_nonnormalized(gm3,x))
scatter(xx3[1,1:2000] , xx3[2 ,1:2000] , xx3[3,1:2000] , ma = 0.1 , m=2, msw=0 , c=:lightrainbow , lab="grad")

# sample..
X_init = SamplingHS3.sample_unif_hr( bb3_G , bb3_h , zeros(3) , 400 , 2000 )
ns_result = SamplingHS3.nested_sampling( X_init , f_cost_3d , bb3_G , bb3_h , SamplingHS3.NSConfig(8,0.8,SamplingHS3.StopCriterion(max_time=20.)) )


ns_x_sequence = reduce( (x,y) -> [x y] , map( x -> x.X_sorted_out , ns_result.ns_iter_data ) )
scatter( ns_x_sequence[1,:] , ns_x_sequence[2,:] , ns_x_sequence[3,:] , ma = 0.75 , zcolor=(1:size(ns_x_sequence,2)) , m=2, msw=0 , c=:lightrainbow , lab="grad" )
scatter( ns_x_sequence[1,end-3000:end] , ns_x_sequence[2,end-3000:end] , ns_x_sequence[3,end-3000:end] , ma = 0.75 , zcolor=(1:2001) , m=2, msw=0 , c=:lightrainbow , lab="grad" )

ssc_data = SamplingHS3.sim_shell_compression(ns_result.ns_iter_data)
ei_result = SamplingHS3.eval_integral(400.,ns_result.ns_iter_data)
rs_result = SamplingHS3.resample_ns_data( ns_result.ns_iter_data )

scatter!(rs_result[1][1,:],rs_result[1][2,:],rs_result[1][3,:])


#binned_a = SamplingHS3.bin2d( ns_x_sequence[[1,2],:] , ones(size(ns_x_sequence,2)) , 40,40 , x -> sum(x) )
binned_a = SamplingHS3.bin2d( ns_x_sequence[[1,2],1:3000] , ones(3000) , 40,40 , x -> sum(x) )
heatmap(binned_a)
binned_b = SamplingHS3.bin2d( ns_x_sequence[[1,3],1:3000] , ones(3000) , 40,40 , x -> sum(x) )
heatmap(binned_b)
