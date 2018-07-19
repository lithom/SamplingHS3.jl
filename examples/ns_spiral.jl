

function plot_3d_view( x_3d::Array{Float64,2} ; xf = [] )
    if(isempty(xf)); xf = ones(size(x_3d)); end
    gr()
    p1 = scatter( x_3d[1,:] , x_3d[2 ,:] , x_3d[3,:] , zcolor = xf , ma = 0.2 , m=3, msw=0 , c=:lightrainbow , leg=:none)
    p2 = scatter( x_3d[1,:] , x_3d[2 ,:]  , zcolor = xf , ma = 0.2 , m=3, msw=0 , c=:lightrainbow , leg=:none)
    p3 = scatter( x_3d[1,:] , x_3d[3 ,:]  , zcolor = xf , ma = 0.2 , m=3, msw=0 , c=:lightrainbow , leg=:none )
    p4 = scatter( x_3d[2,:] , x_3d[3 ,:]  , zcolor = xf , ma = 0.2 , m=3, msw=0 , c=:lightrainbow , leg=:none )

    display(plot(p1,p2,p3,p4,layout=(2,2)))
end


# dist-to-spiral function
spiral_3d = reduce( hcat , [ [sin(0.5*z);cos(0.5*z);z] for z in linspace(-20.,20.,60)] )


li_spir   = SamplingHS3.create_line(spiral_3d)

#f_cost_spir_3d = (x::Vector{Float64}) -> maximum( [0.5 ; min( SamplingHS3.dist2line(li_spir,x[:,:] + [2.5;2.5;0] )[1] , SamplingHS3.dist2line(li_spir,x[:,:] - [2.5;2.5;0] )[1] ) ] )
f_cost_spir_3d = (x::Vector{Float64}) -> begin sdi = SamplingHS3.dist2line(li_spir,x[:,:])[1]; (sdi<0.5)? 0.5-(0.5-sdi)^4 : sdi end

# polytope
xlim = [-20.,20.]; ylim = [-20.,20.]; zlim = [-20.,20.];
bb_G =[eye(3);-eye(3)]
bb_h =[xlim[2];ylim[2];zlim[2];-xlim[1];-ylim[1];-zlim[1]]



using Plots
Plots.gr()
Plots.plot(spiral_3d[1,:],spiral_3d[2,:],spiral_3d[3,:])

X_init = SamplingHS3.unifsamp_box( [-20 20;-20 20;-20 20] , 400 )
ns_result_01 = SamplingHS3.nested_sampling( X_init[:,:] , f_cost_spir_3d , bb_G , bb_h , SamplingHS3.NSConfig(6,0.75,SamplingHS3.StopCriterion(max_time=30.0),SamplingHS3.MCMCConfig()) )

## Plot the animated sampling process (just because it's cool)
x_live_1 = SamplingHS3.reconstruct_live_samples(ns_result_01.ns_iter_data)
x_live  = x_live_1[1]
xf_live = x_live_1[2]
for zi in 1:length(xf_live)
    #xf_zi = 1:length(xf_live) )
    xf_zi = xf_live[zi]
    display( plot_3d_view( x_live[zi] , xf = xf_zi ) )
end
