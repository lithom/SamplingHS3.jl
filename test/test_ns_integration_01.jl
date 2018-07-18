using SamplingHS3

struct gm_test
    gm::SamplingHS3.GaussMixture
    gm_p::Vector{SamplingHS3.GaussParam}

    G::Array{Float64,2}
    h::Array{Float64,1}
    vol::Float64
end

function test_gauss_mixture( means , sigmas , bounds )
    gmx = Vector{SamplingHS3.GaussParam}(length(means))
    for zi=1:length(means)
        gmx[zi] = SamplingHS3.GaussParam(means[zi],sigmas[zi])
    end
    bb1 = SamplingHS3.get_bounding_box(bounds[:,1],bounds[:,2])
    G   = bb1[1]
    h   = bb1[2]
    vol = prod( bounds[:,2]-bounds[:,1] )
    return gm_test( SamplingHS3.GaussMixture(gmx,ones(length(gmx))),gmx,G,h , vol)
end

function run_gm_test( gmt , n_live , n_steps_decorr , sfactor , time_s )
    d = size(gmt.G,2)
    f_cost_gm = (x) -> -log(SamplingHS3.pdf_gaussmixture_nonnormalized(gmt.gm,x)+1.e-100)
    X_init = SamplingHS3.sample_unif_hr( gmt.G , gmt.h , zeros(d) , n_live , 2000 )
    ns_result = SamplingHS3.nested_sampling( X_init , f_cost_gm , gmt.G , gmt.h , SamplingHS3.NSConfig(n_steps_decorr,sfactor,SamplingHS3.StopCriterion(max_time=time_s),SamplingHS3.MCMCConfig()) )

    ns_x_sequence = reduce( (x,y) -> [x y] , map( x -> x.X_sorted_out , ns_result.ns_iter_data ) )

    # test reconstruction of live samples:
    ns_x_livesamples = SamplingHS3.reconstruct_live_samples(ns_result.ns_iter_data)
    # test resampling:
    rs_result = SamplingHS3.resample_ns_data( ns_result.ns_iter_data )

    #scatter( ns_x_sequence[1,:] , ns_x_sequence[2,:] , ns_x_sequence[3,:] , ma = 0.75 , zcolor=(1:size(ns_x_sequence,2)) , m=2, msw=0 , c=:lightrainbow , lab="grad" )
    ei_result = SamplingHS3.eval_integral(gmt.vol,ns_result.ns_iter_data)
    return ei_result[1]
end

gm_test_a = test_gauss_mixture( [[4;4;4],[-4;-4;-4]] , [eye(3),eye(3)] , repmat([-10 10],3,1) )
gm_test_a_results = Vector{Float64}(10)
for zi=1:length(gm_test_a_results)
    ival = run_gm_test( gm_test_a , 200 , 6 , 0.8 , 2. )
    gm_test_a_results[zi] = ival
end
