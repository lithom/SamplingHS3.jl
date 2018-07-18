
using SamplingHS3

## Create gaussian mixture
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
# bounding box
bb_G = [eye(2);-eye(2)]
bb_h = [10.;10;10;10]

# quick plot as sanity check
using Plots
xx  = SamplingHS3.draw_gaussmixture( gm2 , 100000 )
SamplingHS3.gr_hist_2d( xx , nbins=40 )

## Use SamplingHS3

# Input:
# Density function: (x) -> SamplingHS3.pdf_gaussmixture_nonnormalized(gm2,x)
# Bounding box, given by bb_G,bb_h


# Sampling Steps:
# 1. generate uniform samples over polytope
X_init = SamplingHS3.sample_unif_hr( bb_G , bb_h , zeros(2) , 1000 , 2000 )

# 2. create "cost" function, i.e. for a pdf given
#    by f(x), we use f_cost(x) = -log(f(x))
f_cost = (x) -> -log(SamplingHS3.pdf_gaussmixture_nonnormalized(gm2,x))

# Optional: evaluate cost of initial samples (can be supplied to the ns function)
# x_init_cost = [  f_cost(X_init[:,zi]) for zi in (1:2000) ]


# 3. Run NS:
# Required parameters are:
# (obviously): initial samples, cost function, polytope description
# (and):       the NSConfig struct
ns_result = SamplingHS3.nested_sampling( X_init[:,:] , f_cost , bb_G , bb_h , SamplingHS3.NSConfig(8,0.95,SamplingHS3.StopCriterion(max_time=10.0),SamplingHS3.MCMCConfig()) )

# 4. Plot evolution of live samples during ns run (Just because it looks cool..)
x_live_samples_1 = SamplingHS3.reconstruct_live_samples(ns_result.ns_iter_data)
x_live  = x_live_samples_1[1]
xf_live = x_live_samples_1[2]
for zi in 1:length(x_live)
    #zcolor_i = 1:size(x_live[zi],2)
    zcolor_i = xf_live[zi]
    display( scatter(x_live[zi][1,:] , x_live[zi][2,:]  , ma = 0.75 , zcolor=zcolor_i , m=2, msw=0 , c=:lightrainbow) )
end

# 5. Generate an actual sample of the probability distribution from the ns run:
