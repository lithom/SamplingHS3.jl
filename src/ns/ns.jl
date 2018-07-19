
"""
DecorrelationConfig contains config for the decorrelation steps
...
# Fields:
  - `steps::Int64`: number of decorrelation steps
  - `cov_estimate::Array{Float64,2}`: estimated (local) covariance,
  - `p_hr::Float64`: probability HR Steps
  - `p_mh::Float64`: probability MH Steps
...
"""
struct DecorrelationConfig
    steps::Int64
    cov_estimate::Array{Float64,2} # estimated (local) covariance,
    p_hr::Float64 # probability HR Steps
    p_mh::Float64 # probability MH Steps
end

"""
DecorrelationData
  collects data on performed decorrelations..
"""
struct DecorrelationData
    xv::Array{Float64,2} # all steps performed, chronologically (roughly)
    steps_success::Int64
    steps_fail::Int64
    steps::Vector{Int64} # +1 HR success, +2 MH success, -1 HR fail, -2 MH fail

    xx::Array{Float64,2} # on the way collected samples..
    xx_fdens::Array{Float64,1} # on the way collected samples..
    xx_fcost::Array{Float64,1} # on the way collected samples..


    sdata_hr::HRSamplingTrackingData
    sdata_mh::MHSamplingTrackingData
end

function DecorrelationData()
    return DecorrelationData( Array{Float64}(0,0) , Int64(0) , Int64(0) , Vector{Int64}() , Array{Float64}(0,0) , Array{Float64,1}() , Array{Float64,1}(), HRSamplingTrackingData() , MHSamplingTrackingData() )
end
function Base.hcat(a::DecorrelationData,b::DecorrelationData)
    return DecorrelationData( [a.xv b.xv] , a.steps_success+b.steps_success , a.steps_fail+b.steps_fail , [a.steps;b.steps] , [a.xx b.xx] , [a.xx_fdens ; b.xx_fdens] , [a.xx_fcost ; b.xx_fcost] , [ a.sdata_hr b.sdata_hr] , [ a.sdata_mh b.sdata_mh] )
end

#function DecorrelationData( dc1 , dc2 )
#    return DecorrelationData( [ dc1.xv dc2.xv ] , [dc1.steps_success ,dc2.steps_success ] , [dc1.steps_fail ,dc2.steps_fail ] , [dc1.xx ,dc2.xx ]  , [dc1.xxf ,dc2.xxf ]  )
#end

struct DecorrelationResult
    x::Vector{Float64}
    xfdens::Float64
    xfcost::Float64
    dd::DecorrelationData
end
function DecorrelationResult()
    return DecorrelationResult([NaN],NaN,NaN,DecorrelationData())
end


"""
CovEstimConfig configures the (local) covariance estimation
...
# Fields:
  - `n_considered::Int64`: the number of recently performed steps that we consider for estimating the covariance
  - `blowup_abs::Float64`: "Haario" regularization: add ( covestim_blowup_abs * diag(ones(d)) ) to the estimated cov. matrix
  - `blowup_rel_maxsv::Float64`: "Haario" rel regularization: add  ( covestim_blowup_abs * (diag(ones(d))) * norm(cov) ) , i.e. add relative to the frobenius norm ( max singular value)
  - `blowup_rel_minsv::Float64`: "Haario" rel regularization: add  ( covestim_blowup_abs * (diag(ones(d))) * min(svd(cov)[2]) ) , i.e. add relative to the smallest singular value
...


# Explanation about covestim_blowup parameters:
  blowup_abs is a standard procedure to regularize covariance estimation in adaptive mcmc. However, during
  ns, the estimated covariance will decrease, therefore, evetually, absolute blowup will prevent adaptiveness
  of the mcmc procedure.
  Therefore, we add the possibility to add blowup relative to singular values of the covariance matrix.
  The rel_maxsv blowup factor also limits the maximum non-isotropy of the estimated covariance in an intuitive way.
"""
struct CovEstimConfig
    n_considered::Int64 # the number of recently performed steps that we consider for estimating the covariance
    blowup_abs::Float64 # "Haario" regularization: add   ( covestim_blowup_abs * diag(ones(d)) ) to the estimated cov. matrix
    blowup_rel_maxsv::Float64 # "Haario" rel regularization: add  ( covestim_blowup_abs * (diag(ones(d))) * norm(cov) )        , i.e. add relative to the frobenius norm ( max singular value)
    blowup_rel_minsv::Float64 # "Haario" rel regularization: add  ( covestim_blowup_abs * (diag(ones(d))) * min(svd(cov)[2]) ) , i.e. add relative to the smallest singular value
end

# create a rather conservative estimation procedure by default
function CovEstimConfig(; n_considered=-1 , blowup_abs=NaN , blowup_rel_maxsv=NaN , blowup_rel_minsv=NaN)
    if(n_considered<0); n_considered = 200; end
    if(isnan(blowup_abs)); blowup_abs = 0.0001; end
    if(isnan(blowup_rel_maxsv)); blowup_rel_maxsv = 0.005; end # this also limits the ratio of small/large singular values to rouglhy 200
    if(isnan(blowup_rel_minsv)); blowup_rel_minsv = 0.1;   end # this also limits the ratio of small/large singular values
    return CovEstimConfig(n_considered,blowup_abs,blowup_rel_maxsv,blowup_rel_minsv)
end

"""
MCMCConfig configures the behavior of the mcmc decorrelation procedure
...
# Fields:
  - `p_hr::Float64`: probability of hr steps
  - `p_mh::Float64` : probability of mh steps
  - `cov_estim::CovEstimConfig` : covariance estimation configuration
  - `cov_initial::Array{Float64,2}` : initial covariance. Here, [NaN] is an alias for eye(dim)
...

"""
struct MCMCConfig
    p_hr::Float64
    p_mh::Float64
    cov_estim::CovEstimConfig
    cov_initial::Array{Float64,2} # if [NaN], then we take eye(dim)
end

function MCMCConfig()
    return MCMCConfig(0.5,0.5,CovEstimConfig(),reshape([NaN],1,1))
end

"""
NSConfig configures a run of the ns procedure
...
# Fields:
  - `steps_decorrelation::Int64`: the number of decorrelation steps
  - `shrinking_ratio::Float64` : the ratio by which we shrink in every iteration
  - `stop::StopCriterion` : the stop criterion
  - `mcmc::MCMCConfig` : configuration of the mcmc (decorrelation) procedure
...

"""
struct NSConfig
    steps_decorrelation::Int64
    shrinking_ratio::Float64
    stop::StopCriterion
    mcmc::MCMCConfig
    #decorrelation::DecorrelationConfig
    #save_live_samples::Bool
end

function NSConfig(; steps_decorrelation=-1 , shrinking_ratio=NaN , stop=StopCriterion() , mcmc=MCMCConfig() )#, save_live_samples=false)
    if(steps_decorrelation<0); steps_decorrelation=8; end
    if(isnan(shrinking_ratio)); shrinking_ratio=0.9; end
    return NSConfig(10,0.9,stop,mcmc)
end

struct NSIterData
    n_live_samples_a::Int64 # num samples before iteration
    n_live_samples_b::Int64 # num samples after removing sorted out samples
    n_live_samples_c::Int64 # num samples after replaing sorted out samples
    threshold::Float64      # cutoff threshold in this iteration (i.e. max fcost of samples "b" )
    X_sorted_out::Array{Float64} # sorted out samples
    Xf_sorted_out::Vector{Float64} # fvals of sorted out samples
    dd::Vector{DecorrelationData} # contains the data about the decorrelation steps
    ecov_raw::Array{Float64,2}
    ecov_regularized::Array{Float64,2}
end

struct NestedSamplingResult
    #X_sequence::Array{Float64,2}
    #X_sequence_f::Array{Float64,1}
    ns_iter_data::Array{NSIterData}
    #decorrelation_data::Array{DecorrelationData,1}
end

#function nested_sampling( X_init::Array{Float64,2} , f::Any , G::Array{Float64,2} , h::Array{Float64,1} , conf::NSConfig  ; X_f_init::Vector{Float64}=[] )

"""
  nested_sampling( X_init::Array{Float64,2} , f::Any , G::Array{Float64,2} , h::Array{Float64,1} , conf::NSConfig ; X_f_init::Vector{Float64}=[NaN] )

returns NestedSamplingResult

this function runs the nested sampling procedure
"""
function nested_sampling( X_init::Array{Float64,2} , f::Any , G::Array{Float64,2} , h::Array{Float64,1} , conf::NSConfig ; X_f_init::Vector{Float64}=[NaN] )

    output_level = 2

    X = copy(X_init)
    n_samples = size(X,2)
    d         = size(X,1)

    X_f::Array{Float64,1} = Array{Float64,1}()
    # init X_f if not given
    if(isnan(X_f_init[1]))
        X_f = Vector{Float64}(size(X_init,2))
        for zi in 1:size(X,2); X_f[zi] = f(X[:,zi]); end
    else
        X_f = copy(X_f_init)
    end

    # this data structure keeps track of the performed steps
    # create one for each step..
    #dd::Vector{DecorrelationData} = Vector{DecorrelationData}()

    itd::Vector{NSIterData}       = Vector{NSIterData}()

    # here we save the performed steps..
    Xv::Array{Float64,2}          = Array{Float64}(d,0)

    # this keeps track of all generated samples
    #X_sequence::Array{Float64,2}   = Array{Float64}(d,0)
    #X_sequence_f::Array{Float64,1} = Array{Float64}(0)

    # start the ns
    stop_ns = false
    t_start = time()
    iteration = 0
    while( !stop_ns )
        iteration += 1

        # profiling data
        iterdata_n_samples_a = size(X,2)

        # save live samples
        #X_saved_live_samples = Array{Float64,2}()
        #if(conf.save_live_samples)
        #    X_saved_live_samples   = X[:,:]
        #    X_saved_live_samples_f = X_f[:]
        #end


        # 1. sort out samples below threshold given by shrinking_ratio:
        println(X_f)
        f_cut = quantile( X_f , conf.shrinking_ratio)

        b_ok  = X_f .<= f_cut
        if(all(b_ok)); b_ok[ findmin(X_f)[2] ] = 0; end # at least sort out one.. :)

        #X_sequence   = [ X_sequence    X[ : , .~b_ok ] ]
        #X_sequence_f = [ X_sequence_f ; X_f[     .~b_ok ] ]
        iterdata_X_sorted_out   = X[ : , .~b_ok ]
        iterdata_Xf_sorted_out  = X_f[ .~b_ok ]

        #iterdata_fval_sorted_out = X_f[ .~b_ok ]

        X     = X[ : , b_ok ]
        X_f   = X_f[ b_ok ]

        iterdata_n_samples_b = size(X,2)


        # 2. select starting point samples
        n_replace = n_samples - size(X,2)
        idx_startingpoints = rand( 1:size(X,2) , n_replace )
        X0 = X[ : , idx_startingpoints ]

        # print iteration data..
        print( @sprintf("It %04d | %6.3f | %04d \n",iteration,f_cut,size(X,2)) )


        # 3. decorrelate
        # 3.1 estimate covariance..
        ecov_raw = eye(d)
        ecov     = NaN * ones(d,d)
        if( size(Xv,2)>4*d )
            if(size(Xv,2)>conf.mcmc.cov_estim.n_considered)
                ecov_raw = cov(Xv[ : , (end-conf.mcmc.cov_estim.n_considered+1):end ]' )
            else
                ecov_raw = cov(Xv')
            end
            sv_ecov = svd(ecov_raw)[2]
            ecov = ecov_raw + diagm(ones(d))*conf.mcmc.cov_estim.blowup_abs + diagm(ones(d)) * maximum(sv_ecov) * conf.mcmc.cov_estim.blowup_rel_maxsv + diagm(ones(d)) * minimum(sv_ecov) * conf.mcmc.cov_estim.blowup_rel_minsv

            if(output_level>=2)
                print("\necov_raw:\n$(ecov_raw)\n")
                print("ecov_regularized:\n$(ecov)\n")
                print("ecov_raw_sv: $(sv_ecov')\n")
                #print("ecov : blowup_abs: $(conf.mcmc.cov_estim.blowup_abs) , blowup_rel_max_sv: $(conf.mcmc.cov_estim.blowup_rel_maxsv*maximum(sv_ecov))) , blowup_rel_min_sv: $(conf.mcmc.cov_estim.blowup_rel_minsv*minimum(sv_ecov)) \n")
                print( @sprintf("ecov : blowup_abs: %6.3f  , blowup_rel_maxsv: %6.3f   , blowup_rel_minsv: %6.3f  \n", (conf.mcmc.cov_estim.blowup_abs) , (conf.mcmc.cov_estim.blowup_rel_maxsv*maximum(sv_ecov)) , (conf.mcmc.cov_estim.blowup_rel_minsv*minimum(sv_ecov)) ) )
            end
        else
            if( isnan(conf.mcmc.cov_initial[1]) )
                ecov = ecov_raw[:,:]
            else
                ecov = conf.mcmc.cov_initial[:,:]
            end
        end

        # create decorrelation data vector
        dd::Vector{DecorrelationData} = Vector{DecorrelationData}()

        threshold::Float64 = f_cut
        dc = DecorrelationConfig( conf.steps_decorrelation , ecov , conf.mcmc.p_hr , conf.mcmc.p_mh )
        sh_dcr::Array{DecorrelationResult,1} = Array{DecorrelationResult}(size(X0,2))
        # parallel / not parallel
        if(false)
            # 3.2 decorrelate samples..
            sh_dcr_af = Array{Future}(size(X0,2))
            for zi in 1:size(X0,2)
                sh_dcr_af[zi] = @spawn decorrelate( X0[:,zi] , f, threshold , G, h , dc )
            end
            for zi in 1:size(X0,2)
                sh_dcr[zi] = fetch(sh_dcr_af[zi])
            end
        else
            for zi in 1:size(X0,2)
                sh_dcr[zi] = run_mcmc( X0[:,zi] , NaN, Void(), f, threshold , G, h , dc )#decorrelate( X0[:,zi] , f, threshold , G, h , dc )
                #sh_dcr[zi] = decorrelate( X0[:,zi] , f, threshold , G, h , dc )
                #if(sh_dcr[zi].)
            end
        end

        X_new   = Array{Float64}(d,0)
        X_new_f = Array{Float64}(0)
        Xv_new  = Array{Float64}(d,0)
        for zi in 1:length(sh_dcr)
            push!(dd,sh_dcr[zi].dd)
            if(sh_dcr[zi].dd.steps_success>1)
                X_new   = [ X_new sh_dcr[zi].x ]
                Xv_new  = [ Xv_new  sh_dcr[zi].dd.xv ]
                X_new_f = [ X_new_f ; sh_dcr[zi].xfcost ]
            end
        end
        X   = [X     X_new   ]
        Xv  = [Xv    Xv_new  ]
        X_f = [X_f ; X_new_f ]

        iterdata_n_samples_c = size(X,2)

        push!(itd,NSIterData(iterdata_n_samples_a,iterdata_n_samples_b,iterdata_n_samples_c,f_cut,iterdata_X_sorted_out,iterdata_Xf_sorted_out,dd,ecov_raw,ecov))

        # check if we stop:
        if( time()-t_start>conf.stop.max_time)
            stop_ns = true

            # put remaining samples into the last iteration data
            push!(itd,NSIterData(iterdata_n_samples_c,0,0,-Inf,X,X_f,dd,ecov_raw,ecov))
        else
        end
    end

    #return NestedSamplingResult( [ X_sequence X ] , [ X_sequence_f;X_f ] , itd , dd )
    return NestedSamplingResult( itd  )
end



"""
  decorrelate( x_::Array{Float64,1} , f , threshold::Float64 , G::Array{Float64,2} , h::Array{Float64,1} , conf::DecorrelationConfig)

decorrelates the given samples inside {x|f(x)<threshold} by running mcmc steps

# Note
if all steps fail, then the returned value of x (xf) will be NaN

function decorrelate( x_::Array{Float64,1} , f , threshold::Float64 , G::Array{Float64,2} , h::Array{Float64,1} , conf::DecorrelationConfig)
    d = length(x_)

    dd_xx   = NaN * ones(d,conf.steps+1)
    dd_xx_f = NaN * ones(conf.steps+1)

    x  = copy(x_)
    xf = NaN

    hrconf = HRStepConfig()
    tdata_hr  = HRSamplingTrackingData()
    tdata_mh  = MHSamplingTrackingData()
    steps_data = Vector{Int64}(conf.steps)

    # config:
    mh_max_tries = 8

    cnt_successes = 0
    cnt_fails = 0
    rv_all = mvnrnd(x,conf.cov_estimate,conf.steps)
    step_types_all = sample_proportional([conf.p_hr;conf.p_mh],conf.steps)

    for zi in 1:conf.steps
        if(step_types_all[zi]==1)
            # HR STEP
            x_candidate_a = x - 4*rv_all[:,zi]  # times 4, because we want to be outside of viable space, if possible
            x_candidate_b = x + 4*rv_all[:,zi]
            #print("HRStep:\n")
            #print([x , x_candidate_a , x_candidate_b ])
            hr_result = hr_step( Void() , f , threshold , G , h , x , x_candidate_a, x_candidate_b , hrconf ,  tdata_hr  )
            if( hr_result[1]>0 )
                cnt_successes += 1
                steps_data[zi] = +1
                x  = hr_result[2]
                xf = hr_result[3]
                dd_xx[:,cnt_successes] = x
                dd_xx_f[cnt_successes] = xf
            else
                cnt_fails += 1
                steps_data[zi] = -1
            end
        end
        if(step_types_all[zi]==2)
            #MH STEP
            mh_result = mh_step( Void() , f , threshold , G , h , x , NaN , conf.cov_estimate , mh_max_tries , tdata_mh  )
            if(mh_result[1])
                cnt_successes += 1
                steps_data[zi] = +2
                x  = mh_result[2]
                xf = mh_result[4]
                dd_xx[:,cnt_successes] = x
                dd_xx_f[cnt_successes] = xf
            else
                cnt_fails += 1
                steps_data[zi] = -2
            end
        end
    end

    if(cnt_successes >= 1)
        # extract the successful steps data:
        dd_xx   = dd_xx[:,1:cnt_successes]
        dd_xx_f = dd_xx_f[1:cnt_successes]
        # extract xv data:
        xv = dd_xx - [ x_ dd_xx[:,1:(size(dd_xx,2)-1)] ]
        # create dc:
        dd = DecorrelationData(xv , cnt_successes , cnt_fails , steps_data , dd_xx,[NaN],dd_xx_f,tdata_hr,tdata_mh)
        return DecorrelationResult(x,NaN,xf,dd)
    else
        return DecorrelationResult(x,NaN,NaN,DecorrelationData( zeros(d,0) , cnt_successes , cnt_fails ,zeros(d,0),[NaN],tdata_hr,tdata_mh ))
    end
end
"""
