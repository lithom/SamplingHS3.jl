


"""
  run_mcmc( x_::Array{Float64,1} , f , threshold::Float64 , G::Array{Float64,2} , h::Array{Float64,1} , conf::DecorrelationConfig)

runs mcmc steps which sample from the probability distribution f_density, restricted to
the space {x|f_cost(x)<threshold}

# Note
if all steps fail, then the returned value of x (xf) will be NaN
"""
function run_mcmc( x_::Array{Float64,1} , x_f_density::Float64 , f_density , f_cost , threshold::Float64 , G::Array{Float64,2} , h::Array{Float64,1} , conf::DecorrelationConfig)
    d = length(x_)

    dd_xx       = NaN * ones(d,conf.steps+1)
    dd_xx_fdens = NaN * ones(conf.steps+1)
    dd_xx_fcost = NaN * ones(conf.steps+1)

    x  = copy(x_)

    use_density = !isa(f_density,Void)

    xfdens::Float64 = NaN
    if(use_density)
        if(isnan(x_f_density))
            xfdens = f_density(x)
        else
            xfdens = x_f_density
        end
    end

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

    xfcost::Float64 = NaN

    for zi in 1:conf.steps
        if(step_types_all[zi]==1)
            # HR STEP
            x_candidate_a = x - 4*rv_all[:,zi]  # times 4, because we want to be outside of viable space, if possible
            x_candidate_b = x + 4*rv_all[:,zi]
            #print("HRStep:\n")
            #print([x , x_candidate_a , x_candidate_b ])
            hr_result = hr_step( f_density , f_cost , threshold , G , h , x , x_candidate_a, x_candidate_b , hrconf ,  tdata_hr  )
            if( hr_result[1]>0 )
                cnt_successes += 1
                steps_data[zi] = +1
                x  = hr_result[2]
                # NOTE: following part will work only after upgrades to hr_step function..
                xfdens = hr_result[3]
                xfcost = hr_result[4]
                dd_xx[:,cnt_successes] = x
                dd_xx_fdens[cnt_successes] = xfdens
                dd_xx_fcost[cnt_successes] = xfcost
            else
                cnt_fails += 1
                steps_data[zi] = -1
            end
        end
        if(step_types_all[zi]==2)
            #MH STEP
            mh_result = mh_step( f_density , f_cost , threshold , G , h , x , xfdens , conf.cov_estimate , mh_max_tries , tdata_mh  )
            if(mh_result[1])
                cnt_successes += 1
                steps_data[zi] = +2
                x  = mh_result[2]
                xfdens = mh_result[3]
                xfcost = mh_result[4]
                dd_xx[:,cnt_successes] = x
                dd_xx_fdens[cnt_successes] = xfdens
                dd_xx_fcost[cnt_successes] = xfcost
            else
                cnt_fails += 1
                steps_data[zi] = -2
            end
        end
    end

    if(cnt_successes >= 1)
        # extract the successful steps data:
        dd_xx   = dd_xx[:,1:cnt_successes]
        dd_xx_fdens = dd_xx_fdens[1:cnt_successes]
        dd_xx_fcost = dd_xx_fcost[1:cnt_successes]
        # extract xv data:
        xv = dd_xx - [ x_ dd_xx[:,1:(size(dd_xx,2)-1)] ]
        # create dc:
        dd = DecorrelationData(xv , cnt_successes , cnt_fails , steps_data , dd_xx,dd_xx_fdens,dd_xx_fcost,tdata_hr,tdata_mh)
        return DecorrelationResult(x,xfdens,xfcost,dd)
    else
        print("mcmc failed..\n")
        print("step_types: \n")
        print(step_types)
        return DecorrelationResult(x,NaN,NaN,DecorrelationData( zeros(d,0) , cnt_successes , cnt_fails , steps_data ,zeros(d,0),[NaN],[NaN],tdata_hr,tdata_mh ))
    end
end


struct AdaptiveMCMCResult
    x::Vector{Float64}
    xfdens::Float64
    xfcost::Float64
    dd::DecorrelationData
    ecov_raw::Dict{Int,Array{Float64,2}}
    ecov_regularized::Dict{Int,Array{Float64,2}}
end


function run_adaptive_mcmc( steps::Int , x_::Array{Float64,1} , x_f_density::Float64 , f_density , f_cost , threshold::Float64 , G::Array{Float64,2} , h::Array{Float64,1} , conf_mcmc::MCMCConfig ; steps_between_cov_estim::Int = 5)

    output_level = 1

    d    = length(x_)
    xi   = copy(x_)
    xi_f = x_f_density


    # create decorrelation data vector
    dd::Vector{DecorrelationData} = Vector{DecorrelationData}()

    # create structures to save ecov info:
    data_ecov_raw = Dict{Int,Array{Float64,2}}()
    data_ecov_regularized = Dict{Int,Array{Float64,2}}()

    # Performed steps
    Xv = zeros(length(xi),0)

    # Init decorr_result struct
    sh_dcr::DecorrelationResult = DecorrelationResult()

    #
    total_steps = 0
    while(total_steps < steps)
        # 3. decorrelate
        # 3.1 estimate covariance..
        ecov_raw = eye(d)
        ecov     = NaN * ones(d,d)
        if( size(Xv,2)>4*d )
            if(size(Xv,2)>conf_mcmc.cov_estim.n_considered)
                ecov_raw = cov(Xv[ : , (end-conf_mcmc.cov_estim.n_considered+1):end ]' )
            else
                ecov_raw = cov(Xv')
            end
            sv_ecov = svd(ecov_raw)[2]
            ecov = ecov_raw + diagm(ones(d))*conf_mcmc.cov_estim.blowup_abs + diagm(ones(d)) * maximum(sv_ecov) * conf_mcmc.cov_estim.blowup_rel_maxsv + diagm(ones(d)) * minimum(sv_ecov) * conf_mcmc.cov_estim.blowup_rel_minsv

            if(output_level>=2)
                print("\necov_raw:\n$(ecov_raw)\n")
                print("ecov_regularized:\n$(ecov)\n")
                print("ecov_raw_sv: $(sv_ecov')\n")
                #print("ecov : blowup_abs: $(conf.mcmc.cov_estim.blowup_abs) , blowup_rel_max_sv: $(conf.mcmc.cov_estim.blowup_rel_maxsv*maximum(sv_ecov))) , blowup_rel_min_sv: $(conf.mcmc.cov_estim.blowup_rel_minsv*minimum(sv_ecov)) \n")
                print( @sprintf("ecov : blowup_abs: %6.3f  , blowup_rel_maxsv: %6.3f   , blowup_rel_minsv: %6.3f  \n", (conf_mcmc.cov_estim.blowup_abs) , (conf_mcmc.cov_estim.blowup_rel_maxsv*maximum(sv_ecov)) , (conf_mcmc.cov_estim.blowup_rel_minsv*minimum(sv_ecov)) ) )
            end
        else
            ecov = ecov_raw[:,:]
        end
        data_ecov_raw[total_steps]          = ecov_raw
        data_ecov_regularized[total_steps]  = ecov


        #threshold::Float64 = threshold
        steps_i            = min( steps-total_steps , steps_between_cov_estim )
        dc = DecorrelationConfig( steps_i , ecov , conf_mcmc.p_hr , conf_mcmc.p_mh )
        #sh_dcr = run_mcmc( xi , NaN, f_density , f_cost, threshold , G, h , dc )
        sh_dcr = run_mcmc( xi , xi_f , f_density , f_cost, threshold , G, h , dc )
        total_steps += sh_dcr.dd.steps_success

        # add decorrelation data
        push!(dd,sh_dcr.dd)
        # add performed steps
        if(sh_dcr.dd.steps_success>1)
            Xv  = [Xv    sh_dcr.dd.xv ]
        end

        xi   = sh_dcr.x
        xi_f = sh_dcr.xfdens
        #xi_fcost = sh_dcr.xfcost
    end

    # collapse the dds into one..
    dd_collapsed = reduce( Base.hcat , dd )
    
    return AdaptiveMCMCResult( sh_dcr.x , sh_dcr.xfdens , sh_dcr.xfcost , dd_collapsed , data_ecov_raw , data_ecov_regularized )
end


        #sh_dcr::Array{DecorrelationResult,1} = Array{DecorrelationResult}(size(X0,2))
        # parallel / not parallel
        #if(false)
        #    # 3.2 decorrelate samples..
        #    sh_dcr_af = Array{Future}(size(X0,2))
        #    for zi in 1:size(X0,2)
        #        sh_dcr_af[zi] = @spawn decorrelate( X0[:,zi] , f, threshold , G, h , dc )
        #    end
        #    for zi in 1:size(X0,2)
        #        sh_dcr[zi] = fetch(sh_dcr_af[zi])
        #    end
        #else
        #    for zi in 1:size(X0,2)
        #    sh_dcr[zi] = run_mcmc( X0[:,zi] , NaN, Void(), f, threshold , G, h , dc )#decorrelate( X0[:,zi] , f, threshold , G, h , dc )
        #        #sh_dcr[zi] = decorrelate( X0[:,zi] , f, threshold , G, h , dc )
        #        #if(sh_dcr[zi].)
        #    end
        #end

        #X_new   = Array{Float64}(d,0)
        #X_new_f = Array{Float64}(0)
        #Xv_new  = Array{Float64}(d,0)
        ##for zi in 1:length(sh_dcr)
        #    push!(dd,sh_dcr[zi].dd)
        #    if(sh_dcr[zi].dd.steps_success>1)
        #        X_new   = [ X_new sh_dcr[zi].x ]
        #        Xv_new  = [ Xv_new  sh_dcr[zi].dd.xv ]
        #        X_new_f = [ X_new_f ; sh_dcr[zi].xfcost ]
        #        total_steps += sh_dcr[zi].dd.steps_success
        #    end
        ##end
        #X   = [X     X_new   ]
        #Xv  = [Xv    Xv_new  ]
        #X_f = [X_f ; X_new_f ]
    #end
#end
