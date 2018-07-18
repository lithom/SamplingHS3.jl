


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
