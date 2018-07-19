"""
  xrayshoot(A::Array{Float64,2} , b::Array{Float64,1} , x0::Array{Float64,2} , v::Array{Float64,2})

returns: ( lambda_min , lambda_max , x_intersection_min , x_intersection_max )

# xrayshoot considers the system of linear inequalities A*x<b, and assumes
that all columns of x0 lie inside the solution set.
Then this function returns the most negative / most postitive lambdas
such that x0_i + lambda_i * v_i <= A*b holds.


# TEST:
 A1=[eye(2);-eye(2)]; b1=[1;1;1;1];
   lambda_min , lambda_max , x_intersection_min , x_intersection_max = XRayShoot(A1,b1,zeros(2,4),[[1;1],[1;0],[-1;1],[0.1;0.2]]);
   lambda_min , lambda_max , x_intersection_min , x_intersection_max = XRayShoot(A1,b1,zeros(2,1),[[1;1],[1;0],[-1;1],[0.1;0.2]]);
"""
function xrayshoot( A::Array{Float64,2} , b::Array{Float64,1} , x0::Array{Float64,2} , v::Array{Float64,2} )
    # if only one point, then replicate it once for every supplied
    # direction
    if ( (size(x0,2)==1) && (size(v,2)>1) )
        x0 = repmat(x0,1,size(v,2))
    end

    ## 1. compute intersections with constraints:
    AZ = A*v
    S::Array{Float64,2} = broadcast(-,b,A*x0)./AZ
    smax = copy(S)
    smax[AZ.<=0] = Inf
    smax = minimum(smax,1)
    lambda_max = smax;

    smin = S
    smin[AZ.>=0] = -Inf
    #print(smin)
    smin = maximum(smin,1)
    lambda_min = smin

    ## also return the intersecting points:
    #if(nargout >= 3),
    x_intersection_min = x0 + broadcast( * , lambda_min , v )
    x_intersection_max = x0 + broadcast( * , lambda_max , v )

    return lambda_min , lambda_max , x_intersection_min , x_intersection_max
end

mutable struct HRSamplingTrackingData
    X::Vector{Vector{Float64}}
    X_fdens::Vector{Float64}
    X_fcost::Vector{Float64}
    X_viable::Vector{Bool}
    X_opcode::Vector{Int64}

    fevals_expand::Int64
    fevals_bisect::Int64
    fevals_linesearch::Int64
    fevals_fdens::Int64
end
function HRSamplingTrackingData()
    return HRSamplingTrackingData( Vector{Vector{Float64}}() , Vector{Float64}() , Vector{Float64}() , Vector{Bool}() , Vector{Int64}() , 0 ,0 ,0 , 0)
end
function Base.hcat(a::HRSamplingTrackingData,b::HRSamplingTrackingData)
    return HRSamplingTrackingData( [a.X;b.X] , [a.X_fdens;b.X_fdens] , [a.X_fcost;b.X_fcost] , [a.X_viable;b.X_viable] , [a.X_opcode;b.X_opcode] , a.fevals_expand+b.fevals_expand , a.fevals_bisect+b.fevals_bisect , a.fevals_linesearch+b.fevals_linesearch , a.fevals_fdens+b.fevals_fdens )
end


struct HRStepConfig
    bisect_max_it::Int
    bisect_rel_tol::Float64
    bisect_abs_tol::Float64
    linesearch_max_tries::Int
    linedensity_fevals::Int

    track_chords::Int
    track_bisection::Int # 1: viable samples, 2: all samples
    track_linesearch::Int  # 1: viable samples, 2: all samples
end
function HRStepConfig( ; bisect_max_it::Int=Int64(-1) , bisect_rel_tol=NaN , bisect_abs_tol=Inf , linesearch_max_tries::Int=Int64(-1) , linedensity_fevals::Int=Int64(100) ,  track_chords::Int=Int64(0) , track_bisection::Int=Int64(0) , track_linesearch::Int=Int64(0) )
    HRStepConfig(bisect_max_it,bisect_rel_tol,bisect_abs_tol,linesearch_max_tries,linedensity_fevals,track_chords,track_bisection,track_linesearch)
end


function push_hrsample!(data::HRSamplingTrackingData,x::Vector{Float64},xfdens::Float64,xfcost::Float64,xviable::Bool,opcode::Int64)
    push!(data.X,x)
    push!(data.X_fdens,xfdens)
    push!(data.X_fcost,xfcost)
    push!(data.X_viable,xviable)
    push!(data.X_opcode,opcode)
end

"""
  bisection( f , threshold , a , b , max_it , max_rel_tol )

bisection performs bisection of the function of f on the line
from a to b such that the point is viable (i.e. f value is below threshold) and close to a. We assume that
f(b) is viable, f(a) most likely is not.

  f : function which returns as a first return value a
      boolean which describes if the point is viable.

  a,b : specify a line segment. b is the point which is viable, and we try
        to find the point closest to a which is still viable
  max_it : maximum iterations for moving the viable point close to a.
           NOTE: we may require more iterations than max_it if we have to
           move close to b to find a viable point.

# Output:
- x : the discovered boundary
- num_fcalls  : number of function evaluations
- x_viable    : tested points that are viable (without x).
- x_nonviable : tested points that are non-viable.

# TEST:
function
x_result1 = XBisection( x1 -> x1>7 , 5 , 10 , 10 );
x_result2 = XBisection( x1 -> x1>7 , 5 , 1e6 , 6 );
x_result3 = XBisection( x1 -> x1>(1e6-4) , 5 , 1e6 , 5 );


Thomas Liphardt, May 2018
"""
function bisection( f , threshold::Float64 , a::Array{Float64,1} , b::Array{Float64,1} , tdata::HRSamplingTrackingData ; bisect_max_it::Int=-1 , bisect_max_rel_tol::Float64=NaN , bisect_max_abs_tol=Inf , bisect_track_samples::Int64=Int64(0))
    #----- default parameters: --------------------------------
    if(bisect_max_it<0); bisect_max_it = 16; end
    if(isnan(bisect_max_rel_tol)); bisect_max_rel_tol = 1e-3; end
    #----- end default parameters: ----------------------------

    dim = size(a,1)

    #if(nargin<4), max_it      = 10; end
    #if(nargin<4), max_it      = 8; end
    #if(nargin<5), max_rel_tol = 1e-4; end

    n_it::Int64 = 1

    # We parametrize the line via x: (x=0) is a, (x=1) is b
    f_mix = (pa::Array{Float64,1},pb::Array{Float64,1},px::Float64) -> pa + (pb-pa)*px
    xa::Float64 = 0; xb::Float64 = 1

    #x_viable    = zeros(dim,0)
    #x_nonviable = zeros(dim,0)
    found_viable::Bool = false
    viable_fval::Float64   = NaN
    x_best_viable::Float64      = NaN
    x_fval_best_viable::Float64 = NaN

    x::Float64           = NaN

    norm_ab = norm(b-a)

    while( n_it<bisect_max_it && ( (xb-xa)>bisect_max_rel_tol  || (xb-xa)*norm_ab>bisect_max_abs_tol ) )
        #fprintf('\n%10.8f - %10.8f    --> x = %10.8f ',xa,xb,(xa+xb)/2);

        # bisect
        x = (xa+xb)/2.

        # check viability:
        tdata.fevals_bisect += 1
        viable_fval = f( f_mix(a,b,x) )
        viable = viable_fval < threshold

        if(viable)
            found_viable = true
            #x_viable(:,end+1) = f_mix(a,b,x)
            #x_viable = hcat( f_mix(a,b,x) , f_mix(a,b,x) )
            # if viable, go further left:
            x_best_viable = x
            x_fval_best_viable = viable_fval
            xb = x
            if(bisect_track_samples>=1); push_hrsample!(tdata,f_mix(a,b,x_best_viable),NaN,viable_fval,true,65); end
        else
            #x_nonviable(:,end+1) = f_mix(a,b,x)
            #x_nonviable = hcat( x_nonviable , f_mix(a,b,x) )
            # if not viable, go further right:
            xa = x
            if(bisect_track_samples>=2); push_hrsample!(tdata,f_mix(a,b,x),NaN,viable_fval,false,66); end
        end

        n_it = n_it+1
        if(n_it>1e3)
            if( ~f(b) )
                error("\n[ERROR] : XBisection : Point B was non-viable!\n")
            else
                @printf("[Problem] : XBisection : Bisection did not find other point..\nReturn original b..\n")
            end
        end
    end

    # if the last actual point is not viable, then we take the xb point.
    # xb is always viable, that's an assertion.
    x_result::Array{Float64,1} = b
    if(!found_viable)
    else
        x_result = f_mix(a,b,x_best_viable)
    end

    # remove the last point that we added to x_viable:
    #x_viable(:,end) = []
    #x_viable = x_viable(:,setdiff(1:size(x_viable,2),size(x_viable,2)))
    num_fcalls = n_it

    return (found_viable,x_result,x_fval_best_viable,num_fcalls) # :: Tuple{Bool,Vector{Float64},Float64,Integer64}

end


"""
  hr_step_unif!(A::Array{Float64,2} , b::Array{Float64,1} , x::Array{Float64,2} , nsteps::Int64)
hr steps for each column of x will be performed in place
random numbers are generated by the Julia GLOBAL_RNG
"""
function hr_step_unif!(A::Array{Float64,2} , b::Array{Float64,1} , x::Array{Float64,2} , nsteps::Int64)
    (d,n) = size(x)
    for zi=1:nsteps
        rv = randn(d,n)
        #rv = broadcast( / , rv , sum(rv.^2,1) )
        ru = rand(1, n)
        (lambda_min , lambda_max , x_intersection_min , x_intersection_max) = xrayshoot( A , b , x , rv )
        x[:,:] = x_intersection_min + broadcast( * , x_intersection_max-x_intersection_min , ru )
    end
    #return xi
end

function sample_unif_hr( G::Array{Float64,2} , h::Array{Float64,1} , x::Array{Float64,1} , nsamples::Int64 , nsteps::Int64 )
    x2 = repmat(x,1,nsamples)
    hr_step_unif!(G , h , x2 , nsteps)
    return x2
end

"""
  hr_step( f_density , f_cost , threshold::Float64 , G::Array{Float64,2} , h::Array{Float64,1} , x::Array{Float64,1} , x_candidate_a::Array{Float64,1}  , x_candidate_b::Array{Float64,1} , hrconf::HRStepConfig ,  tdata::HRSamplingTrackingData  )

# Returns: ( status , x_next , f_density(x_next) , f_cost(x_next) , n_fevals )
# status:  1: full_success , 2: success, but rbs_left failed, 3: success, but rbs_right failed
          -1: line_search failed , -2 both rbs failed
"""
function hr_step( f_density , f , threshold::Float64 , G::Array{Float64,2} , h::Array{Float64,1} , x::Array{Float64,1} , x_candidate_a::Array{Float64,1}  , x_candidate_b::Array{Float64,1} , hrconf::HRStepConfig ,  tdata::HRSamplingTrackingData  )
    #default parameter
    conf_linesearch_max_tries::Int64 = 20;
    conf_linedensity_fevals::Int64   = 100;
    if(hrconf.linesearch_max_tries>0); conf_linesearch_max_tries=hrconf.linesearch_max_tries; end
    if(hrconf.linedensity_fevals>0); conf_linedensity_fevals=hrconf.linedensity_fevals; end

    rbs_result_l = xrayshoot_and_bisect(f,threshold,G,h,x,x_candidate_a , hrconf , tdata )
    rbs_result_r = xrayshoot_and_bisect(f,threshold,G,h,x,x_candidate_b , hrconf , tdata )
    # ::(bool,x,f(x),num_fcalls)

    xa = rbs_result_l[2]
    xb = rbs_result_r[2]

    if(hrconf.track_chords>0)
        if(rbs_result_l[1]); push_hrsample!(tdata,xa,NaN,rbs_result_l[3],true,17);
        elseif(hrconf.track_chords>=2); push_hrsample!(tdata,xa,NaN,rbs_result_l[3],true,18) end
        if(rbs_result_r[1]); push_hrsample!(tdata,xb,NaN,rbs_result_r[3],true,17);
        elseif(hrconf.track_chords>=2); push_hrsample!(tdata,xb,NaN,rbs_result_r[3],true,18) end
    end

    if( !rbs_result_l[1] && !rbs_result_r[1] )
        print("Warning: bisection double fail..\n")
        # bisection double fail..
        return (-2,x,NaN,rbs_result_l[4]+rbs_result_r[4])
    end

    # line search
    ls_non_uniform = !isa(f_density,Void)
    if(!ls_non_uniform)
        # use linspace over range (we could also use random points..)
        ls_points = shuffle( linspace(1.e-8,1.-(1e-8),conf_linesearch_max_tries) )
        #ls_points = rand(linesearch_max_tries)
        ls_matrix = broadcast( + , xa , (xb-xa) * ls_points' )
    else
        # sample the 1d function:
        dist_border = 1.e-8
        fd_grid = linspace(0+dist_border,1-dist_border,conf_linedensity_fevals)
        fd_xx   = broadcast( + , xa , (xb-xa) * fd_grid' )
        fd_f    = NaN * ones(conf_linedensity_fevals)
        for zi in 1:conf_linedensity_fevals
            fd_f[zi] = f_density(fd_xx[:,zi])
        end
        tdata.fevals_fdens += conf_linedensity_fevals
        # use piecewise lin. function approx, and integrate
        fd_fi_pos = ( fd_grid[1:end-1]+fd_grid[2:end] ) ./ 2
        fd_fi     = ( fd_f[1:end-1]+fd_f[2:end] ) ./ 2
        fd_cumsum = cumsum(fd_fi)
        fd_cumsum = fd_cumsum ./ last(fd_cumsum)
        # draw conf_linesearch_max_tries from this distribution:
        r_cs_vals = rand(conf_linesearch_max_tries)

        ls_points = NaN * ones(conf_linesearch_max_tries)
        for zi=1:conf_linesearch_max_tries
            ls_points[zi] = fd_fi_pos[ findfirst( x -> x>r_cs_vals[zi] , fd_cumsum ) ]
        end
        ls_matrix = broadcast( + , xa , (xb-xa) * ls_points' )
    end
    x_next_fdens = NaN
    x_next_fcost = NaN
    x_next       = NaN * ones( length(x) )
    ls_success   = false
    zi           = 0
    for zi in 1:conf_linesearch_max_tries
        x_next_fcost = f(ls_matrix[:,zi])
        if( x_next_fcost < threshold )
            x_next       = ls_matrix[:,zi]
            if(ls_non_uniform); x_next_fdens = f_density(ls_matrix[:,zi]); end # note that we did not yet compute the exact fdens at this point.. thats why we compute it here

            ls_success = true
            break
        else
            if(hrconf.track_linesearch>=2); push_hrsample!(tdata,x_next,x_next_fdens,x_next_fcost,false,34) ;end
        end
    end
    tdata.fevals_linesearch += zi

    if(ls_success)
        success_code = 1
        if(!rbs_result_l[1]); success_code = 2;
        elseif(!rbs_result_r[1]); success_code = 3;
        end
        if(hrconf.track_linesearch>=1); push_hrsample!(tdata,x_next,x_next_fdens,x_next_fcost,true,33) ;end
        return (success_code,x_next,x_next_fdens,x_next_fcost,zi +rbs_result_l[4]+rbs_result_r[4])
    else
        return (-1,x,NaN,NaN,zi+rbs_result_l[4]+rbs_result_r[4])
    end
end

"""
xrayshoot_and_bisect(f , threshold , G::Array{Float64,2} , h::Array{Float64,1} , x::Array{Float64,1} , x_cand::Array{Float64,1} , hrconf::HRStepConfig , tdata::HRSamplingTrackingData )
shoots from x into direction of x_cand, then performs bisection to detect
point x_p with f(x_p) slightly below the threshold
"""
function xrayshoot_and_bisect(f , threshold , G::Array{Float64,2} , h::Array{Float64,1} , x::Array{Float64,1} , x_cand::Array{Float64,1} , hrconf::HRStepConfig , tdata::HRSamplingTrackingData )
    rv          = x_cand - x
    norm_rv     = norm(rv)
    rv_normed   = rv/norm_rv
    (lambda_min_a , lambda_max_a , x_intersection_min_a , x_intersection_max_a) = xrayshoot( G , h , x[:,:] , rv_normed[:,:] )

    lambda_max::Float64 = lambda_max_a[1]
    # check if candidate sample is outside of polytope
    if( lambda_max < norm_rv )
        x_cand = x + (1-1.e-6) * lambda_max * rv_normed
    end

    # check if candidate sample is viable, if so, then move further out
    fevals_expansion = 1
    xc_f        = f(x_cand)
    found_point = false
    reached_boundary = false # only set true if we reach it..
    if( xc_f < threshold )
        #print("expand..\n")
        # try to find viable point by taking points further away..
        xd  = x_cand-x
        nxd = norm(xd)
        multiplier = 1
        reached_boundary = false
        while( xc_f <= threshold && !reached_boundary )
            if( 2*multiplier*nxd < lambda_max)
                multiplier    =  2*multiplier
                x_cand = x + multiplier * nxd * xd
            else
                #print("reached max_expansion..\n")
                x_cand = x + (1-(1.e-6)) * lambda_max * (xd/nxd) # little bit inside
                reached_boundary = true
            end
            fevals_expansion += 1
            tdata.fevals_expand += 1
            xc_f = f(x_cand)
        end
    end
    if(xc_f <= threshold)
        #print("expansion yielded two viable points\n")
        if(!reached_boundary); print("Warning: max expansion did not reach boundary?! That's strange.."); end
        #return (false,NaN*ones(length(x)),NaN,fevals_expansion)
        # ok, then just accept these two points..
        return (true,x_cand,xc_f,fevals_expansion)
    end

    bisect_result = bisection( f , threshold ,  x_cand , x , tdata ; bisect_max_it=hrconf.bisect_max_it , bisect_max_rel_tol=hrconf.bisect_rel_tol , bisect_max_abs_tol=hrconf.bisect_abs_tol , bisect_track_samples = hrconf.track_bisection )
    # ::(bool,x,f(x),num_fcalls)
    return (bisect_result[1],bisect_result[2],bisect_result[3],bisect_result[4]+fevals_expansion)
end
