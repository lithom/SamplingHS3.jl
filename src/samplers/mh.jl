
mutable struct MHSamplingTrackingData
    X::Vector{Vector{Float64}}
    X_fdens::Vector{Float64}
    X_fcost::Vector{Float64}
    X_viable::Vector{Bool}
    ratios_polytope_hits::Vector{Float64}
    fevals_density::Int64
    fevals_cost::Int64
    successes::Int64
end

function MHSamplingTrackingData()
    return MHSamplingTrackingData( Vector{Vector{Float64}}() , Vector{Float64}() , Vector{Float64}() , Vector{Bool}(), Vector{Float64}() , 0 , 0 , 0 )
end


"""
  mh_step( f , threshold , G::Array{Float64,2} , h::Array{Float64,1} , x::Array{Float64,1} , Sigma::Array{Float64,2} , max_tries::Int64  , tdata::MHSamplingTrackingData )

computes the mh step for the uniform distribution

"""
function mh_step( f , threshold::Float64 , G::Array{Float64,2} , h::Array{Float64,1} , x::Array{Float64,1} , Sigma::Array{Float64,2} , max_tries::Int64  , tdata::MHSamplingTrackingData )
    return mh_step( Void() , f , threshold , G::Array{Float64,2} , h::Array{Float64,1} , x::Array{Float64,1} , NaN , Sigma::Array{Float64,2} , max_tries::Int64  , tdata::MHSamplingTrackingData  )
end

"""
function mhstep( f , threshold , G::Array{Float64,2} , hG::Array{Float64,1} , x::Array{Float64,1} , x_density::Float64 , Sigma::Array{Float64,2} , max_tries::Integer64 )
[ x_next , x_next_cost , X_NonViableInPolytope , num_rejected_jumps , num_fevals , accepted ] = MHStep_01( f , threshold , G , h , x , Sigma , max_tries )

 performs a step of metropolis hastings, i.e. it generates proposal
 samples according to Sigma, and then accepts with probability 0 if
 outside the viable space, and accepts with probability (f_density(x_next) / x_f) if inside the
 viable space.

RETURNS: ::Tuple{ Bool , Array{Float64,1} , Float64  , Float64   ,  Integer64 }
                  succ ,       x    ,       fdens(x) ,  fcost(x) ,
"""
function mh_step( f_density , f , threshold , G::Array{Float64,2} , h::Array{Float64,1} , x::Array{Float64,1} , x_density::Float64 , Sigma::Array{Float64,2} , max_tries::Int64  , tdata::MHSamplingTrackingData )

    #dim = numel(x);
    factor_TriesPolytope = 100

    accepted = false
    n_tries  = 0

    # sample the whole batch.. :)
    x_c = mvnrnd(x,Sigma, factor_TriesPolytope * max_tries )

    # sort out the ones which are not in polytope..
    x_c_inPolytope = all( broadcast( - ,G*x_c,h) .< 0 , 1 )

    x_c = x_c[:,x_c_inPolytope[:]]
    #x_c = x_c[:,1:min(max_tries,size(x_c,2))]

    push!( tdata.ratios_polytope_hits , size(x_c,2)*1.0 / (factor_TriesPolytope * max_tries) )

    #num_rejected_jumps = 0max_tries - sum(x_c_inPolytope)

    num_fevals_fdensity = 0
    num_fevals_fcost    = 0

    zi = 1

    x_c_i      = NaN * ones( length(x) )
    x_c_i_cost = NaN

    # check if we sample non-unif.
    nonunif_sampling = !isa( f_density , Void )
    xfdens_next::Float64 = NaN # holds the value of the density function at x_next

    while( ~accepted && zi <= min( size(x_c,2) , max_tries ) )
        x_c_i = x_c[:,zi]
        # if nonunif. then first find next mcmc step
        if(nonunif_sampling)
            xfdens_next = f_density(x_c_i)
            num_fevals_fdensity += 1
            mcmc_step = rand() < ( xfdens_next / x_density )
            #print(" $(xfdens_next) / $(x_density) = $(xfdens_next/x_density) \n")
            if(!mcmc_step)
                zi = zi + 1
                continue
            end
        end

        num_fevals_fcost = num_fevals_fcost + 1
        x_c_i_cost   = f(x_c_i)

        accepted     = x_c_i_cost <= threshold

        zi = zi + 1

        push!(tdata.X   , x_c_i)
        push!(tdata.X_fdens , xfdens_next)
        push!(tdata.X_fcost , x_c_i_cost)
        push!(tdata.X_viable , accepted)
    end

    if(~accepted)
        #X_NonViableInPolytope = x_c(:,1:end)
        #num_rejected_jumps = zi-1
        tdata.fevals_cost     += num_fevals_fcost
        tdata.fevals_density  += num_fevals_fdensity

        x_next             = x
        x_next_cost        = NaN
        return (false,x_next,xfdens_next,x_next_cost,num_fevals_fdensity,num_fevals_fcost) # ::Tuple{ Bool , Array{Float64,1} , Float64 , Integer64 }
    else
        #X_NonViableInPolytope = x_c(:,1:zi-2);
        #num_rejected_jumps = num_rejected_jumps + zi - 1
        tdata.fevals_cost     += num_fevals_fcost
        tdata.fevals_density  += num_fevals_fdensity
        tdata.successes += 1

        x_next             = x_c_i
        x_next_cost        = x_c_i_cost
        return (true,x_next,xfdens_next,x_next_cost,num_fevals_fdensity,num_fevals_fcost) # ::Tuple{ Bool , Array{Float64,1} , Float64 ,  Integer64 }
    end
end
