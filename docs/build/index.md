
#MCMC methods

<a id='SamplingHS3.mh_step' href='#SamplingHS3.mh_step'>#</a>
**`SamplingHS3.mh_step`** &mdash; *Function*.



mh_step( f , threshold , G::Array{Float64,2} , h::Array{Float64,1} , x::Array{Float64,1} , Sigma::Array{Float64,2} , max_tries::Int64  , tdata::MHSamplingTrackingData )

computes the mh step for the uniform distribution


function mhstep( f , threshold , G::Array{Float64,2} , hG::Array{Float64,1} , x::Array{Float64,1} , x_density::Float64 , Sigma::Array{Float64,2} , max_tries::Integer64 ) [ x_next , x_next_cost , X_NonViableInPolytope , num_rejected_jumps , num_fevals , accepted ] = MHStep_01( f , threshold , G , h , x , Sigma , max_tries )

performs a step of metropolis hastings, i.e. it generates proposal  samples according to Sigma, and then accepts with probability 0 if  outside the viable space, and accepts with probability (f_density(x_next) / x_f) if inside the  viable space.

RETURNS: ::Tuple{ Bool , Array{Float64,1} , Float64  , Float64   ,  Integer64 }                   succ ,       x    ,       fdens(x) ,  fcost(x) ,

<a id='SamplingHS3.hr_step' href='#SamplingHS3.hr_step'>#</a>
**`SamplingHS3.hr_step`** &mdash; *Function*.



hr_step( f_density , f_cost , threshold::Float64 , G::Array{Float64,2} , h::Array{Float64,1} , x::Array{Float64,1} , x_candidate_a::Array{Float64,1}  , x_candidate_b::Array{Float64,1} , hrconf::HRStepConfig ,  tdata::HRSamplingTrackingData  )

**Returns: ( status , x_next , f_density(x_next) , f_cost(x_next) , n_fevals )**

**status:  1: full_success , 2: success, but rbs_left failed, 3: success, but rbs_right failed**

```
      -1: line_search failed , -2 both rbs failed
```


#Nested sampling

<a id='SamplingHS3.nested_sampling' href='#SamplingHS3.nested_sampling'>#</a>
**`SamplingHS3.nested_sampling`** &mdash; *Function*.



nested_sampling( X_init::Array{Float64,2} , f::Any , G::Array{Float64,2} , h::Array{Float64,1} , conf::NSConfig ; X_f_init::Vector{Float64}=[NaN] )   returns NestedSamplingResult

this function runs the nested sampling procedure


<a id='Example-usage-1'></a>

# Example usage


Coming soon..

