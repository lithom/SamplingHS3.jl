
using Distributions

# simulates the shell compression process, based on the number of sorted out samples in the ns iterations
# result is relative to volume of original polytope
function sim_shell_compression( itd::Vector{NSIterData} ; repetitions::Int64=-1)
      # default value
      if(repetitions<0); repetitions = 1000; end

      n_a = map( xi -> xi.n_live_samples_a , itd)
      n_b = map( xi -> xi.n_live_samples_b , itd)

      log_xv      = NaN * ones(repetitions,length(n_a)+1)
      log_xv[:,1] = 0.0

      x_betas = NaN * ones(repetitions,length(n_a))

      # draw betas..
      for zi in 1:length(n_a)

            # k = sorted out, n = samples before sorting out
            # then shrinking is  1-beta( k , n+1-k )
            x_betas[:,zi]  = ( 1.0 - rand(Distributions.Beta( n_a[zi]-n_b[zi] , n_b[zi]+1) , repetitions ) )
            log_xv[:,zi+1] = log_xv[:,zi] + log.(x_betas[:,zi])
      end
      x_v    = exp.(log_xv)
      x_mean = mean( x_v , 1 )
      return ( x_mean , x_v , log_xv )
end

"""
  eval_integral( volume::Float64, itd::Vector{NSIterData} ; repetitions::Int64=-1)

returns (Float64,Vector{Float64},Array{Float64,2})
where [1] is the estimated mean, [2] is the sampled integration results
and [3] is the (sampled) integrated parts of the ns iterations

evaluates the integral based on NS data, for the pdf with pdf(x) = exp( - fcost(x))
volume is the volume of the sampled space
"""
function eval_integral( volume::Float64, itd::Vector{NSIterData} ; repetitions::Int64=-1)
      ssc_1 = sim_shell_compression(itd;repetitions=repetitions)

      ssc = volume * ssc_1[2]

      # consider sorted out ns samples as mc sample of the shell
      iv_parts = Array{Float64}(size(ssc,1),length(itd))
      for zi=1:length(itd)
            mc_xf        = mean( exp.( - itd[zi].Xf_sorted_out )  )
            iv_parts[:,zi] = ( ssc[:,zi]-ssc[:,zi+1] ) * mc_xf
      end

      iv_sum = sum(iv_parts,2)
      return (mean(iv_sum),iv_sum,iv_parts)
end

"""
  resample_ns_data( itd::Vector{NSIterData} )

returns Array{Float64,2}
draws a sample of the target distribution, given the output of a ns run
"""
function resample_ns_data( itd::Vector{NSIterData} )
      int_result = eval_integral( 1.0 , itd )
      iv_parts   = int_result[3]

      # now draw proportional to iv_parts from the sorted out samples of the iterations..
      X   = zeros(size(itd[1].X_sorted_out,1),0)
      X_p = zeros(0)
      for zi in 1:length(itd)
            X   = [ X     itd[zi].X_sorted_out ]
            X_p = [ X_p ; ones(size(itd[zi].X_sorted_out,2)) * (mean(iv_parts[:,zi]) / size(itd[zi].X_sorted_out,2) )  ]
      end

      # sample until we hit collision:
      sampled_indeces = Vector{Int64}()
      no_collision = true
      while(no_collision)
            si = SamplingHS3.sample_proportional( X_p , 4 )
            append!( sampled_indeces , si )
            no_collision = length(unique(sampled_indeces))==length(sampled_indeces)
      end
      return ( X[ : ,  sampled_indeces] , X , X_p  )
end
