
function reconstruct_live_samples(itd::Vector{NSIterData})
      X_live   = Vector{Array{Float64,2}}(length(itd))
      X_live_f = Vector{Array{Float64,1}}(length(itd))

      # start with final iterdata, and then go backwards..
      X_live[length(itd)]   = itd[end].X_sorted_out
      X_live_f[length(itd)] = itd[end].Xf_sorted_out

      for zi in ((length(itd)-1):-1:1)
            X_live[zi]   = X_live[zi+1]
            X_live_f[zi] = X_live_f[zi+1]
            # remove the sorted out samples


            # remove the newly generated samples of step zi:
            x_step_i  = reduce( (x,y) -> [x y] , map( x -> x.xx[:,end] , itd[zi].dd ) )
            idx_new = Vector{Int64}()
            for zj in 1:size(x_step_i,2)
                 append!(idx_new , find( all( X_live[zi] .== x_step_i[:,zj] , 1 ) ) )
            end
            X_live[zi]   = X_live[zi][ : , setdiff(collect(1:size(X_live[zi],2)),idx_new) ]
            X_live_f[zi] = X_live_f[zi][ setdiff(collect(1:size(X_live[zi],2)),idx_new) ]
            #xf_step_i = reduce( (x,y) -> [x;y] , map( x -> x.xxf[end] , itd[zi].dd ) )


            # add the sorted out samples
            X_live[zi]     = [ X_live[zi]     itd[zi].X_sorted_out   ]
            X_live_f[zi]   = [ X_live_f[zi] ; itd[zi].Xf_sorted_out  ]
      end

      return (X_live,X_live_f)
end

struct NSRunStatistics
      steps
end

function get_statistics_of_ns_run(itd::Vector{NSIterData})


end
