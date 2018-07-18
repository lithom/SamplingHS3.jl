
struct StopCriterion
    min_time::Float64
    max_time::Float64

    min_fevals::Int64
    max_fevals::Int64
end

function StopCriterion( ; min_time=Inf , max_time=4. , min_fevals=1000_000_000_000 , max_fevals=1000_000_000_000_000  )
    StopCriterion( min_time , max_time , min_fevals , max_fevals )
end
