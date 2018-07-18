# d x 2 matrix
type LineSegment
    ab::Array{Float64,2}
end

function dist2line(line::Vector{LineSegment},x::Array{Float64,2})
    dd = zeros(size(line,1),size(x,2))
    for i=1:size(line,1)
        dd[i,:] = dist2LS(line[i],x)
    end
    #print(dd)
    d  = minimum(dd,1)
    return d[:]
end

function dist2LS(ls::LineSegment,x::Array{Float64,2})
    u = ls.ab[:,2] - ls.ab[:,1]
    v = broadcast( - , x , ls.ab[:,1] )

    uv = RowVector(u) * v
    dec = uv / (norm(u)^2)

    p = broadcast( + , ls.ab[:,1] , broadcast( * , dec , u  ) )
    d = zeros(1,size(x,2))

    #println(u)
    #println(v)
    #println(uv)
    #println(dec)
    #println(p)

    for i=1:size(x,2)
        if (dec[i]<0)
            d[i] = norm(ls.ab[:,1] - x[:,i])
        elseif (dec[i]>1)
            d[i] = norm(ls.ab[:,2] - x[:,i])
        else
                d[i] = norm( p[:,i] - x[:,i] )
        end
    end
    return d[:]
end
