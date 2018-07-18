

function interp1_lin(fx::Array{Float64,1},fy::Array{Float64,1}, xx::Array{Float64,1} )
    if(!issorted(fx)); error("fx must be sorted ascending"); end
    yy = Array{Float64}(length(xx))
    for zi in 1:length(xx)
        xi = xx[zi]
        pxa = findlast( fx .<= xi  , true )
        pxb = findfirst( fx .>= xi  , true )
        if(pxa==pxb)
            yy[zi] = fy[pxa]
        else
            yy[zi] = fy[pxa] + (xi- fx[pxa] ) *  (fy[pxb]-fy[pxa]) / (fx[pxb]-fx[pxa])
        end
    end
    return yy
end



function bin2d( x::Array{Float64,2} , nbx::Int64 , nby::Int64 )
    xb = Array{Float64}(nbx,nby)
    limx = [ minimum(x[1,:]) ; maximum(x[1,:]) ]
    limy = [ minimum(x[2,:]) ; maximum(x[2,:]) ]
    bbx  = linspace(limx[1],limx[2],nbx+1)
    bby  = linspace(limy[1],limy[2],nby+1)
    for zi in 1:nbx
        for zj in 1:nby
            idx_ij = ( x[1,:].>= bbx[zi] ) .& ( x[1,:] .< bbx[zi+1] ) .& ( x[2,:] .>= bby[zj] ) .& ( x[2,:] .< bby[zj+1] )
            xb[zi,zj] = sum(idx_ij)
        end
    end
    return xb
end


# sample x[:,i] has function value y[i]. The y corresponding to y coordinates will be binned
function bin2d( x::Array{Float64,2} , y::Array{Float64,1} , nbx::Int64 , nby::Int64 )
    xb = Array{Vector{Float64}}(nbx,nby)
    limx = [ minimum(x[1,:]) ; maximum(x[1,:]) ]
    limy = [ minimum(x[2,:]) ; maximum(x[2,:]) ]
    bbx  = linspace(limx[1],limx[2],nbx+1)
    bby  = linspace(limy[1],limy[2],nby+1)
    for zi in 1:nbx
        for zj in 1:nby
            idx_ij = ( x[1,:].>= bbx[zi] ) .& ( x[1,:] .< bbx[zi+1] ) .& ( x[2,:] .>= bby[zj] ) .& ( x[2,:] .< bby[zj+1] )
            xb[zi,zj] = Vector{Float64}()
            append!(xb[zi,zj], y[idx_ij] )
        end
    end
    return xb
end

function bin2d( x::Array{Float64,2} , y::Array{Float64,1} , nbx::Int64 , nby::Int64 , op::Any)
    b2dr = bin2d(x::Array{Float64,2} , y::Array{Float64,1} , nbx::Int64 , nby::Int64)
    opx  = zeros(nby,nbx)
    for zi=1:nbx
        for zj=1:nby
            opx[zi,zj] = op(b2dr[zi,zj])
        end
    end
    return opx
end
