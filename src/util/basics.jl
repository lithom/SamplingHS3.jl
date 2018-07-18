
function mvnrnd(mu::Array{T,1},sigma::Array{T,2},cases::Int64) where T <: Real
      TM = chol(sigma)
      r = broadcast( + ,  TM * randn(size(TM,1),cases) , mu )
end

function sumsq(xx::Vector{T}) where T <: Real
      return sum( xx.^2 )
end
function sumsq(xx::Array{T,2}) where T <: Real
      return sum( xx.^2 , 1)
end

using Distributions.pdf
using Distributions.MvNormal

function mvnpdf(mu::Array{T,1},sigma::Array{T,2},x::Array{T,2}) where T <: Real
      return pdf(MvNormal(mu,sigma),x)
end
function mvnpdf(mu::Array{T,1},sigma::Array{T,2},x::Array{T,1}) where T <: Real
      return pdf(MvNormal(mu,sigma),x)
end



#function mvnpdf_kaputt(mu::Array{Float64,1},sigma::Array{Float64,2},x::Array{Float64,2})
#      p = size(x,2)
#      r = chol(sigma)
#      #y = (2*pi)^(-1/2) * exp.(-sumsq(  (r'\ broadcast(-,x,mu)')'  )/2) ./ prod(diag(r))
#      pdf = (2*pi)^(-p/2) * exp(-sumsq(( broadcast(-,x,mu) )/r,2)/2) / prod(diag(r))
#end

struct GaussParam
      mu::Vector{Float64}
      sigma::Array{Float64,2}
end

struct GaussMixture
      mix::Vector{GaussParam}
      weights::Vector{Float64}
end
function GaussMixture(mix::Vector{GaussParam})
      weights = ones(length(mix)) * (1./length(mix))
      return GaussMixture(mix,weights)
end

function draw_gaussmixture(gm::GaussMixture , n::Int64)
      d   = length(gm.mix[1].mu)
      x   = NaN * ones(d,n)
      xsp = sample_proportional( gm.weights , n )
      for zi=1:length(gm.weights)
            x[:,xsp.==zi] = mvnrnd( gm.mix[zi].mu,gm.mix[zi].sigma, sum(xsp.==zi) )
      end
      return x
end


function pdf_gaussmixture_nonnormalized(gm::GaussMixture , x::Array{T,2}) where T <: Real
      p = zeros(size(x,2))
      for zi=1:length(gm.mix)
            p += mvnpdf( gm.mix[zi].mu , gm.mix[zi].sigma , x )
      end
      return p
end
function pdf_gaussmixture_nonnormalized(gm::GaussMixture , x::Array{T,1}) where T <: Real
      return pdf_gaussmixture_nonnormalized(gm,x[:,:])[1]
end



function sample_proportional( lambdas , n )
#fSampleProportional( lambdas , n ), samples n times integers 1:length(lambdas), according
#to the weights provided in lambdas.
    #lambdas = lambdas[:]'
    lambdas = lambdas ./ sum(lambdas)
    csl     = cumsum(lambdas)
    s       = rand(n,1)
    li      = broadcast( >= , csl , s' )
    fmres   = findmax( li.!=0, 1 )
    li       = ind2sub( size(li) , fmres[2][:] )[1]
end

function unifsamp_box( bounds::Array{T,2} , n::Int64 ) where T <: Real
    return broadcast( + , bounds[:,1] , broadcast( * , (bounds[:,2]-bounds[:,1]) , rand(size(bounds,1),n) ) )
end


function get_bounding_box( lb::Vector{T} , ub::Vector{T} ) where T <: Real
      G = [eye(length(lb));-eye(length(lb))]
      h = [ub;-lb]
      return (G,h)
end

function xygrid(xlim,ylim,xres,yres)
      grid_x = [ xlim[1]+ (0.5+i)*((xlim[2]-xlim[1])/(xres+1)) for i=1:xres, j=1:yres ]
      grid_y = [ ylim[1]+ (0.5+j)*((ylim[2]-ylim[1])/(yres+1)) for i=1:xres, j=1:yres ]
      return (grid_x,grid_y)
end
