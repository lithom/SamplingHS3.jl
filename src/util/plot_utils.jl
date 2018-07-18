


using Plots
function gr_hist_2d(xx::Array{Float64,2} ; nbins=20)
      gr()
      Plots.histogram2d( xx[1,:] , xx[2,:] , nbins=nbins )
end

function plot_heatmap_2d_function( f::Any , xlim , ylim , xres , yres )
    xygrids = xygrid(xlim,ylim,xres,yres)
    zgrid   = NaN*ones(xres,yres)
    for zi=1:length(zgrid); zgrid[zi] = f_ec([xygrids[1][zi];xygrids[2][zi]] ); end
    using Plots
    gr()
    Plots.heatmap(zgrid)
end
