
using OsmoticShockModel


exp_bounds = [ 0 3 ; -1 3 ; -1 2.5 ]'

fparam = (pp) -> OsmoticShockModel.SimData( [0,1.0,2.0,2.5,3.0,6.0,7.0,8.0,12.0,14.0] , exp(pp[1]) , exp(pp[2]),exp(pp[3]),exp(pp[4]),exp(pp[5]),exp(pp[6]) )
fpcomp =  (pp) -> fparam( [-1;0.;pp[1];pp[2];3;pp[3]] )


# simulate reference trajectory:
xr_ref_compact = [.5,2.5,1]
sd_ref = fpcomp(xr_ref_compact)
#ed_ref = OsmoticShockModel.create_oshock_exp_data( sd_ref , collect(linspace(0,10,100)) , 0.005 , 0.02 )
ed_ref = OsmoticShockModel.create_oshock_exp_data( sd_ref , collect(linspace(0,10,100)) , 0.001 , 0.01 )
plot_exp(ed_ref)

ssr_a = OsmoticShockModel.eval_param2ssr( ed_ref , fpcomp([2.;2.;1.]) , 0.001 )

xi = [2,2,1]
for zi=1:10000

end
