#include("SamplingHS3.jl")
import SamplingHS3


A1=[eye(2);-eye(2)]; b1=[1.;1.;1.;1.];

v0 = [1 1 -1 0.1 ; 1 0 1 0.2 ]
v1 = [1 1 -1 0.1 ; 1 0 1 0.2 ]
x0 = [0;0]

lambda_min , lambda_max , x_intersection_min , x_intersection_max = SamplingHS3.xrayshoot(A1,b1,zeros(2,4),v0)
# lambda_min , lambda_max , x_intersection_min , x_intersection_max = XRayShoot(A1,b1,zeros(2,1),[[1;1],[1;0],[-1;1],[0.1;0.2]]);


## TEST BISECTION..
#include("SamplingHS3.jl")
import SamplingHS3
bs_b = [1.]
bs_a = [14.]
bs_f = x1 -> sum(x1.^2.)[1]
tracker = SamplingHS3.HRSamplingTrackingData()
bs_result = SamplingHS3.bisection( bs_f , 4. ,  bs_a , bs_b , tracker )
@test( bs_f( bs_result[2] )<4. )




using SamplingHS3



## SMALL TEST:


G_test         = [eye(2);-eye(2)]
h_test         = [10.;10;10;10]
f_test         = x1 -> sqrt.(sum(x1.^2,1))[1]
threshold_test = 4.0
x              = zeros(2)
Sigma_test     = 0.1*eye(2)

tdata  = SamplingHS3.HRSamplingTrackingData()
hrconf = SamplingHS3.HRStepConfig()
rsb_result =  SamplingHS3.xrayshoot_and_bisect(f_test , threshold_test , G_test , h_test , x , [1.;1] , hrconf , tdata )

n_test         = 10_000
#X_Collected    = zeros(2,n_test)
X_Collected    = Vector{Vector{Float64}}()
tic()
for zi in 1:n_test
    push!(X_Collected,x)
    rand_ab     = randn(2)
    hrstep_result = SamplingHS3.hr_step( Void() , f_test , threshold_test , G_test , h_test , x , x+rand_ab  , x-rand_ab , hrconf , tdata )
    x = hrstep_result[2]
end
toc()


xx_x_collected = reduce( (x,y) -> [x y] ,X_Collected)

if(false)
    using Plots
    Plots.plot(  xx_x_collected[1,1:400],xx_x_collected[2,1:400])
end

tdata   = SamplingHS3.HRSamplingTrackingData()
hrconf  = SamplingHS3.HRStepConfig(;track_chords=2,track_bisection=2,track_linesearch=2)
n_test         = 100_000
#X_Collected    = zeros(2,n_test)
X_Collected    = Vector{Vector{Float64}}()
tic()
for zi in 1:n_test
    push!(X_Collected,x)
    rand_ab     = randn(2)
    hrstep_result = SamplingHS3.hr_step( Void() , f_test , threshold_test , G_test , h_test , x , x+rand_ab  , x-rand_ab , hrconf , tdata )
    x = hrstep_result[2]
end
toc()
x_all_tracked = reduce( (x,y) -> [x y] , tdata.X)
x_all_bisect_a = x_all_tracked[ : , tdata.X_opcode.==65 ]
x_all_bisect_b = x_all_tracked[ : , tdata.X_opcode.==66 ]

if(false)
    Plots.plot()
    Plots.scatter!( x_all_bisect_a[1,1:2000] , x_all_bisect_a[2,1:2000] )
    Plots.scatter!( x_all_bisect_b[1,1:2000] , x_all_bisect_b[2,1:2000] )
end

## Test with dist2line function..
target_line = Vector{SamplingHS3.LineSegment}()
push!(target_line,SamplingHS3.LineSegment([1 2;1 1]))
push!(target_line,SamplingHS3.LineSegment([2 2;1 0]))
push!(target_line,SamplingHS3.LineSegment([2 3;0 0]))


f_test_2 = (x) -> SamplingHS3.dist2line(target_line,x[:,:])[1]
threshold_test_2 = 0.25
x = [1.;1]
tdata   = SamplingHS3.HRSamplingTrackingData()
hrconf  = SamplingHS3.HRStepConfig(;track_chords=2,track_bisection=2,track_linesearch=2)
X_Collected    = Vector{Vector{Float64}}()
tic()
for zi in 1:n_test
    push!(X_Collected,x)
    rand_ab     = randn(2)
    hrstep_result = SamplingHS3.hr_step( Void() , f_test_2 , threshold_test_2 , G_test , h_test , x , x+rand_ab  , x-rand_ab , hrconf , tdata )
    x = hrstep_result[2]
end
toc()
x_all_tracked = reduce( (x,y) -> [x y] , tdata.X)
x_all_bisect_a = x_all_tracked[ : , tdata.X_opcode.==65 ]
x_all_bisect_b = x_all_tracked[ : , tdata.X_opcode.==66 ]

if(false)
    Plots.plot()
    Plots.scatter!( x_all_bisect_a[1,1:2000] , x_all_bisect_a[2,1:2000] )
    Plots.scatter!( x_all_bisect_b[1,1:2000] , x_all_bisect_b[2,1:2000] )
end
