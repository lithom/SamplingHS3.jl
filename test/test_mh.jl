using SamplingHS3



## SMALL TEST:

G_test         = [eye(2);-eye(2)]
h_test         = [10.;10;10;10]
f_test         = x1 -> sqrt(sum(x1.^2,1))[1]
threshold_test = 4.0
x              = zeros(2)
Sigma_test     = 4.0*eye(2)
X_Collected    = zeros(2,0)

mh_tracking_data = SamplingHS3.MHSamplingTrackingData()

for zi in 1:10_000
    X_Collected = hcat( X_Collected , x )
    mhstep_result =  SamplingHS3.mh_step(f_test,threshold_test,G_test,h_test[:],x[:],Sigma_test,10,mh_tracking_data)
    x = mhstep_result[2]
end

if(false)
    using Plots
    Plots.plot(X_Collected[1,:],X_Collected[2,:])
end
