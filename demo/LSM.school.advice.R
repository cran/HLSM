#library(HLSM)
#data(schoolsAdviceData)

niter=10
fit=LSM(Y=School9Network, senderCov=School9NodeCov, receiverCov=School9NodeCov, edgeCov=School9EdgeCov, niter=niter)


