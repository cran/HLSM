#library(HLSM)
#data(schoolsAdviceData)

niter=100

#Single Network LSM
fit=LSM(Y=School9network, senderCov=School9NodeCov, receiverCov=School9NodeCov, edgeCov=School9EdgeCov, niter=niter)

#Random Regression Coefficients#


random.fit = HLSMrandomEF(Y = ps.advice.mat, senderCov = ps.node.df, receiverCov=ps.node.df, edgeCov=ps.edge.df, dd = 2,niter = niter)
summary(random.fit)
names(random.fit)

fixed.fit = HLSMfixedEF(Y = ps.advice.mat, senderCov = ps.node.df, receiverCov=ps.node.df, edgeCov=ps.edge.df, dd = 2,niter = niter)

summary(fixed.fit)
names(fixed.fit)

