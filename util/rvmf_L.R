# generate samples of correlations from distribution induced by von Mises-Fisher distribution vMF(kappa, mu)

# random number generator of Mises-Fisher distribution with one component
library(movMF)
vmfrnd = function(n, theta){
	if(is.vector(theta)) theta=t(theta)
	l=dim(theta)[1]; d=dim(theta)[2]
	
	if(l==1){
		if(d==1) return(matrix(1,n,1))
		else return(rmovMF(n,theta))
	} else{
		if(n==l){
			output=matrix(NA,n,d)
			for(i in 1:n) output[i,]=return(rmovMF(1,theta[i,]))
		} else stop('n does not match the number of rows in theta!')
	}
}

# set random seed
seed_NO=2017
set.seed(seed_NO)
# set working directory
setwd('/Users/LANZI/Statistics/Projects/state-space/BayesCOV/code/iWishart')
# retrieve parameters
setting <- read.table('./tmp/setting.csv',header=FALSE,sep=',')
n <- setting$V1; kappa <- setting$V2
mu <- as.matrix(read.table('./tmp/mu.csv',header=FALSE,sep=','))

# compute the dimensions
sq_d=dim(mu)[1]
d=sq_d*(sq_d+1)/2

# generate vmf random vectors for each row of L
L = matrix(NA,n,d)
for(i in 1:sq_d){
	idx_i = (1+i*(i-1)/2):(i*(i+1)/2)
	L[,idx_i] = vmfrnd(n,kappa*mu[i,1:i])
	}
# write to file
write.table(L,file='./tmp/rv_rmf.csv',sep=',',col.names=FALSE,row.names=FALSE,qmethod='double')
