# generate samples of correlations from distribution induced by Bingham distribution bing(M, Z)

# random number generator of Bingham distribution
library(Directional)
library(rstiefel)
bingrnd = function(n, A, pkg='rstiefel'){
	switch(pkg,
		Directional=return(rbingham(n,A)), # based on rejection sampling
		rstiefel={
			d=dim(A)[1]; v=matrix(NA,n,d)
			v[1,]=rbing.Op(A,diag(d))[,1]
			if(n>1) for(i in 2:n) v[i,]=rbing.vector.gibbs(A,v[i-1,]) # based on Gibbs sampling
			return(v)
		},
		stop('Package not available!')
	)
}

# set random seed
seed_NO=2017
set.seed(seed_NO)
# set working directory
setwd('/Users/LANZI/Statistics/Projects/state-space/BayesCOV/code/iWishart')
# retrieve parameters
setting <- read.table('./tmp/setting.csv',header=FALSE,sep=',')
n <- setting$V1; zeta <- setting$V2
mu <- as.matrix(read.table('./tmp/mu.csv',header=FALSE,sep=','))

# compute the dimensions
sq_d=dim(mu)[1]
d=sq_d*(sq_d+1)/2

# generate bingham random vectors for each row of L
L = matrix(NA,n,d)
for(i in 1:sq_d){
	idx_i = (1+i*(i-1)/2):(i*(i+1)/2)
	L[,idx_i] = bingrnd(n,zeta*mu[i,1:i]*diag(i))
	}
# write to file
write.table(L,file='./tmp/rv_bing.csv',sep=',',col.names=FALSE,row.names=FALSE,qmethod='double')