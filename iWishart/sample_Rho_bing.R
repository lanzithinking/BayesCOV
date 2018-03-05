# sample correlations from distribution induced by Bingham distribution bing(M=I, Z=diag(0,...0,zeta))
# rbingham(n,A=M*Z*t(M))

library(Directional)
library(rstiefel)

sample_Rho_bing = function(n,zeta=1,mu,rho_sub){
	rho_ind=unique(sort(rho_sub))
	L=vector('list',max(rho_ind))
	for(i in rho_ind) L[[i]]=bingrnd(n,zeta*mu[i,1:i]*diag(i))
	l_rho=dim(rho_sub)[1]
	Rho=matrix(NA,n,l_rho)
	for(j in 1:l_rho){
		d=min(rho_sub[j,])
		Rho[,j]=rowSums(L[[rho_sub[j,1]]][,1:d,drop=FALSE]*L[[rho_sub[j,2]]][,1:d,drop=FALSE])
	}
	return(Rho)
}

# random number generator of Bingham distribution
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