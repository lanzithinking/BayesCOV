# sample correlations from distribution induced by von Mises-Fisher distribution vMF(kappa, mu)

library(movMF)

sample_Rho_vmf = function(n,kappa=1,mu,rho_sub){
	rho_ind=unique(sort(rho_sub))
	L=vector('list',max(rho_ind))
	for(i in rho_ind) L[[i]]=vmfrnd(n,kappa*mu[i,1:i])
	l_rho=dim(rho_sub)[1]
	Rho=matrix(NA,n,l_rho)
	for(j in 1:l_rho){
		d=min(rho_sub[j,])
		Rho[,j]=rowSums(L[[rho_sub[j,1]]][,1:d,drop=FALSE]*L[[rho_sub[j,2]]][,1:d,drop=FALSE])
	}
	return(Rho)
}

# random number generator of Mises-Fisher distribution with one component
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