# This is a 1-d slice sampling

slice <- function(logf, x0, a=-Inf, b=Inf, w = 1, m = 5){ # this is used to generate one sample
	
	logy <- logf(x0) - rexp(1)
			
	# Stepping out to obtain the [L, R] range
	u = runif(1);
	L = x0 - w*u;
	R = L + w;
	v = runif(1);
	J = floor(m*v);
	K = (m-1) - J;
	
	# make sure [L, R] is within [a, b]
	L <- max(L, a)
	R <- min(R, b)
	
	
	while (J>0 && logy< logf(L)) {
		L = L - w;
		L <- max(L, a)
		J = J - 1;
	}
	#browser()
	
	while (K>0 && logy<logf(R)) {
		R = R+w;
		R <- min(R, b)
		K = K-1;
	}
	
	# Shrinkage to obtain a sample
	u = runif(1);
	x1 <- L + u*(R-L);
	while (logy > logf(x1)) {
		if (x1 < x0) L <- x1
		else R <- x1
		
		u = runif(1);
		x1 = L + u*(R-L);
	}
	
	return(x1)
}
