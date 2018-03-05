# BayesCOV
Flexible Bayesian Dynamic Modeling of Covariance and Correlation Matrices

Modeling correlation (and covariance) matrices is a challenging problem due to the large dimensionality and positive-definiteness constraint.
In this paper, we propose a novel Bayesian framework based on modeling the correlations as products of vectors on unit spheres. 
The covariance matrix is then modeled by using its decomposition into correlation and variance matrices. 
This approach allows us to induce flexible prior distributions for covariance matrices (that go beyond the commonly used inverse-Wishart prior) by specifying a wide range of distributions on spheres (e.g. the squared-Dirichlet distribution). 
For modeling real-life spatio-temporal processes with complex dependence structures, we extend our method to dynamic cases and introduce unit-vector Gaussian process priors in order to capture the evolution of correlation among multiple time series. 
To handle the intractability of the resulting posterior, we introduce the adaptive $\Delta$-Spherical Hamiltonian Monte Carlo. 
Using an example of Normal-Inverse-Wishart problem, a simulated periodic process, and an analysis of local field potential (LFP) 
data (collected from the hippocampus of rats performing a complex sequence memory task), we demonstrate the validity and effectiveness of our proposed framework for (dynamic) modeling covariance and correlation matrices.

See link [link](https://arxiv.org/abs/1711.02869) to the paper.