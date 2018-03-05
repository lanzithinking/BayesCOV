# generate samples of correlations from distribution induced by Bingham distribution bing(M, Z)

seed_NO=2017
set.seed(seed_NO)

setwd('/Users/LANZI/Statistics/Projects/state-space/BayesCOV/code/iWishart')

source('./sample_Rho_bing.R')

setting <- read.table('./tmp/setting.csv',header=FALSE,sep=',')
n <- setting$V1; zeta <- setting$V2
mu <- as.matrix(read.table('./tmp/mu.csv',header=FALSE,sep=','))
rho_sub <- as.matrix(read.table('./tmp/rho_sub.csv',header=FALSE,sep=','))

Rho = sample_Rho_bing(n,zeta,mu,rho_sub)
write.table(Rho,file='./tmp/Rho.csv',sep=',',col.names=FALSE,row.names=FALSE,qmethod='double')