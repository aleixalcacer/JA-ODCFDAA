library(R.matlab)
library(MASS)

da=readMat('Fl15.mat')



cxyz=c();

   

g=c(rep(1,12), rep(2,6),rep(3,10));



nrep=28;

for (i in 1:nrep){
    
		da$Fl[,,1]

    r=which(da$Fl[,i,]$x != 0)  
    s=which(da$Fl[,i,]$y != 0)
    t=which(da$Fl[,i,]$z != 0)
		
		x=da$Fl[,i,]$x
		y=da$Fl[,i,]$y
		z=da$Fl[,i,]$z
    
    cxyz=rbind(cxyz, cbind(t(x[r]) ,t(y[s]) ,t(z[t]) )   )
    
}




library(archetypes)
source("LstepArchetypesMod.R")

X=cxyz

#Find archetypes and screeplot 
library(DDoutlier)
nada=4
norep=20
set.seed(1234)
ai=LstepArchetypesMod(X, 1:nada, nrep=norep)

#Draw the RSS graph in order to find the elbow (number of initial archetypes)
screeplot(ai)

#Choose the number of initial archetypes
e=2

vi=2
vf=2
#AA+ k-NN 

for (i in vi:vf)
{
pei=rep(0,dim(X)[1])
for (j in e:nada){

a=bestModel(ai[[j]])

pei=pei+KNN_SUM(a$alphas, k=i) 

}
}



#####Converting into binary labels: 
vi=2
vf=2
h=1
pii=matrix(0,nrow=(vf-vi+1),ncol=dim(X)[1])
for (i in vi:vf){

pei=rep(0,dim(X)[1])
for (j in e:nada){

a=bestModel(ai[[j]])

pei=pei+KNN_SUM(a$alphas, k=i) 

}
b=boxplot(pei)
pii[h,]=as.integer(is.element(pei,b$out))
h=h+1
}

#voting: It will be an outlier if it is voted more than 50%
#Draw a table with the number of outliers 
#table(Predictions = as.integer(apply(pii,2,sum)>(((vf-vi+1))*0.5)), TrueLabels = label)

as.integer(apply(pii,2,sum)>(((vf-vi+1))*0.5))



#POLR
m <- polr(ordered(g) ~ pei, Hess=TRUE)
summary(m)
exp(coef(m))
 


