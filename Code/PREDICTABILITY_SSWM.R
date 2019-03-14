PREDICTABILITY_SSWM<-function(FITNESS,x){
    
  #FITNESS: fitness vector of length 2^x, each corresponding to a given genotype.
  #x: number of mutations considered.
  library('gtools')
  
  ### Step1: genotypes
  #source("https://bioconductor.org/biocLite.R")
  #biocLite("OncoSimulR")
  genotypes=OncoSimulR:::generate_matrix_genotypes(x)## generates the genotype space[requires "OncoSimulR" package--two lines above]
  indx<-matrix(0,nrow=2^x,ncol=1)## indexing the genotypes for easier retrival
  for (k in 1:(2^x)){for (j in 1:x){indx[k,1]=indx[k,1]+2^(j-1)*genotypes[k,j]}}
  
  ### Step2: Pathway Probabilities
  PERM<-permutations(x,x)## all x! possible permutations (mutational pathways)
  Prob<-numeric(dim(PERM)[1])## pathway probabilities
  
  
  TOT<-0 ##(the normalization factor)
  for (i in 1:dim(PERM)[1]){# for each pathway
    TEMP1=1;# temporarily stores the pathway probability [later needs to be normalized]
    vec=PERM[i,]#pathway [specific permutation of the original vector]
    GENO=matrix(0,nrow=(x+1),ncol=x)# vector of length (x+1) storing the genotypes associated with each step of the pathway
    for (j in 1:x){for (k in (j+1):(x+1)){GENO[k,(vec[j])]=1}}
    GENO_indx=matrix(0,nrow=(x+1),ncol=1)# storing the indices of the (x+1) genotypes in the genotype space.
    for (j in 1:(x+1)){for (k in 1:x){GENO_indx[j,1]=GENO_indx[j,1]+2^(k-1)*GENO[j,k]}}
    fitness=matrix(0,nrow=(x+1),ncol=1)# fitness vector for the (x+1) genotypes associated with the given pathway
    for (j in 1:(x+1)){fitness[j]=FITNESS[which(indx==GENO_indx[j])]}# retrieving the fitness from the global fitness vector
    flag=0;
    for (j in 2:(x+1)){if (fitness[j]<fitness[(j-1)]){flag=1}}# if the fitness monotonically increases along the pathway, flag remains as 0, otherwise it will become 1 
    if (flag==0){# if flag remains zero (i.e. pathway is accessible)
      for (j in 1:x){
        SN=which(GENO[j,]==0)# possible remaining mutations in the j-th step 
        N=length(SN)
        S=fitness[(j+1)]-fitness[j]# slective coefficient of the j-th step [the numerator of the equation (7) in the main text]
        T=0;
        for (k in 1:N){# checking the genotypes belonging to the exit set
          ggeno=GENO[j,]
          ggeno[(SN[k])]=1
          ggeno_indx=0
          for (l in 1:x){ggeno_indx=ggeno_indx+2^(l-1)*ggeno[l]}
          fitness2=FITNESS[which(indx==ggeno_indx)]
          S1=fitness2-fitness[j]
          if (S1>0){T=T+S1}# sum of the selecive coefficient of the genotypes in the exit set [the denominator of the equation (7) in the main text]
        }
        TEMP1=TEMP1*(S/T)#the product in equation (7)
      }
      Prob[i]=TEMP1# probability of the i-th pathway
    }
  }
  GG=sum(Prob,na.rm=TRUE)#normalization factor (equation 8 in the main text)
  for (i in 1:dim(PERM)[1]){
    Prob[i]=Prob[i]/GG
    if (sum(Prob[i],na.rm=TRUE)>0){TOT=TOT-Prob[i]*log(Prob[i])}#entropy
  }
  Pred=1-(TOT/log(factorial(x)))#predictability
  return(Pred)
}
