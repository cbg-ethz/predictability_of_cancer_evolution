PREDICTABILITY_CBN<-function(DAG,LAMBDA,x){
  #DAG: matrix representing the DAG of restrictions. (nrow>=1 , ncol=2), the last row is always c(0,0)). See step 4 (data preprocessing in the README file)
  #LAMBDA: matrix of Lambda values produced by the CBN model (nrow=(x+1),ncol=1), the first row always equals 1.
  #x: number of mutations considered.
  library('gtools')
  
  ### Step1: genotypes
  #source("https://bioconductor.org/biocLite.R")
  #biocLite("OncoSimulR")
  genotypes=OncoSimulR:::generate_matrix_genotypes(x)## generates the genotype space[requires "OncoSimulR" package--two lines above]
  indx<-matrix(0,nrow=2^x,ncol=1)## indexing the genotypes for easier retrival
  for (k in 1:(2^x)){for (j in 1:x){indx[k,1]=indx[k,1]+2^(j-1)*genotypes[k,j]}}
  
  ### Step2: Allowed genotypes according to the DAG of restrictions (DAG)
  allowed_set<-ALLOWED(genotypes,DAG,x)
  
  ### Step3: Pathway Probabilities
  PERM<-permutations(x,x)## all x! possible permutations (mutational pathways)
  Prob<-numeric(dim(PERM)[1])## pathway probabilities
  TOT<-0 ## (the normalization factor)
  for (i1 in 1:dim(PERM)[1]){
    TEMP1<-1 #temporarily stores the probability (later needs to be normalized)
    vec<-PERM[i1,]# (i1)-th pathway
    GENO<-matrix(0,nrow=(x+1),ncol=x) # constructing the possible x+1 genotypes, each corresponding to a given step of the given mutational pathway.  
    for (j1 in 1:x){for (k1 in (j1+1):(x+1)){GENO[k1,(vec[j1])]<-1}}
    
    flag=1;
    for (j1 in 1:(x+1)){
      index=0 # index of the (j1-th) genotype in the mutational pathway.
      for (k1 in 1:x){index<-index+2^(k1-1)*GENO[j1,k1]} 
      final_index<-which(indx==index) #finding the index of the genotype in the genotype space
      if (allowed_set[final_index]==0){flag<-0}# to check whether each genotype visited is allowed. 
    }
    
    for (j1 in 1:x){
      set1<-which(GENO[j1,]==1)
      set2<-which(GENO[j1+1,]==1)
      indx_lambda<-setdiff(set2,set1)# the index of the (j1-th) genotype in the pathway to retrieve the corresponding Lambda value
      SN<-which(GENO[j1,]==0)# the set of genotypes with one additional mutation than the current genotype
      ###################################################
      SN2<-numeric(length(SN))
      for (kaka in 1:length(SN)){
        TEMP_index=0;
        for (kk in 1:x){TEMP_index<-TEMP_index+2^(kk-1)*GENO[j1,kk]}
        TEMP_index<-TEMP_index+2^(SN[kaka]-1)
        FINAL_index<-which(indx==TEMP_index);
        
        SN2[kaka]<-allowed_set[FINAL_index]
      }
      SNN<-SN[which(SN2==1)]# THE EXIT SET: the set of (allowed) genotypes with one additional mutation than the current genotype
      T<-sum(LAMBDA[(SNN+1),1]);# sum of the lambdas of the exit set [The denominator of the equation (10) in the main text]
      ###################################################
      S<-LAMBDA[(indx_lambda+1),1]# the lambda of the (j1-th mutation) [The numerator of the equation (10) in the main text]
      TEMP1<-TEMP1*(S/T)# The multiplication in the equation (10) in the main text
    }
    if (flag==0){TEMP1<-0}# If the pathway is infeasible, its probability will be zero.
    Prob[i1]<-TEMP1 #pathway probability
  }
  
  ### Step4: Quantifying the predictability
  GG=sum(Prob,na.rm=TRUE)#normalization factor (equation 11 in the main text)
  for (i1 in 1:dim(PERM)[1]){
    Prob[i1]=Prob[i1]/GG #Normalizing pathway probability such that the sum of all x! pathway probabilities become 1. 
    if (sum(Prob[i1],na.rm=TRUE)>0){TOT=TOT-Prob[i1]*log(Prob[i1])} #Computing the entropy of cancer progression for the given fitness landscape 
  }  
  
  Pred=1-(TOT/log(factorial(x)))## computing the predictability
  return(Pred)
}
