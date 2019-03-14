ALLOWED <- function(genotypes,DAG,x){
  # This function determines the feasibility of a set of genotypes according to a given DAG of restrictions.
  # genotypes: the set of genotypes to be analyzed.
  # DAG: matrix representing the DAG of restrictions. (nrow>=1 , ncol=2), the last row is always c(0,0)). See step 4 (data preprocessing in the README file)
  # x: the number of mutations considered.
  vec<-matrix(1,nrow=(2^x),ncol=1)
  D<-dim(DAG)[1]
  if (D>1){
    for (i in 1:(2^x)){
      for (j in 1:(D-1)){
        xx<-DAG[j,1]
        yy<-DAG[j,2]
        if ((genotypes[i,yy]==1)&&(genotypes[i,xx]==0)){vec[i]<-0}## if the mutation ordering is not respected the genotype is labelled as infeasible.
      }
    }
  }
  return(vec)
}
