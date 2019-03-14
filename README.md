# predictability_of_cancer_evolution
The following is a pipeline for estimating the predictability of cancer evolution based either on i) Conjunctive Bayesian Networks (CBNs) or ii) Fitness Landscapes:

# i) Based on Conjunctive Bayesian Networks:
## Step 0: Downloading the CT-CBN software.
CBN model has been developed by Prof. Niko Beerenwinkel group and the software is free to download via the link below:
https://www.bsse.ethz.ch/cbg/software/ct-cbn.html.  
Moreover, the R functions in the step 5 require to install the Bioconductor package "OncoSimulR" (see https://bioconductor.org/packages/release/bioc/html/OncoSimulR.html). 

## Step 1: Preparing the genotype file.
Each line in a genotype file represents a genotype, a binary vector of a given length, each element of which correponds to a given mutation, and is 1 if mutation exists and zero otherwise. The first line of the genotype file is a single number indicating the number of genotypes existing in the genotype file. Furthermore, note that the first column of the genotype file is all one (See an example of a genotype file in DATA/genotype.txt in this repository) 

## Step 2: Generating an initial DAG of restrictions using CT-CBN.
In this step, starting from an empty poset, using CT-CBN we generate an initial DAG of restrictions (see the ReadMe file of the CT-CBN (https://www.bsse.ethz.ch/cbg/software/ct-cbn.html.) for more details)

## Step 3: Generating the final DAG of restrictions and Lambda values using H-CBN.
In this step, starting from the DAG of restrictions generated in step 2, using H-CBN (with 10000 steps of simulated annealing and T=1), we generate the final DAG of restrictions (see the ReadMe file of the CT-CBN (https://www.bsse.ethz.ch/cbg/software/ct-cbn.html.) for more details).

## Step 4: Data Preprocessing.
Make sure to remove the first and the last line of the final DAG file (with .poset extension).   
Moreover, it is necessary to add a line with two zeros "0 0" to the end of the DAG file. 

## Step 5: Quantifying the predictability.
Based on the DAG file (with .poset extension) and the LAMBDA file (with .lambda extension), using "CODE/PREDICTABILITY_CBN.R" function, which depends on the "CODE/ALLOWED.R" function, predictability can be computed.

### Notes on the PREDICTABILITY_CBN.R function
#### Inputs
Inputs of the PREDICTABILITY_CBN.R function are as follows:  
i) DAG: a matrix representing the DAG of restrictions. (nrow>=1 , ncol=2), the last row is always c(0,0)). See step 4 (data preprocessing)  
ii) LAMBDA: matrix of Lambda values produced by the CBN model (nrow=(x+1),ncol=1), the first row always equals 1.  
iii) x: number of mutations considered.  
#### Output
The output of the PREDICTABILITY_CBN.R function is simply a number between 0 and 1, which is the estimated predictability.
#### Example
DAG<-read.table('DATA/EXAMPLE.poset')
LAMBDA<-read.table('DATA/EXAMPLE.lambda')
x<-4
pred<-PREDICTABILITY_CBN(DAG,LAMBDA,x)% the output equals 0.59.





# ii) Based on Fitness landscapes:

## Step 1: preparing the fitness alndscape.
By assigning a fitness value to each of the genotypes in a genotype space comprising 2^x genotypes (x is the number of genes (or mutations)), we will have a fitness landscape that can be used for quantifying the predictability of cancer evolution. 
In our simulations, we used binary genotypes of length 7, which results in a genotype space comprising 2^7=128 total genotypes (see genotypes.rds in this repository). We used 100 representable fitness landscape (see FitnessLandcape_Representable.rds in this repository) and 111 non-representable ones (see FitnessLandcape_NON_Representable.rds in this repository). Each row of the above .rds files corresponds to a given genotype (in the genotypes.rds file) and each column corresponds to a given fitness landscape. 

## Step 2: Quantifying the predictability.
Based on the fitness landscape file in step 1, we can quantify the predictability of evolution using PREDICTABILITY_SSWM.R in this folder.






