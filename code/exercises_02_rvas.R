# Source: https://github.com/joellembatchou/SISG2022_Association_Mapping/blob/master/code/Session07_exercises.R

library(tidyverse)
library(data.table)
library(glue)

library(BEDMatrix)
library(SKAT)
library(ACAT)

## Define paths to data files/tools
# - helpful R commands: getwd(), list.files(), file.exists()
DIR_DATA = 'data/' # or simple '.'
D1_BFILE = glue('{DIR_DATA}/rv_geno_chr1')
D1_PHENO = glue('{DIR_DATA}/rv_pheno.txt')

PLINK = 'plink2'

# Q1: Extract a subset of rare markers and store a new dataset
val_maf_max = NULL # put the value of maxium of MAF
cmd = glue(" {PLINK} --bfile {D1_BFILE} ",
    " --max-maf {val_maf_max} --maj-ref force ",
    " --make-bed --out rv_geno_chr1_subset ")
system(cmd)

## Q2: Load genetic data into R
# - Compute Minor Allele Frequencies (MAF) for each variant. What is the range of MAF values?
# - Check the missingness rate? Is it acceptable?
bfile = 'rv_geno_chr1_subset'

# genotypes are still stored in files 
# and can be loaded into RAM by small chunks 
g = BEDMatrix(bfile, simple_names = TRUE) 

# genotypes are converted to a regular R matrix and all data is in RAM
# !NB! not recommended for real large datasets
G = as.matrix(g)

str(G)

# Hint: apply colMeas function to G and specify the argument na.rm = TRUE
# - Recall the MAF formula: sum(g) / (2*N) = mean(g) / 2

# write your code to compute MAF of all columns in G 
mafs = NULL
range(mafs) 

macs = colSums(G, na.rm = TRUE)
range(macs)

## Q3: Read phenotypes into R. Check the distribution
phen = fread(D1_PHENO)
str(phen)

ggplot(phen, aes(Pheno)) + geom_histogram()

## Q4: Run single-variant association analysis by Plink
# - Make a Volcano plot (beta vs. -log10(p)). What method, Burden, SKAT or ACAT, do you expect to perform best?
cmd = glue(" {PLINK} --bfile rv_geno_chr1_subset ",
    " --pheno {D1_PHENO} --pheno-name Pheno ",
    " --glm allow-no-covars --out gwas_sv ")
system(cmd)

## read association results
f_assoc = NULL # put the path to output file from Plink here
assoc_plink = fread(f_assoc) %>% as_tibble
str(assoc_plink)

## plot effect size vs -log10(p)
ggplot(assoc_plink, aes(BETA, -log10(P))) + geom_point() + 
  geom_vline(xintercept = 0, linetype = 3)

## Q5: Run a weighted burden test

# build a burden score
wts = dbeta(mafs, 1, 25)
burden = G * wts # that is the column-wise multiplication in R
burden[1:5, 1:3]

# statistics on the burden score
score = rowSums(burden)
table(score > 1)

# burden test
fit = lm(phen$Pheno ~ score)
summary(fit)

## Q6: Run SKAT
# - Run SKAT-O with rho = 1 and rho = 0. What two tests are these?
# - RUN SKAT-O with the optimal rho parameter
null_skat = SKAT_Null_Model(phen$Pheno ~ 1 , out_type = "C")
res_skat = SKAT(G, null_skat)
str(res_skat)

# rho = 1 --> SKAT-O -> ?
val_rho = NULL # put the rho value here
res_skat_rho1 = SKAT(G, null_skat, r.corr = val_rho)
res_skat_rho1$p.value

# rho = 0 --> SKAT-O -> ?
val_rho = NULL # put the rho value here
res_skat_rho1 = SKAT(G, null_skat, r.corr = val_rho)
res_skat_rho1$p.value

# optimal rho
res_skat_o = SKAT(G, null_skat, method = "optimal.adj")
res_skat_o$p.value

## Q7: Run ACAT
# - Check ACAT vs. (corrected) MinP. Does ACAT improve over MinP?
wts_beta = dbeta(mafs, 1, 25)
wts_acat <- wts_beta * wts_beta * mafs * (1 - mafs)

ACAT(assoc_plink$P, weights = wts_acat)

# ACAT vs. MinP
min(assoc_plink$P) * nrow(assoc_plink)

## Q8: Finally, combine the 3 p-values from the three tests, Burden, SKAT-O and ACAT.


