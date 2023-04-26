library(tidyverse)
library(data.table)
library(glue)

library(BEDMatrix)
library(SKAT)
library(ACAT)

## Paths to data files/tools
DIR_MAIN = '~/git/variani/UNAM2023_Association_Mapping/'
DIR_DATA = glue('{DIR_MAIN}/data/')
D1_BFILE = glue('{DIR_DATA}/rv_geno_chr1')
D1_PHENO = glue('{DIR_DATA}/rv_pheno.txt')

DIR_TOOLS = glue('{DIR_MAIN}/tools/')
PLINK = glue('{DIR_TOOLS}/plink2')

# Extract a subset of rare markers
cmd = glue(" {PLINK} --bfile {D1_BFILE} ",
    " --max-maf 0.01 --maj-ref force ",
    " --make-bed --out rv_geno_chr1_subset ")
system(cmd)

## load genotypes into R
bfile = 'rv_geno_chr1_subset'

# genotypes are still stored in files 
# and can be loaded into RAM by small chunks 
g = BEDMatrix(bfile, simple_names = TRUE) 

# genotypes are converted to a regular R matrix and all data is in RAM
# !NB! not recommended for real large datasets
G = as.matrix(g)

str(G)

mafs = colMeans(G, na.rm = TRUE) / 2
range(mafs) 

macs = colSums(G, na.rm = TRUE)
range(macs)

## read phenotypes
phen = fread(D1_PHENO)
str(phen)

ggplot(phen, aes(Pheno)) + geom_histogram()

## single-variant association analysis by Plink
cmd = glue(" {PLINK} --bfile rv_geno_chr1_subset ",
    " --pheno {D1_PHENO} --pheno-name Pheno ",
    " --glm allow-no-covars --out gwas_sv ")
system(cmd)

## read association results
assoc_plink = fread('gwas_sv.Pheno.glm.linear') %>% as_tibble
str(assoc_plink)

## plot effect size vs -log10(p)
ggplot(assoc_plink, aes(BETA, -log10(P))) + geom_point() + 
  geom_vline(xintercept = 0, linetype = 3)

## Weighted burden test

# build a burden score
wts = dbeta(mafs, 1, 25)
burden = G * wts
burden[1:5, 1:3]

# statistics on the burden score
score = rowSums(burden)
table(score > 1)

# burden test
fit = lm(phen$Pheno ~ score)
summary(fit)

## SKAT
null_skat = SKAT_Null_Model(phen$Pheno ~ 1 , out_type = "C")
res_skat = SKAT(G, null_skat)
str(res_skat)

# rho = 1 --> SKAT-O = burden
res_skat_rho1 = SKAT(G, null_skat, r.corr = 1)
res_skat_rho1$p.value

# rho = 0 --> SKAT-O = SKAT or variance component test
res_skat_rho1 = SKAT(G, null_skat, r.corr = 0)
res_skat_rho1$p.value

# optimal rho
res_skat_o = SKAT(G, null_skat, method = "optimal.adj")
res_skat_o$p.value

## ACAT
wts_beta = dbeta(mafs, 1, 25)
wts_acat <- wts_beta * wts_beta * mafs * (1 - mafs)

ACAT(assoc_plink$P, weights = wts_acat)

# ACAT vs. MinP
min(assoc_plink$P) * nrow(assoc_plink)


