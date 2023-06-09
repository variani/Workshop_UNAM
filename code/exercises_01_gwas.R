# Source: https://github.com/joellembatchou/SISG2022_Association_Mapping/blob/master/code/Session01_exercises.R
library(tidyverse)
library(data.table)
library(glue)

library(qqman)

## Define paths to data files/tools
# - helpful R commands: getwd(), list.files(), file.exists()
DIR_DATA = 'data/' # or simple '.'
D1_PHENO = glue('{DIR_DATA}/sim_rels_pheno.txt')
D1_BFILE = glue('{DIR_DATA}/sim_rels_geno')
D2_PHENO = glue('{DIR_DATA}/bpdata.csv')

PLINK = 'plink2'
REGENIE = 'regenie'

# local functions
fun_lambda = function(pvals) median(qchisq(pvals, df = 1, lower.tail = FALSE)) / qchisq(0.5, df = 1, lower.tail = FALSE)

## Q1: Explore the BP dataset
# - Pick a SNP (snp1, snp2, ..., snp10) and make a boxplot SBP vs three SNP genotypes
#   -- What allele is minor?
# - Convert the selected SNP to a factor variable for further association analyses
bp = fread(D2_PHENO)

# Example plot for snp1 
ggplot(bp, aes(snp1, sbp)) + geom_boxplot()

# format snp3
summary(bp$snp3)

bp = within(bp, {
    snp3 <- factor(snp3, levels = c('CC', 'TC', 'TT')) # ordered genotype groups
    g3 = as.numeric(snp3) # additive coding for SNP: 0, 1 and 2
})
summary(bp$g3)

# Q2: Run an additive association model
fit_add = lm(sbp ~ g3, bp)
summary(fit_add)

# Q3: Adjust by covariates like sex and check if that helps for association
fit_add_cov1 = lm(sbp ~ sex + g3, bp)
summary(fit_add_cov1)

fit_add_cov2 = lm(sbp ~ sex + bmi + g3, bp)
summary(fit_add_cov2)

# Q4: Explore Dominant, Recessive and the full 2-parameter model.
# - Do results differ for your selected SNP?

## Q5: Run GWAS using Plink on sim_rels dataset (related samples)
# - Specify Quality Control filters in Plink: min MAF, genotype missingness level, the HWE test thrshold.
# - Do we control for relatedness?
val_maf = NULL # put some value here
val_geno = NULL # put some value here
val_hwe = NULL # put some value here
cmd = glue(" {PLINK} --bfile {D1_BFILE} ",
    " --pheno {D1_PHENO} --pheno-name Pheno ",
    " --maf {val_maf} --geno {val_geno} --hwe {val_hwe} ",
    " --glm allow-no-covars --out gwas_plink ")
system(cmd)

# Q5: Explore association results from Plink
# - Make a QQ plot & compute the inflation factor
f_assoc_plink = NULL # put the path to output path here
assoc_plink = fread(f_assoc_plink) %>% as_tibble
str(assoc_plink)

# qq plot 
qq(assoc_plink$P)

# inflation factor lambda
fun_lambda(assoc_plink$P)

## Q6: Run Step 1 Regenie
# - Does Regenie controls for relatedness? 
# - What other benefits of Regenei Step 1 compared to the Plink analysis?
cmd = glue(" {REGENIE} --bed {D1_BFILE} ",
    " --phenoFile {D1_PHENO} --phenoCol Pheno ",
    " --qt ",
    " --step 1 --loocv --bsize 1000 ", 
    " --out gwas_regenie ")
system(cmd)

## Q7: Run Step 2 Regenie
f_step1 = NULL # put the path of the step1 file here
cmd = glue(" {REGENIE} --bed {D1_BFILE} ",
    " --phenoFile {D1_PHENO} --phenoCol Pheno ",
    " --qt ",
    " --pred {f_step1} ",
    " --step 2 --bsize 400 ",
    " --out gwas_regenie ")
system(cmd)

# Q8: Explore association results from Regenie
# - Is the inflation under control? 
# - Is there any significan association?
f_assoc_regenie = NULL # put the path to the output Regenie file
assoc_regenie = fread(f_assoc_regenie) %>% as_tibble
str(assoc_regenie)

# convert LOG10P values reported by Regenie to raw P values
assoc_regenie = within(assoc_regenie, P <- 10^(-LOG10P))
range(assoc_regenie$P)

# qq plot
qq(assoc_regenie$P)

# inflation factor lambda
fun_lambda(assoc_regenie$P)
