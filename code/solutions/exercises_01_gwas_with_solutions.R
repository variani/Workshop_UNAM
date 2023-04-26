# Source: https://github.com/joellembatchou/SISG2022_Association_Mapping/blob/master/code/Session01_exercises.R
library(tidyverse)
library(data.table)
library(glue)

library(qqman)

## Define paths to data files/tools
# - helpful R commands: getwd(), list.files(), file.exists()
DIR_MAIN = '~/git/variani/UNAM2023_Association_Mapping/'
DIR_DATA = glue('{DIR_MAIN}/data/')
D1_PHENO = glue('{DIR_DATA}/sim_rels_pheno.txt')
D1_BFILE = glue('{DIR_DATA}/sim_rels_geno')
D2_PHENO = glue('{DIR_DATA}/bpdata.csv')

DIR_TOOLS = glue('{DIR_MAIN}/tools/')
PLINK = glue('{DIR_TOOLS}/plink2')
REGENIE = '/usr/local/Caskroom/miniconda/base/envs/regenie/bin/regenie'

# local functions
fun_lambda = function(pvals) median(qchisq(pvals, df = 1, lower.tail = FALSE)) / qchisq(0.5, df = 1, lower.tail = FALSE)

## Q1: Explore the BP dataset
# - Pick a SNP (snp1, snp2, ..., snp10) and make a boxplot SBP vs three SNP genotypes
#   -- What allele is minor?
# - Convert the selected SNP to a factor variable for further association analyses
bp = fread(D2_PHENO)

# plot snp1 vs. bp
ggplot(bp, aes(snp1, sbp)) + geom_boxplot()

# format snp3
summary(bp$snp3)

bp = within(bp, {
    snp3 <- factor(snp3, levels = c('CC', 'TC', 'TT'))
    g3 = as.numeric(snp3)
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
cmd = glue(" {PLINK} --bfile {D1_BFILE} ",
    " --pheno {D1_PHENO} --pheno-name Pheno ",
    " --maf 0.01 --geno 0.1 --hwe 1e-10 ",
    " --glm allow-no-covars --out gwas_plink ")
system(cmd)

# Q5: Explore association results from Plink
# - Make a QQ plot & compute the inflation factor
assoc_plink = fread('gwas_plink.Pheno.glm.linear') %>% as_tibble
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
cmd = glue(" {REGENIE} --bed {D1_BFILE} ",
    " --phenoFile {D1_PHENO} --phenoCol Pheno ",
    " --qt ",
    " --pred gwas_regenie_pred.list ",
    " --step 2 --bsize 400 ",
    " --out gwas_regenie ")
system(cmd)

# Q8: Explore association results from Regenie
# - Is the inflation under control? 
# - Is there any significan association?
assoc_regenie = fread('gwas_regenie_Pheno.regenie') %>% as_tibble
str(assoc_regenie)

assoc_regenie = within(assoc_regenie, P <- 10^(-LOG10P))
range(assoc_regenie$P)

# qq plot
qq(assoc_regenie$P)

# inflation factor lambda
fun_lambda(assoc_regenie$P)
