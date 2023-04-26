library(tidyverse)
library(data.table)
library(glue)

library(qqman)

## Paths to data files/tools
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

## Part 1.1: BP dataset
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

# additive association model
fit_add = lm(sbp ~ g3, bp)
summary(fit_add)

# adjust by covariates like sex
fit_add_cov1 = lm(sbp ~ sex + g3, bp)
summary(fit_add_cov1)

fit_add_cov2 = lm(sbp ~ sex + bmi + g3, bp)
summary(fit_add_cov2)

## Part 1.2: GWAS on sim_rels dataset
cmd = glue(" {PLINK} --bfile {D1_BFILE} ",
    " --pheno {D1_PHENO} --pheno-name Pheno ",
    " --maf 0.01 --geno 0.1 --hwe 1e-10 ",
    " --glm allow-no-covars --out gwas_plink ")
system(cmd)

# read association results from Plink
assoc_plink = fread('gwas_plink.Pheno.glm.linear') %>% as_tibble
str(assoc_plink)

# qq plot 
qq(assoc_plink$P)

# inflation factor lambda
fun_lambda(assoc_plink$P)

## run Step 1 Regenie
# --maf 0.01 --geno 0.1 --hwe 1e-10
# --extract <plink_QC_pass_snplist>

cmd = glue(" {REGENIE} --bed {D1_BFILE} ",
    " --phenoFile {D1_PHENO} --phenoCol Pheno ",
    " --qt ",
    " --step 1 --loocv --bsize 1000 ", 
    " --out gwas_regenie ")
system(cmd)

## run Step 2 Regenie
cmd = glue(" {REGENIE} --bed {D1_BFILE} ",
    " --phenoFile {D1_PHENO} --phenoCol Pheno ",
    " --qt ",
    " --pred gwas_regenie_pred.list ",
    " --step 2 --bsize 400 ",
    " --out gwas_regenie ")
system(cmd)

# read association results from Regenie
assoc_regenie = fread('gwas_regenie_Pheno.regenie') %>% as_tibble
str(assoc_regenie)

assoc_regenie = within(assoc_regenie, P <- 10^(-LOG10P))
range(assoc_regenie$P)

# qq plot
qq(assoc_regenie$P)

# inflation factor lambda
fun_lambda(assoc_regenie$P)
