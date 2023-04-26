library(tidyverse)
library(data.table)
library(glue)

library(BEDMatrix)
library(bigsnpr)

## Paths to data files/tools
DIR_MAIN = '~/git/variani/UNAM2023_Association_Mapping/'
DIR_DATA = glue('{DIR_MAIN}/data/')
D1_BFILE = glue('{DIR_DATA}/YRI_CEU_ASW_MEX_NAM')
D1_LABS = glue('{DIR_DATA}/Population_Sample_Info.txt')

DIR_TOOLS = glue('{DIR_MAIN}/tools/')
PLINK = glue('{DIR_TOOLS}/plink2')

## Explore the dataset
bmat = BEDMatrix(D1_BFILE, simple_names = TRUE)

dim(bmat)

rownames(bmat) %>% head

# sample labels
labs = fread(D1_LABS) %>% as_tibble

str(labs)

## PCA by Plink
cmd = glue(" {PLINK} --bfile {D1_BFILE} ",
    " --pca 10 --out pca_plink ")
system(cmd)

