# Source: https://github.com/joellembatchou/SISG2022_Association_Mapping/blob/master/code/Session02_exercises.R

library(tidyverse)
library(data.table)
library(glue)

library(BEDMatrix)
library(bigsnpr)
library(hexbin)

## Define paths to data files/tools
# - helpful R commands: getwd(), list.files(), file.exists()
DIR_MAIN = '~/git/variani/UNAM2023_Association_Mapping/'

DIR_DATA = glue('{DIR_MAIN}/data/')
D1_BFILE = glue('{DIR_DATA}/YRI_CEU_ASW_MEX_NAM')
D1_LABS = glue('{DIR_DATA}/Population_Sample_Info.txt')

DIR_TOOLS = glue('{DIR_MAIN}/tools/')
PLINK = glue('{DIR_TOOLS}/plink2')

# Q1. Explore the dataset from the 1,000 Genomes Project + (Ingideneous people in the Americas; NAM) the Human Genome Diversity Project
# https://www.coriell.org/1/NHGRI/Collections/1000-Genomes-Project-Collection/1000-Genomes-Project
# - How many samples/SNPs?
# - How many population groups?
bmat = BEDMatrix(D1_BFILE, simple_names = TRUE)

dim(bmat)

rownames(bmat) %>% head

# sample population labels
labs = fread(D1_LABS) %>% as_tibble

str(labs)

## Q2. Perform PCA by Plink using all SNPs; compute 10 PCs
# - Make PC1 vs PC2 scatterplot. Interpret PC1/PC2: what genetic ancestries they represent?
# - Check the output file with eigen values. Approximate the proportion of variance explained by each PC.
#   -- How many PCs should we keep so they underline the population structure?
#   -- What structure in data later PCs 9 and 10 show us?
cmd = glue(" {PLINK} --bfile {D1_BFILE} ",
    " --pca 10 --out pca_plink ")
system(cmd)

# scatterplot of PC1 vs PC2
pcs = fread('pca_plink.eigenvec') 
pcs = full_join(pcs, labs, by = 'IID')

ggplot(pcs, aes(PC1, PC2)) + geom_point(aes(color = Population))

ggplot(pcs, aes(PC3, PC4)) + geom_point(aes(color = Population))

ggplot(pcs, aes(PC9, PC10)) + geom_point(aes(color = Population))

# variance captured by PCs
evals = read_lines('pca_plink.eigenval') %>% as.numeric
prop = evals / sum(evals)

## Q3. Repeat the PCA using bigsnpr computing 10 PCs. 
# - Specify SNP pruning with R2 threshold and minimum MAC. Tip: check the defaults of bed_autoSVD function
# - Check again PC1 vs PC2 plot. Is it similar to that from Plink?
# - Check the loading plots. What special patterns for the later PCs? Is that reflected in the PC9 vs PC10 plot?
bedfile = glue('{D1_BFILE}.bed')
bed = bed(bedfile)
pca = bed_autoSVD(bed, thr.r2 = 0.2, min.mac = 20, k = 10)

# plot PC1 vs PC2
plot(pca, type = "scores", scores = 1:2) + aes(color = labs$Population)

# check SNP loadgings
plot(pca, type = "loadings", loadings = 1:10, coeff = 0.4)



