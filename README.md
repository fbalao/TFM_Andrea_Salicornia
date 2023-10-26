# Transcriptome alteration in **Salicornia ramossisima** after exposure to multiple combined stresses
This folder contains the R code employed in the final master's thesis of the student Andrea Martín Díaz. October call of the academic year 2022/23 MADOBIS

# Check strandedness of RNAseq

check_strandedness --gtf ../../salicornia_genomes/Salicornia_europaea_2023/SAEU.gtf  --transcripts ../../salicornia_genomes/Salicornia_europaea_2023/SAEU_CDS.fasta --reads_1 HH3_1CWS_1_c.fq.gz --reads_2 HH3_1CWS_2_c.fq.gz

# STAR mapping RNA reads

## set conda enviroment
conda activate  STAR_mapping

## Index genome
### set variables
ASM=/mnt/nfs_storage9tb/export/data/andrea_salicornia/salicornia_genomes/Salicornia_europaea_2023/SAEU_assembly_Full.fasta

### indexing
STAR --runMode genomeGenerate --genomeSAindexNbases 13 \
--genomeDir /mnt/nfs_storage9tb/export/data/andrea_salicornia/salicornia_genomes/Salicornia_europaea_2023/STAR_map/star_index \
--genomeFastaFiles $ASM \
--sjdbGTFfile /mnt/nfs_storage9tb/export/data/andrea_salicornia/salicornia_genomes/Salicornia_europaea_2023/SAEU.gtf \
--sjdbOverhang 150 \
--runThreadN 23

### STAR mapping
cd /mnt/nfs_storage9tb/export/data/andrea_salicornia/salicornia_analysis/salicornia_cleanreads/
FILES=$(find . -name "*1_c.fq.gz" | rev | cut -c11- | rev)
for i in $FILES
do
STAR --runThreadN 21 \
--genomeDir /mnt/nfs_storage9tb/export/data/andrea_salicornia/salicornia_genomes/Salicornia_europaea_2023/STAR_map/star_index \
--readFilesCommand zcat \
--readFilesIn $i\_1_c\.fq.gz $i\_2_c\.fq.gz \
--outFileNamePrefix $i\_STAR \
--outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts \
--outFilterScoreMinOverLread 0.1 \
--outFilterMatchNminOverLread 0.1 \
--outFilterMismatchNmax 5
done

# Make the count table for R 

paste *ReadsPerGene.out.tab | grep -v "N_" | awk '{printf "%s\t", $1}{for (i=2;i<=NF;i+=2) printf "%s\t", $i; printf "\n" }' > tmpfile
sed -e "1igene_name\t$(ls *ReadsPerGene.out.tab | tr '\n' '\t' | sed 's/ReadsPerGene.out.tab//g')" tmpfile | cut -f1-97 > raw_counts_salicornia_matrix.txt
sed -i s/_STAR//g raw_counts_salicornia_matrix.txt

# Import genes counts into R :

``` r
list_paths<-list.files("/home/fbalao/Nextcloud/TFMAndrea/salicornia_europea/count_tables", full.names = T )

list_names<-substr(list.files("/home/fbalao/Nextcloud/TFMAndrea/salicornia_europea/count_tables"), 1, 8)

df_list<-lapply(list_paths, function(x){read.table(x,skip=4, row.names=1)$V2})

names(df_list)<-list_names

counts<-do.call(cbind.data.frame, df_list)
row.names(counts)<-row.names(read.table(list_paths[1],skip=4, row.names=1))
dim(counts)
 write.table(counts,"raw_counts_salicornia_matrix.txt")
```

## Import counts :

``` r
sample<-read.table("C:/Users/Usuario/Desktop/TFM/salmon/salmon_results/salico_samples2.txt", sep=";", header=T, as.is=T)
counts <- read.table("raw_counts_salicornia_matrix.txt", header = T, row.names = 1)
```

# PCA:

## Install Packages:

``` r
install_github("mixOmicsTeam/mixOmics")
install.packages("ggplot2")
```

## Libraries:

``` r
library(devtools)
library(mixOmics)
library(ggplot2)
```

``` r
mypca.rna = mixOmics::pca(t(counts), ncomp = 10, center = TRUE, scale = FALSE)

plot(mypca.rna)
```

``` r
mypca.rna = pca(t(counts), ncomp = 4, center = TRUE , scale = FALSE)
pca <- pca(t(counts), ncomp = 2, center = TRUE , scale = FALSE)
data <- data.frame(pca$variates)
muestras=rownames(t(counts))
water= sample[,8]

ggplot(data, aes(x = X.PC1, y = X.PC2, color = water)) +
  geom_point() +
  geom_text(aes(label = muestras), size = 3) +  
  labs(x = "PC1 (36%)", y = "PC2 (19%)", title = "") +
  theme_bw() +
  guides(color = guide_legend(title = NULL, override.aes = list(size = 4))) 
```

# Differential expression with EdgeR

## Install Packages:

``` r
BiocManager::install("edgeR")
install.packages("xfun")
```

## Libraries:

``` r
library(xfun)
library(edgeR)
```

## Model:

``` r
y <- DGEList(counts)
keep <- rowSums(cpm(y) >= 1) >= 3 # Filter counts less than 3
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y) # # Normalize by TMM
targets<- sample[,-c(1,2,5)]
targets$salinity<-as.factor(targets$salinity)
Group <- factor(paste(targets$co2,targets$temp,targets$salinity, targets$bact, targets$water, sep="."))
design <- model.matrix(~0+Group)
colnames(design) <- levels(Group)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
```

## Making the contrasts:

``` r
my.contrasts <- makeContrasts(
  CO2 = ( H.M.1.C.WW+H.H.1.C.WW+H.M.3.C.WW+H.H.3.C.WW+H.M.1.I.WW+H.H.1.I.WW+H.M.3.I.WW+H.H.3.I.WW+H.M.1.C.WS+H.H.1.C.WS+H.M.3.C.WS+H.H.3.C.WS+H.M.1.I.WS+H.H.1.I.WS+H.M.3.I.WS+H.H.3.I.WS )/16-( A.M.1.C.WW+A.H.1.C.WW+A.M.3.C.WW+A.H.3.C.WW+A.M.1.I.WW+A.H.1.I.WW+A.M.3.I.WW+A.H.3.I.WW+A.M.1.C.WS+A.H.1.C.WS+A.M.3.C.WS+A.H.3.C.WS+A.M.1.I.WS+A.H.1.I.WS+A.M.3.I.WS+A.H.3.I.WS )/16,
  Temp=( A.H.1.C.WW+H.H.1.C.WW+A.H.3.C.WW+H.H.3.C.WW+A.H.1.I.WW+H.H.1.I.WW+A.H.3.I.WW+H.H.3.I.WW+A.H.1.C.WS+H.H.1.C.WS+A.H.3.C.WS+H.H.3.C.WS+A.H.1.I.WS+H.H.1.I.WS+A.H.3.I.WS+H.H.3.I.WS )/16-( A.M.1.C.WW+H.M.1.C.WW+A.M.3.C.WW+H.M.3.C.WW+A.M.1.I.WW+H.M.1.I.WW+A.M.3.I.WW+H.M.3.I.WW+A.M.1.C.WS+H.M.1.C.WS+A.M.3.C.WS+H.M.3.C.WS+A.M.1.I.WS+H.M.1.I.WS+A.M.3.I.WS+H.M.3.I.WS )/16,
  NaCl = ( A.M.3.C.WW+H.M.3.C.WW+A.H.3.C.WW+H.H.3.C.WW+A.M.3.I.WW+H.M.3.I.WW+A.H.3.I.WW+H.H.3.I.WW+A.M.3.C.WS+H.M.3.C.WS+A.H.3.C.WS+H.H.3.C.WS+A.M.3.I.WS+H.M.3.I.WS+A.H.3.I.WS+H.H.3.I.WS )/16-( A.M.1.C.WW+H.M.1.C.WW+A.H.1.C.WW+H.H.1.C.WW+A.M.1.I.WW+H.M.1.I.WW+A.H.1.I.WW+H.H.1.I.WW+A.M.1.C.WS+H.M.1.C.WS+A.H.1.C.WS+H.H.1.C.WS+A.M.1.I.WS+H.M.1.I.WS+A.H.1.I.WS+H.H.1.I.WS )/16,
  Bact= ( A.M.1.I.WW+H.M.1.I.WW+A.H.1.I.WW+H.H.1.I.WW+A.M.3.I.WW+H.M.3.I.WW+A.H.3.I.WW+H.H.3.I.WW+A.M.1.I.WS+H.M.1.I.WS+A.H.1.I.WS+H.H.1.I.WS+A.M.3.I.WS+H.M.3.I.WS+A.H.3.I.WS+H.H.3.I.WS )/16-( A.M.1.C.WW+H.M.1.C.WW+A.H.1.C.WW+H.H.1.C.WW+A.M.3.C.WW+H.M.3.C.WW+A.H.3.C.WW+H.H.3.C.WW+A.M.1.C.WS+H.M.1.C.WS+A.H.1.C.WS+H.H.1.C.WS+A.M.3.C.WS+H.M.3.C.WS+A.H.3.C.WS+H.H.3.C.WS )/16,
  Wat = ( A.M.1.C.WS+H.M.1.C.WS+A.H.1.C.WS+H.H.1.C.WS+A.M.3.C.WS+H.M.3.C.WS+A.H.3.C.WS+H.H.3.C.WS+A.M.1.I.WS+H.M.1.I.WS+A.H.1.I.WS+H.H.1.I.WS+A.M.3.I.WS+H.M.3.I.WS+A.H.3.I.WS+H.H.3.I.WS )/16-( A.M.1.C.WW+H.M.1.C.WW+A.H.1.C.WW+H.H.1.C.WW+A.M.3.C.WW+H.M.3.C.WW+A.H.3.C.WW+H.H.3.C.WW+A.M.1.I.WW+H.M.1.I.WW+A.H.1.I.WW+H.H.1.I.WW+A.M.3.I.WW+H.M.3.I.WW+A.H.3.I.WW+H.H.3.I.WW )/16,
  CO2_Temp = ((A.H.1.C.WS+
A.H.3.C.WS+A.H.1.I.WS+A.H.3.I.WS+A.H.1.C.WW+A.H.3.C.WW+A.H.1.I.WW+A.H.3.I.WW)/8-(H.H.1.C.WS+H.H.3.C.WS+
H.H.1.I.WS+H.H.3.I.WS+H.H.1.C.WW+H.H.3.C.WW+H.H.1.I.WW+H.H.3.I.WW)/8)-((A.M.1.C.WS+A.M.3.C.WS+A.M.1.I.WS+
A.M.3.I.WS+A.M.1.C.WW+A.M.3.C.WW+A.M.1.I.WW+A.M.3.I.WW)/8-(H.M.1.C.WS+H.M.3.C.WS+H.M.1.I.WS+H.M.3.I.WS+
H.M.1.C.WW+H.M.3.C.WW+H.M.1.I.WW+H.M.3.I.WW)/8),
 CO2_Bact = (( A.M.1.I.WW+A.H.1.I.WW+A.M.3.I.WW+A.H.3.I.WW+A.M.1.I.WS+A.H.1.I.WS+A.M.3.I.WS+A.H.3.I.WS)/8-( H.M.1.I.WW+H.H.1.I.WW+H.M.3.I.WW+H.H.3.I.WW+H.M.1.I.WS+H.H.1.I.WS+H.M.3.I.WS+H.H.3.I.WS )/8)-(( A.M.1.C.WW+A.H.1.C.WW+A.M.3.C.WW+A.H.3.C.WW+A.M.1.C.WS+A.H.1.C.WS+A.M.3.C.WS+A.H.3.C.WS )/8-( H.M.1.C.WW+H.H.1.C.WW+H.M.3.C.WW+H.H.3.C.WW+H.M.1.C.WS+H.H.1.C.WS+H.M.3.C.WS+H.H.3.C.WS )/8),
  CO2_NaCl= (( A.M.3.C.WW+A.H.3.C.WW+A.M.3.I.WW+A.H.3.I.WW+A.M.3.C.WS+A.H.3.C.WS+A.M.3.I.WS+A.H.3.I.WS )/8-( H.M.3.C.WW+H.H.3.C.WW+H.M.3.I.WW+H.H.3.I.WW+H.M.3.C.WS+H.H.3.C.WS+H.M.3.I.WS+H.H.3.I.WS )/8)-(( A.M.1.C.WW+A.H.1.C.WW+A.M.1.I.WW+A.H.1.I.WW+A.M.1.C.WS+A.H.1.C.WS+A.M.1.I.WS+A.H.1.I.WS )/8-( H.M.1.C.WW+H.H.1.C.WW+H.M.1.I.WW+H.H.1.I.WW+H.M.1.C.WS+H.H.1.C.WS+H.M.1.I.WS+H.H.1.I.WS )/8),
  CO2_Wat= (( A.M.1.C.WS+A.H.1.C.WS+A.M.3.C.WS+A.H.3.C.WS+A.M.1.I.WS+A.H.1.I.WS+A.M.3.I.WS+A.H.3.I.WS )/8-( H.M.1.C.WS+H.H.1.C.WS+H.M.3.C.WS+H.H.3.C.WS+H.M.1.I.WS+H.H.1.I.WS+H.M.3.I.WS+H.H.3.I.WS )/8)-(( A.M.1.C.WW+A.H.1.C.WW+A.M.3.C.WW+A.H.3.C.WW+A.M.1.I.WW+A.H.1.I.WW+A.M.3.I.WW+A.H.3.I.WW )/8-( H.M.1.C.WW+H.H.1.C.WW+H.M.3.C.WW+H.H.3.C.WW+H.M.1.I.WW+H.H.1.I.WW+H.M.3.I.WW+H.H.3.I.WW )/8),
  Temp_NaCl = (( A.M.3.C.WW+H.M.3.C.WW+A.M.3.I.WW+H.M.3.I.WW+A.M.3.C.WS+H.M.3.C.WS+A.M.3.I.WS+H.M.3.I.WS )/8-( A.H.3.C.WW+H.H.3.C.WW+A.H.3.I.WW+H.H.3.I.WW+A.H.3.C.WS+H.H.3.C.WS+A.H.3.I.WS+H.H.3.I.WS )/8)-(( A.M.1.C.WW+H.M.1.C.WW+A.M.1.I.WW+H.M.1.I.WW+A.M.1.C.WS+H.M.1.C.WS+A.M.1.I.WS+H.M.1.I.WS )/8-( A.H.1.C.WW+H.H.1.C.WW+A.H.1.I.WW+H.H.1.I.WW+A.H.1.C.WS+H.H.1.C.WS+A.H.1.I.WS+H.H.1.I.WS )/8),
  Temp_Bact = (( A.M.1.I.WW+H.M.1.I.WW+A.M.3.I.WW+H.M.3.I.WW+A.M.1.I.WS+H.M.1.I.WS+A.M.3.I.WS+H.M.3.I.WS )/8-( A.H.1.I.WW+H.H.1.I.WW+A.H.3.I.WW+H.H.3.I.WW+A.H.1.I.WS+H.H.1.I.WS+A.H.3.I.WS+H.H.3.I.WS )/8)-(( A.M.1.C.WW+H.M.1.C.WW+A.M.3.C.WW+H.M.3.C.WW+A.M.1.C.WS+H.M.1.C.WS+A.M.3.C.WS+H.M.3.C.WS )/8-( A.H.1.C.WW+H.H.1.C.WW+A.H.3.C.WW+H.H.3.C.WW+A.H.1.C.WS+H.H.1.C.WS+A.H.3.C.WS+H.H.3.C.WS )/8),
  Temp_Wat = (( A.M.1.C.WS+H.M.1.C.WS+A.M.3.C.WS+H.M.3.C.WS+A.M.1.I.WS+H.M.1.I.WS+A.M.3.I.WS+H.M.3.I.WS )/8-( A.H.1.C.WS+H.H.1.C.WS+A.H.3.C.WS+H.H.3.C.WS+A.H.1.I.WS+H.H.1.I.WS+A.H.3.I.WS+H.H.3.I.WS )/8)-(( A.M.1.C.WW+H.M.1.C.WW+A.M.3.C.WW+H.M.3.C.WW+A.M.1.I.WW+H.M.1.I.WW+A.M.3.I.WW+H.M.3.I.WW )/8-( A.H.1.C.WW+H.H.1.C.WW+A.H.3.C.WW+H.H.3.C.WW+A.H.1.I.WW+H.H.1.I.WW+A.H.3.I.WW+H.H.3.I.WW )/8),
  NaCl_Bact = (( A.M.1.I.WW+H.M.1.I.WW+A.H.1.I.WW+H.H.1.I.WW+A.M.1.I.WS+H.M.1.I.WS+A.H.1.I.WS+H.H.1.I.WS )/8-( A.M.3.I.WW+H.M.3.I.WW+A.H.3.I.WW+H.H.3.I.WW+A.M.3.I.WS+H.M.3.I.WS+A.H.3.I.WS+H.H.3.I.WS )/8)-(( A.M.1.C.WW+H.M.1.C.WW+A.H.1.C.WW+H.H.1.C.WW+A.M.1.C.WS+H.M.1.C.WS+A.H.1.C.WS+H.H.1.C.WS )/8-( A.M.3.C.WW+H.M.3.C.WW+A.H.3.C.WW+H.H.3.C.WW+A.M.3.C.WS+H.M.3.C.WS+A.H.3.C.WS+H.H.3.C.WS )/8),
  NaCl_Wat = (( A.M.1.C.WS+H.M.1.C.WS+A.H.1.C.WS+H.H.1.C.WS+A.M.1.I.WS+H.M.1.I.WS+A.H.1.I.WS+H.H.1.I.WS )/8-( A.M.3.C.WS+H.M.3.C.WS+A.H.3.C.WS+H.H.3.C.WS+A.M.3.I.WS+H.M.3.I.WS+A.H.3.I.WS+H.H.3.I.WS )/8)-(( A.M.1.C.WW+H.M.1.C.WW+A.H.1.C.WW+H.H.1.C.WW+A.M.1.I.WW+H.M.1.I.WW+A.H.1.I.WW+H.H.1.I.WW )/8-( A.M.3.C.WW+H.M.3.C.WW+A.H.3.C.WW+H.H.3.C.WW+A.M.3.I.WW+H.M.3.I.WW+A.H.3.I.WW+H.H.3.I.WW )/8),
Bact_Wat = (( A.M.1.C.WS+H.M.1.C.WS+A.H.1.C.WS+H.H.1.C.WS+A.M.3.C.WS+H.M.3.C.WS+A.H.3.C.WS+H.H.3.C.WS )/8-( A.M.1.I.WS+H.M.1.I.WS+A.H.1.I.WS+H.H.1.I.WS+A.M.3.I.WS+H.M.3.I.WS+A.H.3.I.WS+H.H.3.I.WS )/8)-(( A.M.1.C.WW+H.M.1.C.WW+A.H.1.C.WW+H.H.1.C.WW+A.M.3.C.WW+H.M.3.C.WW+A.H.3.C.WW+H.H.3.C.WW )/8-( A.M.1.I.WW+H.M.1.I.WW+A.H.1.I.WW+H.H.1.I.WW+A.M.3.I.WW+H.M.3.I.WW+A.H.3.I.WW+H.H.3.I.WW )/8),
levels=design)
```

## Average effect:

### Packages:

``` r
#install.packages("dplyr")
#install.packages("tidyr")
```

### Libraries:

``` r
library(dplyr)
library(tidyr)
library(magrittr)
```

### Function to get the topgenes:

``` r
topTags2<- function(glmfit,lfc=log2(1.5), fdr=0.05 ){
  z<-topTags(glmfit, n=dim(glmfit$table)[1])
  z2<-z$table
  z3<-z2[abs(z2$logFC)>lfc & z2$FDR<fdr,]
  z3
}
```

### Average effect of the CO2:

``` r
qlf_co2 <- glmQLFTest(fit, contrast=my.contrasts[,"CO2"])
summary(decideTests(qlf_co2, lfc=log2(1.5)))
```

#### Top genes:

``` r
toptagsCo2 = topTags2(qlf_co2)
dim(toptagsCo2)
```

#### Plot:

``` r
topco2 <- rownames(topTags(qlf_co2,sort.by = "logFC", p.value = 0.05)) # order by logFC and select the significant p.value
logCPMco2<-cpm(y)[topco2,]
```

``` r
logCPM<-as.data.frame(logCPMco2)
logCPM$gene = row.names(logCPM)
d = logCPM %>% gather(Sample, logCPM, -gene)
d$group<-rep(sample[,3], each=10)  # 10 is the  number of genes
geneI<-'Seu_jg2401' # Select one gene
geneOfInterest = d %>% filter(gene == geneI)
ggplot(geneOfInterest, aes(x=group, y=logCPM)) + geom_violin(aes(fill=group)) +
  theme_bw() +  labs(title=geneI) + xlab('Treatment') +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
```

### Average effect of the Temp:

``` r
qlf_Temp <- glmQLFTest(fit, contrast=my.contrasts[,"Temp"])
summary(decideTests(qlf_Temp, lfc=log2(1.5)))
```

#### Top genes:

``` r
toptagsTemp = topTags2(qlf_Temp)
dim(toptagsTemp)
```

#### Plot:

``` r
topTemp <- rownames(topTags(qlf_Temp,sort.by = "logFC", p.value = 0.05)) # order by logFC and select the significant p.value
logCPMTemp<-cpm(y)[topTemp,]
```

``` r
logCPM<-as.data.frame(logCPMTemp)
logCPM$gene = row.names(logCPM)
d = logCPM %>% gather(Sample, logCPM, -gene)
d$group<-rep(sample[,4], each=10)  # 10 is the  number of genes
geneI<-'Seu_jg2033' # Select one gene
geneOfInterest = d %>% filter(gene == geneI)
ggplot(geneOfInterest, aes(x=group, y=logCPM)) + geom_violin(aes(fill=group)) +
  theme_bw() +  labs(title=geneI) + xlab('Treatment') +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
```

### Average effect of the NaCl:

``` r
qlf_NaCl <- glmQLFTest(fit, contrast=my.contrasts[,"NaCl"])
summary(decideTests(qlf_NaCl, lfc=log2(1.5)))
```

#### Top genes:

``` r
toptagsNaCl = topTags2(qlf_NaCl)
dim(toptagsNaCl)
```

#### Plot:

``` r
topNaCl <- rownames(topTags(qlf_NaCl,sort.by = "logFC", p.value = 0.05)) # order by logFC and select the significant p.value
logCPMNaCl<-cpm(y)[topNaCl,]
```

``` r
logCPM<-as.data.frame(logCPMNaCl)
logCPM$gene = row.names(logCPM)
d = logCPM %>% gather(Sample, logCPM, -gene)
d$group<-as.factor(rep(sample[,6], each=10))  # 10 is the  number of genes
geneI<-'Seu_jg14388' # Select one gene
geneOfInterest = d %>% filter(gene == geneI)
ggplot(geneOfInterest, aes(x=group, y=logCPM)) + geom_violin(aes(fill=group)) +
  theme_bw() +  labs(title=geneI) + xlab('Treatment') +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
```

### Average effect of the Bact:

``` r
qlf_Bact <- glmQLFTest(fit, contrast=my.contrasts[,"Bact"])
summary(decideTests(qlf_Bact, lfc=log2(1.5)))
```

#### Top genes:

``` r
toptagsBact = topTags2(qlf_Bact)
dim(toptagsBact)
```

#### Plot:

``` r
topBact <- rownames(topTags(qlf_Bact,sort.by = "logFC", p.value = 0.05)) # order by logFC and select the significant p.value
logCPMBact<-cpm(y)[topBact,]
```

``` r
logCPM<-as.data.frame(logCPMBact)
logCPM$gene = row.names(logCPM)
d = logCPM %>% gather(Sample, logCPM, -gene)
d$group<-rep(sample[,7], each=10)  # 10 is the  number of genes
geneI<-'Seu_jg4696' # Select one gene
geneOfInterest = d %>% filter(gene == geneI)
ggplot(geneOfInterest, aes(x=group, y=logCPM)) + geom_violin(aes(fill=group)) +
  theme_bw() +  labs(title=geneI) + xlab('Treatment') +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
```

### Average effect of the Wat:

``` r
qlf_Wat <- glmQLFTest(fit, contrast=my.contrasts[,"Wat"])
summary(decideTests(qlf_Wat, lfc=log2(1.5)))
```

#### Top genes:

``` r
toptagsWat = topTags2(qlf_Wat)
dim(toptagsWat)
```

#### Plot:

``` r
topWat <- rownames(topTags(qlf_Wat,sort.by = "logFC", p.value = 0.05)) # order by logFC and select the significant p.value
logCPMWat<-cpm(y)[topWat,]
```

``` r
logCPM<-as.data.frame(logCPMWat)
logCPM$gene = row.names(logCPM)
d = logCPM %>% gather(Sample, logCPM, -gene)
d$group<-rep(sample[,8], each=10)  # 10 is the  number of genes
geneI<-'Seu_jg16393' # Select one gene
geneOfInterest = d %>% filter(gene == geneI)
ggplot(geneOfInterest, aes(x=group, y=logCPM)) + geom_violin(aes(fill=group)) +
  theme_bw() +  labs(title=geneI) + xlab('Treatment') +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
```

### Average effect of the CO2_Temp:

``` r
qlf_CO2_Temp <- glmQLFTest(fit, contrast=my.contrasts[,"CO2_Temp"])
summary(decideTests(qlf_CO2_Temp, lfc=log2(1.5)))
```

#### Top genes:

``` r
toptagsCO2_Temp = topTags2(qlf_CO2_Temp)
dim(toptagsCO2_Temp)
```

#### Plot:

``` r
topCO2_Temp <- rownames(topTags(qlf_CO2_Temp,sort.by = "logFC", p.value = 0.05)) # order by logFC and select the significant p.value
logCPMCO2_Temp<-cpm(y)[topCO2_Temp,]
```

``` r
logCPM<-as.data.frame(logCPMCO2_Temp)
logCPM$gene = row.names(logCPM)
d = logCPM %>% gather(Sample, logCPM, -gene)
d$group<-rep(interaction(sample[,3],sample[,4]), each=10)  # 10 is the  number of genes
geneI<-'Seu_jg16187' # Select one gene
geneOfInterest = d %>% filter(gene == geneI)
ggplot(geneOfInterest, aes(x=group, y=logCPM)) + geom_violin(aes(fill=group)) +
  theme_bw() +  labs(title=geneI) + xlab('Treatment') +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
```

### Average effect of the CO2_NaCl:

``` r
qlf_CO2_NaCl <- glmQLFTest(fit, contrast=my.contrasts[,"CO2_NaCl"])
summary(decideTests(qlf_CO2_NaCl, lfc=log2(1.5)))
```

#### Top genes:

``` r
toptagsCO2_NaCl = topTags2(qlf_CO2_NaCl)
dim(toptagsCO2_NaCl)
```

#### Plot:

``` r
topCO2_NaCl <- rownames(topTags(qlf_CO2_NaCl,sort.by = "logFC", p.value = 0.05)) # order by logFC and select the significant p.value
logCPMCO2_NaCl<-cpm(y)[topCO2_NaCl,]
```

``` r
logCPM<-as.data.frame(logCPMCO2_NaCl)
logCPM$gene = row.names(logCPM)
d = logCPM %>% gather(Sample, logCPM, -gene)
d$group<- interaction(sample[,4],sample[,7])  # 10 is the  number of genes
geneI<-'Seu_jg9924' # Select one gene
geneOfInterest = d %>% filter(gene == geneI)
ggplot(geneOfInterest, aes(x=group, y=logCPM)) + geom_violin(aes(fill=group)) +
  theme_bw() +  labs(title=geneI) + xlab('Treatment') +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
```

### Average effect of the CO2_Bact:

``` r
qlf_CO2_Bact <- glmQLFTest(fit, contrast=my.contrasts[,"CO2_Bact"])
summary(decideTests(qlf_CO2_Bact, lfc=log2(1.5)))
```

#### Top genes:

``` r
toptagsCO2_Bact = topTags2(qlf_CO2_Bact)
dim(toptagsCO2_Bact)
```

#### Plot:

``` r
topCO2_Bact <- rownames(topTags(qlf_CO2_Bact,sort.by = "logFC", p.value = 0.05)) # order by logFC and select the significant p.value
logCPMCO2_Bact<-cpm(y)[topCO2_Bact,]
```

``` r
logCPM<-as.data.frame(logCPMCO2_Bact)
logCPM$gene = row.names(logCPM)
d = logCPM %>% gather(Sample, logCPM, -gene)
d$group<-interaction(sample[,4],sample[,7], each=10)  # 10 is the  number of genes
geneI<-'Seu_jg26252' # Select one gene
geneOfInterest = d %>% filter(gene == geneI)
ggplot(geneOfInterest, aes(x=group, y=logCPM)) + geom_violin(aes(fill=group)) +
  theme_bw() +  labs(title=geneI) + xlab('Treatment') +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
```

### Average effect of the CO2_Wat:

``` r
qlf_CO2_Wat <- glmQLFTest(fit, contrast=my.contrasts[,"CO2_Wat"])
summary(decideTests(qlf_CO2_Wat, lfc=log2(1.5)))
```

#### Top genes:

``` r
toptagsCO2_Wat = topTags2(qlf_CO2_Wat)
dim(toptagsCO2_Wat)
```

#### Plot:

``` r
topCO2_Wat <- rownames(topTags(qlf_CO2_Wat,sort.by = "logFC", p.value = 0.05)) # order by logFC and select the significant p.value
logCPMCO2_Wat<-cpm(y)[topCO2_Wat,]
```

``` r
logCPM<-as.data.frame(logCPMCO2_Wat)
logCPM$gene = row.names(logCPM)
d = logCPM %>% gather(Sample, logCPM, -gene)
d$group<-rep(interaction(sample[,3],sample[,8]), each=10)  # 10 is the  number of genes
geneI<-'Seu_jg10859' # Select one gene
geneOfInterest = d %>% filter(gene == geneI)
ggplot(geneOfInterest, aes(x=group, y=logCPM)) + geom_violin(aes(fill=group)) +
  theme_bw() +  labs(title=geneI) + xlab('Treatment') +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
```

### Average effect of the Temp_NaCl:

``` r
qlf_Temp_NaCl <- glmQLFTest(fit, contrast=my.contrasts[,"Temp_NaCl"])
summary(decideTests(qlf_Temp_NaCl, lfc=log2(1.5)))
```

#### Top genes:

``` r
toptagsTemp_NaCl = topTags2(qlf_Temp_NaCl)
dim(toptagsTemp_NaCl)
```

#### Plot:

``` r
topTemp_NaCl <- rownames(topTags(qlf_Temp_NaCl,sort.by = "logFC", p.value = 0.05)) # order by logFC and select the significant p.value
logCPMTemp_NaCl<-cpm(y)[topTemp_NaCl,]
```

``` r
# This would be the code but since there are no genes it cannot be represented:
logCPM<-as.data.frame(logCPMTemp_NaCl)
logCPM$gene = row.names(logCPM)
d = logCPM %>% gather(Sample, logCPM, -gene)
d$group<-interaction(sample[,4],sample[,6])  # 10 is the  number of genes
geneI<-'Seu_jg8665' # Select one gene
geneOfInterest = d %>% filter(gene == geneI)
ggplot(geneOfInterest, aes(x=group, y=logCPM)) + geom_violin(aes(fill=group)) +
  theme_bw() +  labs(title=geneI) + xlab('Treatment') +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
```

### Average effect of the Temp_Bact:

``` r
qlf_Temp_Bact <- glmQLFTest(fit, contrast=my.contrasts[,"Temp_Bact"])
summary(decideTests(qlf_Temp_Bact, lfc=log2(1.5)))
```

#### Top genes:

``` r
toptagsTemp_Bact = topTags2(qlf_Temp_Bact)
dim(toptagsTemp_Bact)
```

#### Plot:

``` r
topTemp_Bact <- rownames(topTags(qlf_Temp_Bact,sort.by = "logFC", p.value = 0.05)) # order by logFC and select the significant p.value
logCPMTemp_Bact<-cpm(y)[topTemp_Bact,]
```

``` r
logCPM<-as.data.frame(logCPMTemp_Bact)
logCPM$gene = row.names(logCPM)
d = logCPM %>% gather(Sample, logCPM, -gene)
d$group<-rep(interaction(sample[,4],sample[,7]), each=10)  # 10 is the  number of genes
geneI<-'Seu_jg23107' # Select one gene
geneOfInterest = d %>% filter(gene == geneI)
ggplot(geneOfInterest, aes(x=group, y=logCPM)) + geom_violin(aes(fill=group)) +
  theme_bw() +  labs(title=geneI) + xlab('Treatment') +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
```

### Average effect of the Temp_Wat:

``` r
qlf_Temp_Wat <- glmQLFTest(fit, contrast=my.contrasts[,"Temp_Wat"])
summary(decideTests(qlf_Temp_Wat, lfc=log2(1.5)))
```

#### Top genes:

``` r
toptagsTemp_Wat = topTags2(qlf_Temp_Wat)
dim(toptagsTemp_Wat)
```

#### Plot:

``` r
topTemp_Wat <- rownames(topTags(qlf_Temp_Wat,sort.by = "logFC", p.value = 0.05)) # order by logFC and select the significant p.value
logCPMTemp_Wat<-cpm(y)[topTemp_Wat,]
```

``` r
logCPM<-as.data.frame(logCPMTemp_Wat)
logCPM$gene = row.names(logCPM)
d = logCPM %>% gather(Sample, logCPM, -gene)
d$group<-rep(interaction(sample[,4],sample[,8]), each=10)  # 10 is the  number of genes
geneI<-'Seu_jg19766' # Select one gene
geneOfInterest = d %>% filter(gene == geneI)
ggplot(geneOfInterest, aes(x=group, y=logCPM)) + geom_violin(aes(fill=group)) +
  theme_bw() +  labs(title=geneI) + xlab('Treatment') +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
```

### Average effect of the NaCl_Bact:

``` r
qlf_NaCl_Bact <- glmQLFTest(fit, contrast=my.contrasts[,"NaCl_Bact"])
summary(decideTests(qlf_NaCl_Bact, lfc=log2(1.5)))
```

#### Top genes:

``` r
toptagsNaCl_Bact = topTags2(qlf_NaCl_Bact)
dim(toptagsNaCl_Bact)
```

#### Plot:

``` r
topNaCl_Bact <- rownames(topTags(qlf_NaCl_Bact,sort.by = "logFC", p.value = 0.05)) # order by logFC and select the significant p.value
logCPMNaCl_Bact<-cpm(y)[topNaCl_Bact,]
```

``` r
# This would be the code but since there are no genes it cannot be represented:
logCPM<-as.data.frame(logCPMNaCl_Bact)
logCPM$gene = row.names(logCPM)
d = logCPM %>% gather(Sample, logCPM, -gene)
d$group<-rep(interaction(sample[,4],sample[,7]), each=10)  # 10 is the  number of genes
geneI<-'Seu_jg18063' # Select one gene
geneOfInterest = d %>% filter(gene == geneI)
ggplot(geneOfInterest, aes(x=group, y=logCPM)) + geom_violin(aes(fill=group)) +
  theme_bw() +  labs(title=geneI) + xlab('Treatment') +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
```

### Average effect of the NaCl_Wat:

``` r
qlf_NaCl_Wat <- glmQLFTest(fit, contrast=my.contrasts[,"NaCl_Wat"])
summary(decideTests(qlf_NaCl_Wat, lfc=log2(1.5)))
```

#### Top genes:

``` r
toptagsNaCl_Wat = topTags2(qlf_NaCl_Wat)
dim(toptagsNaCl_Wat)
```

#### Plot:

``` r
topNaCl_Wat <- rownames(topTags(qlf_NaCl_Wat,sort.by = "logFC", p.value = 0.05)) # order by logFC and select the significant p.value
logCPMNaCl_Wat<-cpm(y)[topNaCl_Wat,]
```

``` r
logCPM<-as.data.frame(logCPMNaCl_Wat)
logCPM$gene = row.names(logCPM)
d = logCPM %>% gather(Sample, logCPM, -gene)
d$group<-rep(interaction(sample[,4],sample[,7]), each=10)  # 10 is the  number of genes
geneI<-'Seu_jg27310' # Select one gene
geneOfInterest = d %>% filter(gene == geneI)
ggplot(geneOfInterest, aes(x=group, y=logCPM)) + geom_violin(aes(fill=group)) +
  theme_bw() +  labs(title=geneI) + xlab('Treatment') +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
```

### Average effect of the Bact_Wat:

``` r
qlf_Bact_Wat <- glmQLFTest(fit, contrast=my.contrasts[,"Bact_Wat"])
summary(decideTests(qlf_Bact_Wat, lfc=log2(1.5)))
```

#### Top genes:

``` r
toptagsBact_Wat = topTags2(qlf_Bact_Wat)
dim(toptagsBact_Wat)
```

#### Plot:

``` r
topBact_Wat <- rownames(topTags(qlf_Bact_Wat,sort.by = "logFC", p.value = 0.05)) # order by logFC and select the significant p.value
logCPMBact_Wat<-cpm(y)[topBact_Wat,]
```

``` r
# This would be the code but since there are no genes it cannot be represented:
logCPM<-as.data.frame(logCPMBact_Wat)
logCPM$gene = row.names(logCPM)
d = logCPM %>% gather(Sample, logCPM, -gene)
d$group<-interaction(sample[,4],sample[,7]) # 10 is the  number of genes
geneI<-'Seu_jg14599' # Select one gene
geneOfInterest = d %>% filter(gene == geneI)
ggplot(geneOfInterest, aes(x=group, y=logCPM)) + geom_violin(aes(fill=group)) +
  theme_bw() +  labs(title=geneI) + xlab('Treatment') +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
```

## Volcano plots:

### CO2

``` r
CO2_volcano = as.data.frame(topTags(qlf_co2, n= 18246))

ggplot(CO2_volcano, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = ifelse(logFC > log2(1.5), ifelse(FDR < 0.05  , "red", "grey"), ifelse(logFC < -log2(1.5), ifelse(FDR < 0.05  , "blue", "grey"), "grey")))) +
  xlab("logFC") +
  ylab("-log10(FDR)") +
  scale_color_identity() +
  theme_minimal()
```

### Temp

``` r
Temp_volcano = as.data.frame(topTags(qlf_Temp, n= 18246))

ggplot(Temp_volcano, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = ifelse(logFC > log2(1.5), ifelse(FDR < 0.05  , "red", "grey"), ifelse(logFC < -log2(1.5), ifelse(FDR < 0.05  , "blue", "grey"), "grey")))) +
  xlab("logFC") +
  ylab("-log10(FDR)") +
  scale_color_identity() +
  theme_minimal()
```

### NaCl

``` r
NaCl_volcano = as.data.frame(topTags(qlf_NaCl, n= 18246))

ggplot(NaCl_volcano, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = ifelse(logFC > log2(1.5), ifelse(FDR < 0.05  , "red", "grey"), ifelse(logFC < -log2(1.5), ifelse(FDR < 0.05  , "blue", "grey"), "grey")))) +
  xlab("logFC") +
  ylab("-log10(FDR)") +
  scale_color_identity() +
  theme_minimal()
```

### Bact

``` r
bact_volcano = as.data.frame(topTags(qlf_Bact, n= 18246))

ggplot(bact_volcano, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = ifelse(logFC > log2(1.5), ifelse(FDR < 0.05  , "red", "grey"), ifelse(logFC < -log2(1.5), ifelse(FDR < 0.05  , "blue", "grey"), "grey")))) +
  xlab("logFC") +
  ylab("-log10(FDR)") +
  scale_color_identity() +
  theme_minimal()
```

### Wat

``` r
Wat_volcano = as.data.frame(topTags(qlf_Wat, n= 18246))

ggplot(Wat_volcano, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = ifelse(logFC > log2(1.5), ifelse(FDR < 0.05  , "red", "grey"), ifelse(logFC < -log2(1.5), ifelse(FDR < 0.05  , "blue", "grey"), "grey")))) +
  xlab("logFC") +
  ylab("-log10(FDR)") +
  scale_color_identity() +
  theme_minimal()
```

# bar Char:

``` r
tratamientos <- c("CO2","Temp","NaCl","Bact","Wat","CO2-Temp","CO2-NaCl","CO2-Bact","CO2-Wat","Temp-NaCl","Temp-Bact","Temp-Wat","NaCl-Bact","NaCl-Wat","Bact-Wat")
upregulados <- c(557,709,312,0,1695,1180,0,4,453,0,32,231,0,8,0)
downregulated <- c(502,1022,188,4,3405,1642,3,2,424,0,9,287,0,12,0)
datos <- data.frame(tratamientos, upregulados, downregulated)

ggplot(datos, aes(x = tratamientos, y = upregulados, fill = "Upregulados")) +
  geom_bar(stat = "identity") +
  geom_bar(aes(y = -downregulated, fill = "Downregulated"), stat = "identity", position = "identity") +
  geom_text(aes(label = upregulados), vjust = -0.5, color = "black", position = position_stack(vjust = 0.5), size = 3) +
  geom_text(aes(label = downregulated), vjust = 1.5, color = "black", position = position_stack(vjust = -0.5), size = 3) +
  labs(x = "Tratamientos", y = "Number of genes", fill = "") +
  scale_fill_manual(values = c("Upregulados" = "cornflowerblue", "Downregulated" = "cadetblue2")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_y_continuous(labels = scales::comma_format(scale = 0.001))
```

# HeatMap:

## Pack:

``` r
install.packages("pheatmap")
```

## Library:

``` r
library(pheatmap)
```

## CO2 :

``` r
topco2 <- rownames(toptagsCo2)
logCPMco2<-cpm(y)[topco2,]
pheatmap(logCPMco2, scale = "row")
```

## Temp:

``` r
topTemp <- rownames(toptagsTemp)
logCPMTemp<-cpm(y)[topTemp,]
pheatmap(logCPMTemp, scale = "row")
```

## NaCl :

``` r
topNaCl <- rownames(toptagsNaCl)
logCPMNaCl<-cpm(y)[topNaCl,]
pheatmap(logCPMNaCl, scale = "row")
```

## Bact:

``` r
topBact <- rownames(toptagsBact)
logCPMBact<-cpm(y)[topBact,]
pheatmap(logCPMBact, scale = "row")
```

## Wat:

``` r
topWat <- rownames(toptagsWat)
logCPMWat<-cpm(y)[topWat,]
pheatmap(logCPMWat, scale = "row")
```

## CO2-Temp:

``` r
topCO2_Temp <- rownames(toptagsCO2_Temp)
logCPMCO2_Temp<-cpm(y)[topCO2_Temp,]
pheatmap(logCPMCO2_Temp, scale = "row")
```

## CO2-NaCl :

``` r
topCO2_NaCl <- rownames(toptagsCO2_NaCl)
logCPMCO2_NaCl<-cpm(y)[topCO2_NaCl,]
pheatmap(logCPMCO2_NaCl, scale = "row")
```

## CO2-Bact :

``` r
topCO2_Bact <- rownames(toptagsCO2_Bact)
logCPMCO2_Bact<-cpm(y)[topCO2_Bact,]
pheatmap(logCPMCO2_Bact, scale = "row")
```

## CO2-Wat :

``` r
topco2_wat <- rownames(toptagsCO2_Wat)
logCPMco2_wat<-cpm(y)[topco2_wat,]
pheatmap(logCPMco2_wat, scale = "row")
```

## Temp-NaCl :

``` r
# This would be the code but since there are no genes it cannot be represented:
topTemp_NaCl <- rownames(toptagsTemp_NaCl)
logCPMTemp_NaCl<-cpm(y)[topTemp_NaCl,]
pheatmap(logCPMTemp_NaCl, scale = "row")
```

## Temp-Bact :

``` r
# This would be the code but since there are no genes it cannot be represented :
topTemp_Bact <- rownames(topTemp_Bact)
logCPMTemp_Bact<-cpm(y)[topTemp_Bact,]
pheatmap(logCPMTemp_Bact, scale = "row")
```

## Temp-Wat :

``` r
# This would be the code but since there are no genes it cannot be represented:
topTemp_Wat <- rownames(topTemp_Wat)
logCPMTemp_Wat<-cpm(y)[topTemp_Wat,]
pheatmap(logCPMTemp_Wat, scale = "row")
```

## NaCl-Bact :

``` r
# This would be the code but since there are no genes it cannot be represented:
topNaCl_Bact <- rownames(toptagsNaCl_Bact)
logCPMNaCl_Bact<-cpm(y)[topNaCl_Bact,]
pheatmap(logCPMNaCl_Bact, scale = "row")
```

## NaCl-Wat :

``` r
topNaCl_Wat <- rownames(toptagsNaCl_Wat)
logCPMNaCl_Wat<-cpm(y)[topNaCl_Wat,]
pheatmap(logCPMNaCl_Wat, scale = "row")
```

## Bact-Wat :

``` r
# This would be the code but since there are no genes it cannot be represented:
topBact_Wat <- rownames(toptagsBact_Wat)
logCPMBact_Wat<-cpm(y)[topBact_Wat,]
pheatmap(logCPMBact_Wat ,scale= "row")
```

# Venn :

## Pack :

``` r
install.packages("ggVennDiagram")
```

## Libraries :

``` r
library('ggVennDiagram')
```

## All contrast :

``` r
condicion<- list(A= rownames(toptagsCo2),
B=rownames(toptagsTemp),
C=rownames(toptagsBact),
D=rownames(toptagsNaCl),
E= rownames(toptagsWat))

ggVennDiagram(condicion, label= "count",label_alpha = 0, category.names = c("CO2","Temp","Bact","NaCl","Wat"),set_size = 2)+
scale_fill_gradient(low= "lightpink", high= "lightblue")+ ggtitle("All conditions")
```

## CO2-Temp :

``` r
condicion_CO2_Temp<- list(A= rownames(toptagsCo2),
                 B=rownames(toptagsTemp),
                 C=rownames(toptagsCO2_Temp))
ggVennDiagram(condicion_CO2_Temp, label= "count",label_alpha = 0, category.names = c("CO2","Temp","CO2_Temp"), set_size = 2)+
  scale_fill_gradient(low= "lightpink", high= "lightblue",)+ ggtitle("CO2-Temperature") 
```

## CO2-NaCl:

``` r
condicion_CO2_NaCl<- list(A= rownames(toptagsCo2),
                 B=rownames(toptagsNaCl),
                 C=rownames(toptagsCO2_NaCl))
ggVennDiagram(condicion_CO2_NaCl,label= "count",label_alpha = 0, category.names = c("CO2","NaCl","CO2_NaCl"), set_size=2 )+
  scale_fill_gradient(low= "lightpink", high= "lightblue",)+ ggtitle("CO2-NaCl") 
```

## CO2-Bact :

``` r
 condicion_CO2_Bact<- list(A= rownames(toptagsCo2),
                 B=rownames(toptagsBact),
                 C=rownames(toptagsCO2_Bact))
ggVennDiagram(condicion_CO2_Bact,label= "count",label_alpha = 0, category.names = c("CO2","Bact","CO2_Bact"), set_size = 2)+
  scale_fill_gradient(low= "lightpink", high= "lightblue",) + ggtitle("CO2-Bacterium")
```

## CO2-Wat:

``` r
condicion_CO2_Wat<- list(A= rownames(toptagsCo2),
                 B=rownames(toptagsWat),
                 C=rownames(toptagsCO2_Wat))
ggVennDiagram(condicion_CO2_Wat,label= "count",label_alpha = 0, category.names = c("CO2","Wat","CO2_Wat"), set_size = 2)+
  scale_fill_gradient(low= "lightpink", high= "lightblue",) + ggtitle("CO2-Water")
```

## Temp-NaCl:

``` r
condicion_Temp_NaCl<- list(A= rownames(toptagsTemp),
                 B=rownames(toptagsNaCl),
                 C=rownames(toptagsTemp_NaCl))
ggVennDiagram(condicion_Temp_NaCl,label= "count",label_alpha = 0, category.names = c("Temp","NaCl","T_NaCl"),set_size = 2)+
  scale_fill_gradient(low= "lightpink", high= "lightblue",) + ggtitle("Temperature-NaCl")
```

## Temp-Bact:

``` r
condicion_Temp_Bact<- list(A= rownames(toptagsTemp),
                 B=rownames(toptagsBact),
                 C=rownames(toptagsTemp_Bact))
ggVennDiagram(condicion_Temp_Bact,label= "count",label_alpha = 0, category.names = c("Temp","Bact","Temp_Bact"), set_size = 2)+
  scale_fill_gradient(low= "lightpink", high= "lightblue",)+ ggtitle("Temperature-Bacterium")
```

## Temp-Wat:

``` r
condicion_Temp_Wat<- list(A= rownames(toptagsTemp),
                 B=rownames(toptagsWat),
                 C=rownames(toptagsTemp_Wat))
ggVennDiagram(condicion_Temp_Wat,label= "count",label_alpha = 0, category.names = c("Temp","Wat","Temp_Wat"), set_size = 2)+
  scale_fill_gradient(low= "lightpink", high= "lightblue",)+ ggtitle("Temperature-Water")
```

## NaCl-Wat:

``` r
condicion_NaCl_Water<- list(A= rownames(toptagsNaCl),
                 B=rownames(toptagsWat),
                 C=rownames(toptagsNaCl_Wat))

ggVennDiagram(condicion_NaCl_Water,label= "count",label_alpha = 0, category.names = c("NaCl","Wat","NaCl-Wat"),set_size = 2)+
  scale_fill_gradient(low= "lightpink", high= "lightblue",)+ ggtitle("NaCl-Water")
```

## NaCl-Bact:

``` r
condicion_NaCl_Bact<- list(A= rownames(toptagsNaCl),
                 B=rownames(toptagsBact),
                 C=rownames(toptagsNaCl_Bact))

ggVennDiagram(condicion_NaCl_Bact,label= "count",label_alpha = 0, category.names = c("NaCl","Bact","NaCl-Bact"),set_size = 2)+
  scale_fill_gradient(low= "lightpink", high= "lightblue",)+ ggtitle("NaCl-Bacterium")
```

## Bact-Wat:

``` r
condicion_Bact_Wat<- list(A= rownames(toptagsBact),
                 B=rownames(toptagsWat),
                 C=rownames(toptagsBact_Wat))
ggVennDiagram(condicion_Bact_Wat,label= "count",label_alpha = 0, category.names = c("Bact","Wat","Bact-Wat"),set_size = 2)+
  scale_fill_gradient(low= "lightpink", high= "lightblue",)+ ggtitle("Bacterium-Water")
```

# GO :

## Packages:

``` r
BiocManager::install("topGO")
```

## Libraries :

``` r
library(topGO)
```

## Funtion for plots:

``` r
plottopGo<-function(goEnrichment, ntop=20){

goEnrichment$classic <- as.numeric(goEnrichment$classic)
goEnrichment <- goEnrichment[goEnrichment$classic < 0.05,] # filter terms for Fisher p<0.05
goEnrichment <- goEnrichment[,c("GO.ID","Term","classic")]
goEnrichment


ggdata <- goEnrichment[1:ntop,]
ggdata<-ggdata[!duplicated(ggdata$Term),]
ggdata$Term <- factor(ggdata$Term, levels = rev(ggdata$Term)) # fixes order
gg1 <- ggplot(ggdata,
              aes(x = Term, y = -log10(classic), size = -log10(classic), fill = -log10(classic))) +
  
  expand_limits(y = 1) +
  geom_point(shape = 21) +
  scale_size(range = c(1.5,5.5)) +
  scale_fill_continuous(low = 'royalblue', high = 'red4') +
  
  xlab('') + ylab('Enrichment score') +
  labs(
    title = 'GO category',
    subtitle = 'Top 20 terms ordered by p-value',
    caption = 'Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001') +
  
  geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
             linetype = c("dotted", "longdash", "solid"),
             colour = c("black", "black", "black"),
             linewidth = c(0.5, 0.5, 0.5)) +
  
  theme_bw(base_size = 10) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 10, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 10, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 10, face = 'bold', vjust = 1),
    
    axis.text.x = element_text(angle = 0, size = 10, face = 'bold', hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 10, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 10, face = 'bold'),
    axis.title.x = element_text(size = 10, face = 'bold'),
    axis.title.y = element_text(size = 10, face = 'bold'),
    axis.line = element_line(colour = 'black'),
    
    #Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(0.5, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 10, face = "bold"), # Text size
    title = element_text(size = 10, face = "bold")) +
  
  coord_flip()
gg1
}
```

### CO2

``` r
geneID2GOCO2 <- readMappings("Salicornia_europaea_annotation_sint.txt")
ErodiumGenes<-names(geneID2GOCO2) # gene names

genes_CO2 <- rownames(toptagsCo2)
geneListCO2 <- factor(as.integer(ErodiumGenes %in% genes_CO2)) # make a factor object with selected and non-selected genes

names(geneListCO2) <- ErodiumGenes


GOdataCO2 <-new ("topGOdata", 
              ontology = "BP", #change BP, ML or
              allGenes = geneListCO2, 
              nodeSize = 7, 
              annot = annFUN.gene2GO,
              gene2GO = geneID2GOCO2)
              
resultFisCO2 <- runTest(GOdataCO2, algorithm = "classic", statistic = "fisher") # run Fisher's test

allResCO2 <- GenTable(GOdataCO2, classic = resultFisCO2,topNodes = 100) # complete list of enriched GO terms; topNodes = 30 for top30
write.table(allResCO2,file = 'GOCO2.txt', row.names = F, sep="\t",quote = FALSE)
```

#### Plot:

``` r
plottopGo(allResCO2)
```

### Temp

``` r
geneID2GOTemp <- readMappings("Salicornia_europaea_annotation_sint.txt")
ErodiumGenes<-names(geneID2GOTemp) # gene names

genes_Temp <- rownames(toptagsTemp)
geneListTemp <- factor(as.integer(ErodiumGenes %in% genes_Temp)) # make a factor object with selected and non-selected genes

names(geneListTemp) <- ErodiumGenes


GOdataTemp <-new ("topGOdata", 
                  ontology = "BP", #change BP, ML or
                  allGenes = geneListTemp, 
                  nodeSize = 5, 
                  annot = annFUN.gene2GO,
                  gene2GO = geneID2GOTemp)

resultFisTemp <- runTest(GOdataTemp, algorithm = "classic", statistic = "fisher") # run Fisher's test

allResTemp<- GenTable(GOdataTemp, classic = resultFisTemp, topNodes = 100) # complete list of enriched GO terms; c for top30
write.table(allResTemp,file = 'GOTemp.txt', row.names = F, sep="\t",quote = FALSE)
```

#### Plot

``` r
plottopGo(allResTemp)
```

### NaCl

``` r
geneID2GONaCl <- readMappings("Salicornia_europaea_annotation_sint.txt")
ErodiumGenes<-names(geneID2GONaCl) # gene names

genes_NaCl <- rownames(toptagsNaCl)
geneListNaCl <- factor(as.integer(ErodiumGenes %in% genes_NaCl)) # make a factor object with selected and non-selected genes

names(geneListNaCl) <- ErodiumGenes


GOdataNaCl <-new ("topGOdata", 
                  ontology = "BP", #change BP, ML or
                  allGenes = geneListNaCl, 
                  nodeSize = 5, 
                  annot = annFUN.gene2GO,
                  gene2GO = geneID2GONaCl)

resultFisNaCl <- runTest(GOdataNaCl, algorithm = "classic", statistic = "fisher") # run Fisher's test

allResNaCl<- GenTable(GOdataNaCl, classic = resultFisNaCl, topNodes = 100)
write.table(allResNaCl,file = 'GONaCl.txt', row.names = F, sep="\t",quote = FALSE)
```

#### Plot

``` r
plottopGo(allResNaCl)
```

##### NaCl downegulated

``` r
filas_negativoNaCl <- toptagsNaCl[toptagsNaCl$logFC < 0, ]
dim(filas_negativoNaCl)
downregulated = rownames(filas_negativoNaCl)



geneID2GONaCldown <- readMappings("Salicornia_europaea_annotation_sint.txt")
ErodiumGenes<-names(geneID2GONaCldown) # gene names

genes_NaCldown <- downregulated
geneListNaCldown <- factor(as.integer(ErodiumGenes %in% genes_NaCldown)) # make a factor object with selected and non-selected genes

names(geneListNaCldown) <- ErodiumGenes


GOdataNaCldown <-new ("topGOdata", 
                    ontology = "BP", #change BP, ML or
                    allGenes = geneListNaCldown, 
                    nodeSize = 5, 
                    annot = annFUN.gene2GO,
                    gene2GO = geneID2GONaCldown)

resultFisNaCldown <- runTest(GOdataNaCldown, algorithm = "classic", statistic = "fisher") # run Fisher's test

allResNaCldown<- GenTable(GOdataNaCldown, classic = resultFisNaCldown, topNodes = 100) # complete list of enriched GO terms; c for top30
write.table(allResNaCldown,file = 'GONaCldown.txt', row.names = F, sep="\t",quote = FALSE)
```

##### NaCl upregulated

``` r
filas_positivasNaCl <- toptagsNaCl[toptagsNaCl$logFC > 0, ]
dim(filas_positivasNaCl)
Upregulated = rownames(filas_positivasNaCl)


geneID2GONaClup <- readMappings("Salicornia_europaea_annotation_sint.txt")
ErodiumGenes<-names(geneID2GONaClup) # gene names

genes_NaClup <- Upregulated
geneListNaClup <- factor(as.integer(ErodiumGenes %in% genes_NaClup)) # make a factor object with selected and non-selected genes

names(geneListNaClup) <- ErodiumGenes


GOdataNaClup <-new ("topGOdata", 
                  ontology = "BP", #change BP, ML or
                  allGenes = geneListNaClup, 
                  nodeSize = 5, 
                  annot = annFUN.gene2GO,
                  gene2GO = geneID2GONaClup)

resultFisNaClup <- runTest(GOdataNaClup, algorithm = "classic", statistic = "fisher") # run Fisher's test

allResNaClup<- GenTable(GOdataNaClup, classic = resultFisNaClup, topNodes = 100) # complete list of enriched GO terms; c for top30
write.table(allResNaClup,file = 'GONaClup.txt', row.names = F, sep="\t",quote = FALSE)




filas_positivasNaCl <- toptagsNaCl[toptagsNaCl$logFC > 0, ]
dim(filas_positivasNaCl)
Upregulated = rownames(filas_positivasNaCl)


geneID2GONaClup <- readMappings("Salicornia_europaea_annotation_sint.txt")
ErodiumGenes<-names(geneID2GONaClup) # gene names

genes_NaClup <- Upregulated
geneListNaClup <- factor(as.integer(ErodiumGenes %in% genes_NaClup)) # make a factor object with selected and non-selected genes

names(geneListNaClup) <- ErodiumGenes


GOdataNaClup <-new ("topGOdata", 
                  ontology = "BP", #change BP, ML or
                  allGenes = geneListNaClup, 
                  nodeSize = 5, 
                  annot = annFUN.gene2GO,
                  gene2GO = geneID2GONaClup)

resultFisNaClup <- runTest(GOdataNaClup, algorithm = "classic", statistic = "fisher") # run Fisher's test

allResNaClup<- GenTable(GOdataNaClup, classic = resultFisNaClup, topNodes = 100) # complete list of enriched GO terms; c for top30
write.table(allResNaClup,file = 'GONaClup.txt', row.names = F, sep="\t",quote = FALSE)
```

### Bact

``` r
geneID2GOBact <- readMappings("Salicornia_europaea_annotation_sint.txt")
ErodiumGenes<-names(geneID2GOBact) # gene names

genes_Bact <- topBact
geneListBact <- factor(as.integer(ErodiumGenes %in% genes_Bact)) # make a factor object with selected and non-selected genes
names(geneListBact) <- ErodiumGenes


GOdataBact <-new ("topGOdata", 
                  ontology = "BP", #change BP, ML or
                  allGenes = geneListBact, 
                  nodeSize = 5, 
                  annot = annFUN.gene2GO,
                  gene2GO = geneID2GOBact)

resultFisBact <- runTest(GOdataBact, algorithm = "classic", statistic = "fisher") # run Fisher's test

allResBact<- GenTable(GOdataBact, classic = resultFisBact,topNodes = 100)
write.table(allResBact,file = 'GOBact.txt', row.names = F, sep="\t",quote = FALSE)
```

#### Plot

``` r
plottopGo(allResBact)
```

### Wat

``` r
geneID2GOWat <- readMappings("Salicornia_europaea_annotation_sint.txt")
ErodiumGenes<-names(geneID2GOWat) # gene names

genes_Wat <- rownames(toptagsWat)
geneListWat <- factor(as.integer(ErodiumGenes %in% genes_Wat)) # make a factor object with selected and non-selected genes


names(geneListWat) <- ErodiumGenes


GOdataWat <-new ("topGOdata", 
                  ontology = "BP", #change BP, ML or
                  allGenes = geneListWat, 
                  nodeSize = 5, 
                  annot = annFUN.gene2GO,
                  gene2GO = geneID2GOWat)

resultFisWat <- runTest(GOdataWat, algorithm = "classic", statistic = "fisher") # run Fisher's test

allResWat<- GenTable(GOdataWat, classic = resultFisWat ,topNodes = 100)
write.table(allResWat,file = 'GOWat.txt', row.names = F, sep="\t",quote = FALSE)


write.table(allResWat[, c(1, 6)],file = 'GOWatp-valor.txt', row.names = F, sep="\t",quote = FALSE)
```

#### Plot

``` r
plottopGo(allResWat)
```

### Temp-Bact

``` r
geneID2GOTempBact <- readMappings("Salicornia_europaea_annotation_sint.txt")
ErodiumGenes<-names(geneID2GOTempBact) # gene names

genes_TempBact <- rownames(toptagsTemp_Bact)
geneListTempBact <- factor(as.integer(ErodiumGenes %in% genes_TempBact)) # make a factor object with selected and non-selected genes

names(geneListTempBact) <- ErodiumGenes


GOdataTempBact <-new ("topGOdata", 
                 ontology = "BP", #change BP, ML or
                 allGenes = geneListTempBact, 
                 nodeSize = 5, 
                 annot = annFUN.gene2GO,
                 gene2GO = geneID2GOTempBact)

resultFisTempBact <- runTest(GOdataTempBact, algorithm = "classic", statistic = "fisher") # run Fisher's test

allResTempBact <- GenTable(GOdataTempBact, classic = resultFisTempBact,topNodes = 100) # complete list of enriched GO terms; topNodes = 30 for top30
write.table(allResTempBact,file = 'GOTempBact.txt', row.names = F, sep="\t",quote = FALSE)
```

#### Plot:

``` r
plottopGo(allResCO2)
```

### CO2-Temp

``` r
geneID2GOCO2_Temp <- readMappings("Salicornia_europaea_annotation_sint.txt")
ErodiumGenes<-names(geneID2GOCO2_Temp) # gene names

genes_CO2_Temp <- rownames(toptagsCO2_Temp)
geneListCO2_Temp <- factor(as.integer(ErodiumGenes %in% genes_CO2_Temp)) # make a factor object with selected and non-selected genes

names(geneListCO2_Temp) <- ErodiumGenes


GOdataCO2_Temp <-new ("topGOdata", 
                  ontology = "BP", #change BP, ML or
                  allGenes = geneListCO2_Temp, 
                  nodeSize = 5, 
                  annot = annFUN.gene2GO,
                  gene2GO = geneID2GOCO2_Temp)

resultFisCO2_Temp <- runTest(GOdataCO2_Temp, algorithm = "classic", statistic = "fisher") # run Fisher's test

allResCO2_Temp<- GenTable(GOdataCO2_Temp, classic = resultFisCO2_Temp, topNodes = 1000)
write.table(allResCO2_Temp,file = 'GOCO2_Temp.txt', row.names = F, sep="\t",quote = FALSE)
```

#### Plot

``` r
plottopGo(allResCO2_Temp)
```

``` r
save.image("UltimosresultadosT-C_sinred.RData")
```

# Cytoscape network but with only the DE:

## Install Packages:

``` r
#install.packages("WGCNA")
```

## Library:

``` r
library(WGCNA)
```

## Selec the DE:

``` r
Topgenes <- c(rownames(toptagsCo2),rownames(toptagsTemp),rownames(toptagsWat),rownames(toptagsBact),rownames(toptagsCO2_Temp),rownames(toptagsCO2_Bact),rownames(toptagsCO2_Wat),rownames(toptagsCO2_NaCl),rownames(toptagsTemp_NaCl),rownames(toptagsTemp_Wat),rownames(toptagsTemp_Bact),rownames(toptagsNaCl_Bact),rownames(toptagsNaCl_Wat),rownames(toptagsBact_Wat),rownames(toptagsNaCl))

Topgenes<-unique(Topgenes)
n.counts.topgenes <- cpm(y)[Topgenes,]
n.counts.topgenes = t(n.counts.topgenes)
dim(n.counts.topgenes)
```

## umbral:

``` r
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(n.counts.topgenes, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1);
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90)

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1)
```

``` r
net.topgenes = blockwiseModules(n.counts.topgenes, power = 12,
TOMType = "unsigned",networkType = "unsigned", minModuleSize = 30,
reassignThreshold = 0, mergeCutHeight = 0.25,
numericLabels = TRUE, pamRespectsDendro = FALSE,
saveTOMs = TRUE,
saveTOMFileBase = "DE_Salicornia_TOM",
verbose = 3)
```

## Modules:

``` r
 table(net.topgenes$colors)
```

## save :

``` r
moduleLabels = net.topgenes$colors
moduleColors = labels2colors(net.topgenes$colors)
MEs = net.topgenes$MEs;
geneTree = net.topgenes$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
file = "DE_Salicornia_Red.RData")
```

## Dendrograma:

``` r
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net.topgenes$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net.topgenes$dendrograms[[1]], mergedColors[net.topgenes$blockGenes[[1]]],
"Module colors",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
```

``` r
modulos <-unique( moduleColors)
```

``` r
# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(n.counts.topgenes, power = 12,
TOMType = "signed")
# Select modules :
modules = modulos ;
# Select module probes
probes = colnames(n.counts.topgenes)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
weighted = TRUE,
threshold = 0.02,
nodeNames = modProbes,
nodeAttr = moduleColors[inModule])
```

## Hubs:

``` r
colorh <- modulos
hubs= chooseTopHubInEachModule(
   n.counts.topgenes, 
   net.topgenes$colors, 
   power = 12, 
   type = "unsigned")
write.table(hubs,file = 'hubs.txt', row.names = F, sep="\t",quote = FALSE)
```

## Analysis of GO terms for modules and hubs:

### Hubs :

``` r
geneID2GOhub <- readMappings("Salicornia_europaea_annotation_sint.txt")
ErodiumGenes<-names(geneID2GOhub) # gene names


genes_hub <- c("Seu_jg10511", "Seu_jg25347", "Seu_jg22365", "Seu_jg23957", "Seu_jg209", "Seu_jg5477", "Seu_jg2688", "Seu_jg10672", "Seu_jg1849", "Seu_jg2158", "Seu_jg8570", "Seu_jg3126", "Seu_jg20981", "Seu_jg8737", "Seu_jg2401", "Seu_jg3558", "Seu_jg10959", "Seu_jg10115", "Seu_jg10262", "Seu_jg18573", "Seu_jg28392")

geneListhub <- factor(as.integer(ErodiumGenes %in% genes_hub)) # make a factor object with selected and non-selected genes

names(geneListhub) <- ErodiumGenes


GOdatahub <-new ("topGOdata", 
              ontology = "BP", #change BP, ML or
              allGenes = geneListhub, 
              nodeSize = 7, 
              annot = annFUN.gene2GO,
              gene2GO = geneID2GOhub)
              
resultFishub <- runTest(GOdatahub, algorithm = "classic", statistic = "fisher") # run Fisher's test

allReshub <- GenTable(GOdatahub, classic = resultFishub,topNodes = 100) # complete list of enriched GO terms; topNodes = 30 for top30
write.table(allReshub,file = 'GOhub.txt', row.names = F, sep="\t",quote = FALSE)
```

###Module 1 :

``` r
acoloeess=net.topgenes$colors[which(net.topgenes$colors == 1)]
modulo1=names(acoloeess)
```

``` r
geneID2GOmodulo1 <- readMappings("Salicornia_europaea_annotation_sint.txt")
ErodiumGenes<-names(geneID2GOmodulo1) # gene names


genes_modulo1 <- modulo1

geneListmodulo1 <- factor(as.integer(ErodiumGenes %in% genes_modulo1)) # make a factor object with selected and non-selected genes

names(geneListmodulo1) <- ErodiumGenes


GOdatamodulo1 <-new ("topGOdata", 
                 ontology = "BP", #change BP, ML or
                 allGenes = geneListmodulo1, 
                 nodeSize = 7, 
                 annot = annFUN.gene2GO,
                 gene2GO = geneID2GOmodulo1)

resultFismodulo1 <- runTest(GOdatamodulo1, algorithm = "classic", statistic = "fisher") # run Fisher's test

allResmodulo1 <- GenTable(GOdatamodulo1, classic = resultFismodulo1,topNodes = 100) # complete list of enriched GO terms; topNodes = 30 for top30
write.table(allResmodulo1,file = 'GOmodulo1.txt', row.names = F, sep="\t",quote = FALSE)
```

### Module 2 :

``` r
acoloeess=net.topgenes$colors[which(net.topgenes$colors == 2)]
modulo2=names(acoloeess)


geneID2GOmodulo2 <- readMappings("Salicornia_europaea_annotation_sint.txt")
ErodiumGenes<-names(geneID2GOmodulo2) # gene names


genes_modulo2 <- modulo2

geneListmodulo2 <- factor(as.integer(ErodiumGenes %in% genes_modulo2)) # make a factor object with selected and non-selected genes

names(geneListmodulo2) <- ErodiumGenes


GOdatamodulo2 <-new ("topGOdata", 
                     ontology = "BP", #change BP, ML or
                     allGenes = geneListmodulo2, 
                     nodeSize = 7, 
                     annot = annFUN.gene2GO,
                     gene2GO = geneID2GOmodulo2)

resultFismodulo2 <- runTest(GOdatamodulo2, algorithm = "classic", statistic = "fisher") # run Fisher's test

allResmodulo2 <- GenTable(GOdatamodulo2, classic = resultFismodulo2,topNodes = 100) # complete list of enriched GO terms; topNodes = 30 for top30
write.table(allResmodulo2,file = 'GOmodulo2.txt', row.names = F, sep="\t",quote = FALSE)
```

### Module 3 :

``` r
acoloeess=net.topgenes$colors[which(net.topgenes$colors == 3)]
modulo3=names(acoloeess)


geneID2GOmodulo3 <- readMappings("Salicornia_europaea_annotation_sint.txt")
ErodiumGenes<-names(geneID2GOmodulo3) # gene names


genes_modulo3 <- modulo3

geneListmodulo3 <- factor(as.integer(ErodiumGenes %in% genes_modulo3)) # make a factor object with selected and non-selected genes

names(geneListmodulo3) <- ErodiumGenes


GOdatamodulo3 <-new ("topGOdata", 
                     ontology = "BP", #change BP, ML or
                     allGenes = geneListmodulo3, 
                     nodeSize = 7, 
                     annot = annFUN.gene2GO,
                     gene2GO = geneID2GOmodulo3)

resultFismodulo3 <- runTest(GOdatamodulo3, algorithm = "classic", statistic = "fisher") # run Fisher's test

allResmodulo3 <- GenTable(GOdatamodulo3, classic = resultFismodulo3,topNodes = 100) # complete list of enriched GO terms; topNodes = 30 for top30
write.table(allResmodulo3,file = 'GOmodulo3.txt', row.names = F, sep="\t",quote = FALSE)
```

# Creating multi-panel charts

## install packages

``` r
install.packages("grid")
install.packages("gridExtra")
install.packages("ggplotify")
install.packages("imager")
```

## Libraries

``` r
library(grid)
library(gridExtra)
library(ggplotify)
library(imager)
```

## PCA y Barplot

``` r
g1=ggplot(data, aes(x = X.PC1, y = X.PC2, color = water)) +
  geom_point() +
  geom_text(aes(label = muestras), size = 3) +  
  labs(x = "PC1 (36%)", y = "PC2 (19%)", title = "") +
  theme_bw() +
  guides(color = guide_legend(title = NULL, override.aes = list(size = 4))) 

tratamientos <- data.DE[,2]
upregulados <- c(557,709,312,0,1695,1180,0,4,453,0,32,231,0,8,0)
downregulated <- c(502,1022,188,4,3405,1642,3,2,424,0,9,287,0,12,0)
datos <- data.frame(tratamientos, upregulados, downregulated)

g2=ggplot(datos, aes(x = tratamientos, y = upregulados, fill = "Upregulados")) +
  geom_bar(stat = "identity") +
  geom_bar(aes(y = -downregulated, fill = "Downregulated"), stat = "identity", position = "identity") +
  geom_text(aes(label = upregulados), vjust = -0.5, color = "black", position = position_stack(vjust = 0.5), size = 3) +
  geom_text(aes(label = downregulated), vjust = 1.5, color = "black", position = position_stack(vjust = -0.5), size = 3) +
  labs(x = "Tratamientos", y = "Number of genes", fill = "") +
  scale_fill_manual(values = c("Upregulados" = "cornflowerblue", "Downregulated" = "cadetblue2")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_y_continuous(labels = scales::comma_format(scale = 0.001))

letras <- c("A", "B")
grid.arrange(g1, g2, nrow = 1)
grid.text(letras, x = c(0.02, 0.5), y = 0.95, gp = gpar(fontsize = 12, fontface = "bold"))
```

## Water

``` r
g1= ggplot(Wat_volcano, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = ifelse(logFC > log2(1.5), ifelse(FDR < 0.05  , "red", "grey"), ifelse(logFC < -log2(1.5), ifelse(FDR < 0.05  , "blue", "grey"), "grey")))) +
  xlab("logFC") +
  ylab("-log10(FDR)") +
  scale_color_identity() +
  theme_minimal()

g2= plottopGo(allResWat)

g3= ggVennDiagram(condicion_CO2_Wat,label= "count",label_alpha = 0, category.names = c("CO2","Wat","C_W"), set_size = 2)+
  scale_fill_gradient(low= "lightpink", high= "lightblue",) + ggtitle("CO2-Water")
g4= ggVennDiagram(condicion_Temp_Wat,label= "count",label_alpha = 0, category.names = c("Temp","Wat","T_W"), set_size = 2)+
  scale_fill_gradient(low= "lightpink", high= "lightblue",)+ ggtitle("Temperature-Water")
g5= ggVennDiagram(condicion_NaCl_Water,label= "count",label_alpha = 0, category.names = c("NaCl","Wat","N-W"),set_size = 2)+
  scale_fill_gradient(low= "lightpink", high= "lightblue",)+ ggtitle("NaCl-Water")
g6=ggVennDiagram(condicion_Bact_Wat,label= "count",label_alpha = 0, category.names = c("Bact","Wat","B-W"),set_size = 2)+
  scale_fill_gradient(low= "lightpink", high= "lightblue",)+ ggtitle("Bacterium-Water")

heatmap1 = pheatmap(logCPMWat, scale = "row",show_rownames = FALSE, show_colnames = FALSE)
g7 = as.grob(heatmap1)

gl= list(g1,g2,g3,g4,g5,g6,g7)

grid.arrange(
  grobs = gl,
  widths = c(1, 1, 1, 1),
  layout_matrix = rbind(c(1, 1,2,2),c(1,1,2,2),c(7,7,3,4),c(7,7,5,6)))

letras2 = c("A", "C")
letras3 = c("B","D", "E")
letras4= c("F","G")

grid.text(letras2, x = c(0.01, 0.51), y = 0.97, gp = gpar(fontsize = 12, fontface = "bold"))
grid.text(letras3, x = c(0.01, 0.5, 0.75), y = 0.47, gp = gpar(fontsize = 12, fontface = "bold"))
grid.text(letras4, x = c(0.5,0.75), y = 0.25, gp = gpar(fontsize = 12, fontface = "bold"))
```

## Temperature

``` r
g1= ggplot(Temp_volcano, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = ifelse(logFC > log2(1.5), ifelse(FDR < 0.05  , "red", "grey"), ifelse(logFC < -log2(1.5), ifelse(FDR < 0.05  , "blue", "grey"), "grey")))) +
  xlab("logFC") +
  ylab("-log10(FDR)") +
  scale_color_identity() +
  theme_minimal()

g2= plottopGo(allResTemp)

g3= ggVennDiagram(condicion_CO2_Temp, label= "count",label_alpha = 0, category.names = c("CO2","Temp","C_T"), set_size = 2)+
  scale_fill_gradient(low= "lightpink", high= "lightblue",)+ ggtitle("CO2-Temperature") 
g4=  ggVennDiagram(condicion_Temp_NaCl,label= "count",label_alpha = 0, category.names = c("Temp","NaCl","T_N"),set_size = 2)+
  scale_fill_gradient(low= "lightpink", high= "lightblue",) + ggtitle("Temperature-NaCl")
g5=ggVennDiagram(condicion_Temp_Bact,label= "count",label_alpha = 0, category.names = c("Temp","Bact","T_B"), set_size = 2)+
  scale_fill_gradient(low= "lightpink", high= "lightblue",)+ ggtitle("Temperature-Bacterium") 
heatmap1 = pheatmap(logCPMTemp, scale = "row" ,show_rownames = FALSE, show_colnames = FALSE)
g7 = as.grob(heatmap1)

gl= list(g1,g2,g3,g4,g5,g7)

grid.arrange(
  grobs = gl,
  widths = c(1, 1, 1, 1),
  layout_matrix = rbind(c(1, 1,2,2),c(1,1,2,2),c(7,7,3,4),c(7,7,5,NaN)))

letras2 = c("A", "C")
letras3 = c("B","D", "E")
letras4= c("F")

grid.text(letras2, x = c(0.01, 0.53), y = 0.97, gp = gpar(fontsize = 12, fontface = "bold"))
grid.text(letras3, x = c(0.01, 0.51, 0.76), y = 0.47, gp = gpar(fontsize = 12, fontface = "bold"))
grid.text(letras4, x = c(0.51), y = 0.26, gp = gpar(fontsize = 12, fontface = "bold"))
```

## CO2

``` r
g1= ggplot(CO2_volcano, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = ifelse(logFC > log2(1.5), ifelse(FDR < 0.05  , "red", "grey"), ifelse(logFC < -log2(1.5), ifelse(FDR < 0.05  , "blue", "grey"), "grey")))) +
  xlab("logFC") +
  ylab("-log10(FDR)") +
  scale_color_identity() +
  theme_minimal()
g2= plottopGo(allResCO2)

g4= ggVennDiagram(condicion_CO2_NaCl,label= "count",label_alpha = 0, category.names = c("CO2","NaCl","C_N"), set_size=2 )+
  scale_fill_gradient(low= "lightpink", high= "lightblue",)+ ggtitle("CO2-NaCl") 
g5= ggVennDiagram(condicion_CO2_Bact,label= "count",label_alpha = 0, category.names = c("CO2","Bact","C_B"), set_size = 2)+
  scale_fill_gradient(low= "lightpink", high= "lightblue",) + ggtitle("CO2-Bacterium")


heatmap1 = pheatmap(logCPMco2, scale = "row" ,show_rownames = FALSE, show_colnames = FALSE)
g7 = as.grob(heatmap1)


gl= list(g1,g2,g4,g5,g7)
grid.arrange(
  grobs = gl,
  widths = c(1, 1, 1, 1),
  layout_matrix = rbind(c(1,1,2,2),c(1,1,2,2),c(7,7,4,5),c(7,7,4,5)))

letras2 = c("A", "C")
letras3 = c("B","D", "E")


grid.text(letras2, x = c(0.01, 0.51), y = 0.97, gp = gpar(fontsize = 12, fontface = "bold"))
grid.text(letras3, x = c(0.01, 0.49, 0.74), y = 0.47, gp = gpar(fontsize = 12, fontface = "bold"))
```

## NaCl

``` r
g1=ggplot(NaCl_volcano, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = ifelse(logFC > log2(1.5), ifelse(FDR < 0.05  , "red", "grey"), ifelse(logFC < -log2(1.5), ifelse(FDR < 0.05  , "blue", "grey"), "grey")))) +
  xlab("logFC") +
  ylab("-log10(FDR)") +
  scale_color_identity() +
  theme_minimal()
g2= plottopGo(allResNaCl)

g5= ggVennDiagram(condicion_NaCl_Bact,label= "count",label_alpha = 0, category.names = c("NaCl","Bact","N-B"),set_size = 2)+
  scale_fill_gradient(low= "lightpink", high= "lightblue",)+ ggtitle("NaCl-Bacterium")
heatmap1 = pheatmap(logCPMNaCl, scale = "row" ,show_rownames = FALSE, show_colnames = FALSE)
g7 = as.grob(heatmap1)



x11()
gl= list(g1,g2,g5,g7)
grid.arrange(
  grobs = gl,
  widths = c(1, 1, 1, 1),
  layout_matrix = rbind(c(1,1,2,2),c(1,1,2,2),c(7,7,5,5),c(7,7,5,5)))

letras2 = c("A", "C")
letras3 = c("B","D")


grid.text(letras2, x = c(0.01, 0.51), y = 0.97, gp = gpar(fontsize = 12, fontface = "bold"))
grid.text(letras3, x = c(0.01, 0.51), y = 0.5, gp = gpar(fontsize = 12, fontface = "bold"))
```

## Bact

``` r
g1= ggplot(bact_volcano, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = ifelse(logFC > log2(1.5), ifelse(FDR < 0.05  , "red", "grey"), ifelse(logFC < -log2(1.5), ifelse(FDR < 0.05  , "blue", "grey"), "grey")))) +
  xlab("logFC") +
  ylab("-log10(FDR)") +
  scale_color_identity() +
  theme_minimal()
g2= plottopGo(allResBact)

heatmap1 = pheatmap(logCPMBact, scale = "row", show_rownames = FALSE, show_colnames = FALSE)
g7 = as.grob(heatmap1)

gl= list(g1,g2,g7)
grid.arrange(
  grobs = gl,
  widths = c(1, 1, 1, 1),
  layout_matrix = rbind(c(1,1,2,2),c(1,1,2,2),c(7,7,2,2),c(7,7,2,2)))

letras2 = c("A", "C")
letras3 = c("B")


grid.text(letras2, x = c(0.01, 0.51), y = 0.97, gp = gpar(fontsize = 12, fontface = "bold"))
grid.text(letras3, x = c(0.01), y = 0.5, gp = gpar(fontsize = 12, fontface = "bold"))
```

## Plots

``` r
g2=plottopGo(allResCO2_Temp)
g1=plottopGo(allResTempBact)
letras <- c("A", "B")
grid.arrange(g1, g2, nrow = 1)
grid.text(letras, x = c(0.02, 0.5), y = 0.95, gp = gpar(fontsize = 12, fontface = "bold"))
```

## Plots

``` r
 g1= load.image("C:/Users/Usuario/Desktop/TFM/star T-C/WGCNSalicornia3.png")
 g2= load.image("C:/Users/Usuario/Desktop/TFM/star T-C/cluster2.jpeg")
 x11()
 par(mfrow = c(1, 2))
plot(g1)
plot(g2)

letras <- c("A", "B")
grid.text(letras, x = c(0.02, 0.5), y = 0.95, gp = gpar(fontsize = 12, fontface = "bold"))
```

``` r
save.image("UltimosresultadosT-C.RData")
```
