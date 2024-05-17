# NasalMicrobiome
Micorbiome Immunopolyp



## Introduction

In the IMMUNOPOLYP study we sampled 27 patients with chronic rhinosinusitis treated with dupilumab 300mg s.c. every other week. We collected nasal swabs and stool samples during the course of treatment. We performed 16S shotgun-metagenomics analysis.

## Load packages

Prior to the metagenomics data analysis we load all necessary packages.

```{r}
#| output: false
library(openxlsx)
library(DECIPHER)
library(dada2)
library(readr)
library(seqinr)
library(phyloseq)
library(vegan)
library(ggplot2)
library(viridis)
library(dplyr)
library(gss)
library(coda4microbiome)
library(pheatmap)
library(ggpubr)
library(microbiome)
library(mia)
library(ggrepel)
library(metagenomeSeq)
library(MicrobiomeStat)
library(RColorBrewer)
library(DESeq2)
library(rstatix)
library(tidyr)
```

## Data import

For further information about cleaning data and the preliminary analysis of ASV signals please refer to [DADA2 package](https://www.bioconductor.org/packages/release/bioc/html/dada2.html) and the corresponding tutorials on the [github page](https://benjjneb.github.io/dada2/tutorial.html).
The data provided in this example have already been cleaned and the taxonomy assignement has been conducted.
We start by load the packaged ASV counts without contamination and the ASV taxonomy as well as the sample metadata.

```{r}
#| output: false
# load our ASV count data from the shotgun metagenomic analysis
count_tab <- read.table("//filer300/USERS3001/I0331555/Documents/IMMUNOPOLYP/MIKROPOLYP/ASVs_counts-no-contam.tsv", header=T, row.names=1,
                        check.names=F, sep="\t")

# load our taxnomy table (output of the assinging taxonomy according to DADA2 package)
tax_tab <- as.matrix(read.table("//filer300/USERS3001/I0331555/Documents/IMMUNOPOLYP/MIKROPOLYP/ASVs_taxonomy.tsv", header=T,
                                row.names=1, check.names=F, sep="\t"))
# match ASV to the taxonomy
ASV <- rownames(count_tab)
tax_tab <- cbind(tax_tab, ASV)
tax_tab <- as.data.frame(tax_tab)
rownames(tax_tab) <- tax_tab$ASV
# save as matrix for later
tax_tab <- as.matrix(tax_tab)

# load the sample information data (only nasal samples)
sample1 <- read.xlsx("//filer300/USERS3001/I0331555/Documents/IMMUNOPOLYP/MIKROPOLYP/sample_information_nasal.xlsx")

sample_info_tab <- read.table("//filer300/USERS3001/I0331555/Documents/IMMUNOPOLYP/MIKROPOLYP/sample_information_nasal.tsv", header=T, row.names=1,
                              check.names=F, sep="\t")

# load the sample information for stool samples
sample_info_stool <- read.table("//filer300/USERS3001/I0331555/Documents/IMMUNOPOLYP/MIKROPOLYP/sample_information_stool.tsv", header=T, row.names=1,
                              check.names=F, sep="\t")
sample2 <- read.xlsx("//filer300/USERS3001/I0331555/Documents/IMMUNOPOLYP/MIKROPOLYP/sample_information_stool.xlsx")

# merge all stool sample data to one dataframe
sample_info_stool$SampleID <- rownames(sample_info_stool)
sample_info_stool <- merge(sample_info_stool, sample2, by = "treat")
rownames(sample_info_stool) <- sample_info_stool$SampleID
sample_info_stool$Time <- factor(sample_info_stool$Time, levels = c("healthy", "0", "28", "90", "180"))

# merge all nasal sample data to one dataframe
sample_info_tab$SampleID <- rownames(sample_info_tab)
sample_info_tab <- merge(sample_info_tab, sample1, by = "treat")
rownames(sample_info_tab) <- sample_info_tab$SampleID
sample_info_tab$Time <- factor(sample_info_tab$Time, levels = c("healthy", "0", "28", "90", "180"))

# subset the count table to only nasal microbiome and exclude all stool samples
count_tab_nasal <- count_tab[, colnames(count_tab) %in% sample_info_tab$SampleID]

# subset the count table to only stool microbiome and exclude all stool samples
count_tab_stool <- count_tab[, colnames(count_tab) %in% sample_info_stool$SampleID]
count_tab_stool <- as.data.frame(count_tab_stool)

```

## Designing a phyloseq object

To proceed with the microbiome analysis we use the phyloseq package and design a phyloseq object with the count table, taxonomy table and the sample data.

```{r}
#| output: false
# assign the right tables to the right matrices
count_tab_phy <- otu_table(count_tab_nasal, taxa_are_rows=T)
tax_tab_phy <- tax_table(tax_tab)
sample_info_tab_phy <- sample_data(sample_info_tab)

# merge all table to a S4 object
ASV_physeq <- phyloseq(count_tab_phy, tax_tab_phy, sample_info_tab_phy)
```

We do the same for our stool samples. 

```{r}
#| output: false
# assign the right tables to the right matrices
count_tab_phy_stool <- otu_table(count_tab_stool, taxa_are_rows=T)
tax_tab_phy_stool <- tax_table(tax_tab)
sample_info_tab_phy_stool <- sample_data(sample_info_stool)

# merge all table to a S4 object
ASV_physeq_stool <- phyloseq(count_tab_phy_stool, tax_tab_phy_stool, sample_info_tab_phy_stool)
```

# Alpha Diversity

Because the cleaning of the data has already been done in the preliminary analysis using DADA2 package, we can start straight forward with the analysis of the alpha diversity following the [phyloseq workflow](https://joey711.github.io/phyloseq/plot_richness-examples.html) and the [microbiome worflow](https://microbiome.github.io/tutorials/).

```{r}
#| fig-label: diversity
#| fig-width: 8
#| fig-length: 20
#| fig-cap: "Shannon diversity in nasal microbiome in patients with diffuse type 2 chronic rhinosinusitis under dupilumab treatment and healthy controls."
# create boxplot including shannon diversity with wilcoxan rank test
p.shannon <- boxplot_alpha(ASV_physeq, 
                           index = "shannon",
                           x_var = "Time") +
  stat_compare_means(method = "wilcox", step.increase = 0.06, 
                     comparisons = list(c("healthy", "0"), c("healthy", "28"), c("healthy", "90"), c("healthy", "180"), c("0", "28"), c("0", "90"), c("0", "180"))) +
  scale_fill_viridis(discrete = TRUE, alpha = 1, option="viridis")+
  theme_bw()+
  ylab("Shannon-Diversity")+
  labs(title = "Nasal Microbiome")

p.chao <- boxplot_alpha(ASV_physeq, 
                           index = "chao1",
                           x_var = "Time") +
  stat_compare_means(method = "wilcox", step.increase = 0.06, 
                     comparisons = list(c("healthy", "0"), c("healthy", "28"), c("healthy", "90"), c("healthy", "180"), c("0", "28"), c("0", "90"), c("0", "180"))) +
  scale_fill_viridis(discrete = TRUE, alpha = 1, option="viridis")+
  theme_bw()+
  ylab("Chao1-Diversity")+
  labs(title = "Nasal Microbiome")
```

We do the same analysis for our stool samples. Then we arrange the two plots in one.

```{r}
#| fig-label: diversity
#| fig-width: 8
#| fig-length: 20
#| fig-cap: "Shannon diversity in nasal and gut microbiome in patients with diffuse type 2 chronic rhinosinusitis under dupilumab treatment and healthy controls."

# create boxplot including shannon diversity with wilcoxan rank test
p.shannon_stool <- boxplot_alpha(ASV_physeq_stool, 
                           index = "shannon",
                           x_var = "Time") +
  stat_compare_means(method = "wilcox", step.increase = 0.06, 
                     comparisons = list(c("healthy", "0"), c("healthy", "28"), c("healthy", "90"), c("healthy", "180"), c("0", "28"), c("0", "90"), c("0", "180"))) +
  scale_fill_viridis(discrete = TRUE, alpha = 1, option="viridis")+
  theme_bw()+
  ylab("Shannon-Diversity")+
  labs(title = "Intestinal Microbiome")

p.chao_stool <- boxplot_alpha(ASV_physeq_stool, 
                           index = "chao1",
                           x_var = "Time") +
  stat_compare_means(method = "wilcox", step.increase = 0.06, 
                     comparisons = list(c("healthy", "0"), c("healthy", "28"), c("healthy", "90"), c("healthy", "180"), c("0", "28"), c("0", "90"), c("0", "180"))) +
  scale_fill_viridis(discrete = TRUE, alpha = 1, option="viridis")+
  theme_bw()+
  ylab("Chao1-Diversity")+
  labs(title = "Intestinal Microbiome")

legend <- get_legend(p.shannon)
div <- ggarrange(p.shannon, p.chao, p.shannon_stool, p.chao_stool, nrow = 1, ncol = 4, legend = "none")
div
```

The alpha diversity decreases significantly under dupilumab therapy. 

# Beta diversity

Now we conduct beta-diversity as instructed on the github webpage [course_2021_radboud](https://microbiome.github.io/course_2021_radboud/beta-diversity.html#examples-of-pcoa-with-different-settings).

We use Bray-Curtis distance for this analysis.

```{r}
#| fig-width: 12
#| fig-length: 6
#| fig-cap: "Bray-Curtis analysis nasal"
# with the phyloseq package we design a beta diversity plot
bray <- ordinate(ASV_physeq, method = "PCoA", distance = "bray")
beta1 <- plot_ordination(physeq = ASV_physeq, ordination = bray, color = "Time")+
  geom_point(size = 4)+
  stat_ellipse()+
  scale_color_viridis(discrete = TRUE)+
  theme_bw()

bray_stool <- ordinate(ASV_physeq_stool, method = "PCoA", distance = "bray")
beta2 <- plot_ordination(physeq = ASV_physeq_stool, ordination = bray_stool, color = "Time")+
  geom_point(size = 4)+
  stat_ellipse()+
  scale_color_viridis(discrete = TRUE)+
  theme_bw()

bt <- ggarrange(beta1, beta2, ncol =2)
bt
```


```{r}
#| output: false
# now we calculate the bray distance to test using PERMANOVA testing with 999 permutations.
abrel_bray <- phyloseq::distance(ASV_physeq, method = "bray")
sam <- data.frame(sample_data(ASV_physeq))
abrel_bray <- as.matrix(abrel_bray)

# our permanova test for healthy vs. CRS
sam_CRS <- sam[!is.na(sam$CRS_healthy), ]
abrel_bray_CRS <- abrel_bray[rownames(abrel_bray) %in% sam_CRS$SampleID, colnames(abrel_bray) %in% sam_CRS$SampleID]
permanova <- adonis2(abrel_bray_CRS ~ CRS_healthy, data = sam_CRS)

# our permanova test for healthy vs. d28
sam_CRS_d28 <- sam[!is.na(sam$healtyh_d28), ]
abrel_bray_CRS_d28 <- abrel_bray[rownames(abrel_bray) %in% sam_CRS_d28$SampleID, colnames(abrel_bray) %in% sam_CRS_d28$SampleID]
permanova <- adonis2(abrel_bray_CRS_d28 ~ healtyh_d28, data = sam_CRS_d28)

# our permanova test for healthy vs. d90
sam_CRS_d90 <- sam[!is.na(sam$healthy_d90), ]
abrel_bray_CRS_d90 <- abrel_bray[rownames(abrel_bray) %in% sam_CRS_d90$SampleID, colnames(abrel_bray) %in% sam_CRS_d90$SampleID]
permanova <- adonis2(abrel_bray_CRS_d90 ~ healthy_d90, data = sam_CRS_d90)

# our permanova test for healthy vs. d180
sam_180 <- sam[!is.na(sam$healthy_d180), ]
abrel_bray_180 <- abrel_bray[rownames(abrel_bray) %in% sam_180$SampleID, colnames(abrel_bray) %in% sam_180$SampleID]
permanova2 <- adonis2(abrel_bray_180 ~ healthy_d180, data = sam_180)

# our permanova test for d0 vs. d28
sam_28 <- sam[!is.na(sam$d0_d28), ]
abrel_bray_28 <- abrel_bray[rownames(abrel_bray) %in% sam_28$SampleID, colnames(abrel_bray) %in% sam_28$SampleID]
permanova2 <- adonis2(abrel_bray_28 ~ d0_d28, data = sam_28)

# our permanova test for d0 vs. d90
sam_90 <- sam[!is.na(sam$d0_d90), ]
abrel_bray_90 <- abrel_bray[rownames(abrel_bray) %in% sam_90$SampleID, colnames(abrel_bray) %in% sam_90$SampleID]
permanova2 <- adonis2(abrel_bray_90 ~ d0_d90, data = sam_90)

# our permanova test for d0 vs. d180
sam_180 <- sam[!is.na(sam$d0_180), ]
abrel_bray_180 <- abrel_bray[rownames(abrel_bray) %in% sam_180$SampleID, colnames(abrel_bray) %in% sam_180$SampleID]
permanova2 <- adonis2(abrel_bray_180 ~ d0_180, data = sam_180)
```

```{r}
#| output: false
# we do the same PERMANOVA analysis with the stool probes to test for similarity of the samples.

# now we calculate the bray distance to test using PERMANOVA testing with 999 permutations.
abrel_bray <- phyloseq::distance(ASV_physeq_stool, method = "bray")
sam <- data.frame(sample_data(ASV_physeq_stool))
abrel_bray <- as.matrix(abrel_bray)

# our permanova test for healthy vs. CRS
sam_CRS <- sam[!is.na(sam$CRS_healthy), ]
abrel_bray_CRS <- abrel_bray[rownames(abrel_bray) %in% sam_CRS$SampleID, colnames(abrel_bray) %in% sam_CRS$SampleID]
permanova <- adonis2(abrel_bray_CRS ~ CRS_healthy, data = sam_CRS)

# our permanova test for healthy vs. d28
sam_CRS_d28 <- sam[!is.na(sam$healtyh_d28), ]
abrel_bray_CRS_d28 <- abrel_bray[rownames(abrel_bray) %in% sam_CRS_d28$SampleID, colnames(abrel_bray) %in% sam_CRS_d28$SampleID]
permanova <- adonis2(abrel_bray_CRS_d28 ~ healtyh_d28, data = sam_CRS_d28)

# our permanova test for healthy vs. d90
sam_CRS_d90 <- sam[!is.na(sam$healthy_d90), ]
abrel_bray_CRS_d90 <- abrel_bray[rownames(abrel_bray) %in% sam_CRS_d90$SampleID, colnames(abrel_bray) %in% sam_CRS_d90$SampleID]
permanova <- adonis2(abrel_bray_CRS_d90 ~ healthy_d90, data = sam_CRS_d90)

# our permanova test for healthy vs. d180
sam_180 <- sam[!is.na(sam$healthy_d180), ]
abrel_bray_180 <- abrel_bray[rownames(abrel_bray) %in% sam_180$SampleID, colnames(abrel_bray) %in% sam_180$SampleID]
permanova2 <- adonis2(abrel_bray_180 ~ healthy_d180, data = sam_180)

# our permanova test for d0 vs. d28
sam_28 <- sam[!is.na(sam$d0_d28), ]
abrel_bray_28 <- abrel_bray[rownames(abrel_bray) %in% sam_28$SampleID, colnames(abrel_bray) %in% sam_28$SampleID]
permanova2 <- adonis2(abrel_bray_28 ~ d0_d28, data = sam_28)

# our permanova test for d0 vs. d90
sam_90 <- sam[!is.na(sam$d0_d90), ]
abrel_bray_90 <- abrel_bray[rownames(abrel_bray) %in% sam_90$SampleID, colnames(abrel_bray) %in% sam_90$SampleID]
permanova2 <- adonis2(abrel_bray_90 ~ d0_d90, data = sam_90)

# our permanova test for d0 vs. d180
sam_180 <- sam[!is.na(sam$d0_180), ]
abrel_bray_180 <- abrel_bray[rownames(abrel_bray) %in% sam_180$SampleID, colnames(abrel_bray) %in% sam_180$SampleID]
permanova2 <- adonis2(abrel_bray_180 ~ d0_180, data = sam_180)
```
The permanova test for healthy vs CRS is significant with a F value of 1.864 and p-value 0.011. This means the centroid of the two clusters are not the same and are significantly different from each other.

### nasal microbiome
| Test           | F-value | p-value |
|----------------|---------|---------|
|healthy vs. CRS |1.864    |0.015*  
|healthy vs. d28 |3.067    |0.001***
|healthy vs. d90 |2.941    |0.001***
|healthy vs. d180|2.342    |0.003**
|d0 vs. d28      |0.892    |0.572
|d0 vs. d90      |1.036    |0.373
|d0 vs. d180     |1.391    |0.105

### stool microbiome
| Test           | F-value | p-value |
|----------------|---------|---------|
|healthy vs. CRS |2.816    |0.001*** 
|healthy vs. d28 |2.372    |0.001***
|healthy vs. d90 |2.087    |0.002**
|healthy vs. d180|1.725    |0.015*
|d0 vs. d28      |0.877    |0.691
|d0 vs. d90      |0.864    |0.748
|d0 vs. d180     |0.848    |0.769
