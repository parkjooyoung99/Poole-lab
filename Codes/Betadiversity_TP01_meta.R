
setwd('~/Documents/CU_project/Rotation2/TP01_meta/')
set.seed(1234)
library(dplyr)
library(stringr)
library(ggplot2)
library(phyloseq)
library(DESeq2)
library(ggrepel)

# 1. With taxonomy - Metaphlan ####
## 1. Generate phyloseq data from metaphlan ####
### Taxonomy profile from Metaphlan ####
mphlanin = read.csv('~/Documents/CU_project/Rotation2/TP01_meta/metaphlan/metaphlan_taxonomic_profiles.csv',check.names=FALSE, row.names = 'Taxonomy') 
colnames(mphlanin)  = colnames(mphlanin) %>% word(1,sep = '_') %>% str_remove_all('TP01')
mphlanin= mphlanin[rowSums(mphlanin) > 0, ]
mphlanin= mphlanin[rownames(mphlanin) %>% str_detect('k__Bacteria'),] 
mphlanin = mphlanin[(rownames(mphlanin) %>% str_detect('s__')), ] ## 434  72

### Metatdata ####
input_metadata_raw = read.csv('~/Documents/CU_project/Rotation2/TP01_meta/metadata_collection_compliance_filtered_20220714_NAfilt_jy_241008.csv', row.names = 'X')
input_metadata_raw  = input_metadata_raw %>% filter(timepoint == 'TP01') ; rownames(input_metadata_raw) = input_metadata_raw$participant_id; input_metadata_raw  = input_metadata_raw %>% dplyr::select(AMY1CNfinal,AMY1Groupfinal)
input_metadata_dia = read.csv('~/Documents/CU_project/Rotation2/TP01_meta/All_AMY1CN_q_d_PCR.csv')
input_metadata_dia= input_metadata_dia[(input_metadata_dia$participant_id %>% startsWith('5')) & (input_metadata_dia$participant_id %in% colnames(mphlanin)),]
input_metadata_dia = input_metadata_dia %>% dplyr::select(participant_id,qdPCR_rounded)
rownames(input_metadata_dia) = input_metadata_dia$participant_id; input_metadata_dia$participant_id = NULL
colnames(input_metadata_dia) = 'AMY1CNfinal'
input_metadata_dia$AMY1Groupfinal = 'Notset'
input_metadata = rbind(input_metadata_raw, input_metadata_dia) ## dim:  72  2


### Phyloseq object  ####
## source code from https://www.biostars.org/p/449688/
metaphlanToPhyloseq <- function(
    tax,
    metadat=NULL,
    simplenames=TRUE,
    roundtointeger=FALSE,
    split="|"){
  ### tax is a matrix or data.frame with the table of taxonomic abundances, rows are taxa, columns are samples
  ### metadat is an optional data.frame of specimen metadata, rows are samples, columns are variables
  ### if simplenames=TRUE, use only the most detailed level of taxa names in the final object
  ### if roundtointeger=TRUE, values will be rounded to the nearest integer
  xnames = rownames(tax)
  shortnames = gsub(paste0(".+\\", split), "", xnames)
  if(simplenames){
    rownames(tax) = shortnames
  }
  if(roundtointeger){
    tax = round(tax * 1e4)
  }
  x2 = strsplit(xnames, split=split, fixed=TRUE)
  taxmat = matrix(NA, ncol=max(sapply(x2, length)), nrow=length(x2))
  colnames(taxmat) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")[1:ncol(taxmat)]
  rownames(taxmat) = rownames(tax)
  for (i in 1:nrow(taxmat)){
    taxmat[i, 1:length(x2[[i]])] <- x2[[i]]
  }
  taxmat = gsub("[a-z]__", "", taxmat)
  taxmat = phyloseq::tax_table(taxmat)
  otutab = phyloseq::otu_table(tax, taxa_are_rows=TRUE)
  if(is.null(metadat)){
    res = phyloseq::phyloseq(taxmat, otutab)
  }else{
    res = phyloseq::phyloseq(taxmat, otutab, phyloseq::sample_data(metadat))
  }
  return(res)
}

phyloseqin= metaphlanToPhyloseq(mphlanin, metadat = input_metadata)

## 2. Explore phyloseq object ####
### Taxa information ####
phyloseqin@tax_table

### Relative abundance - remember this is metagenomic data so it is ra not otu ####
phyloseqin@otu_table


## 3. Calculate Beta diversity + Generate Plot ####
meta_ord <- ordinate(physeq = phyloseqin, method = "PCoA",  distance = "bray")
plot_ordination(physeq = phyloseqin, ordination = meta_ord)

coordinates_df = meta_ord$vectors %>% as.data.frame()  ## 72 49
colnames(coordinates_df) = paste0(colnames(coordinates_df) %>% str_replace_all('Axis[.]','PC_'),' [',round(meta_ord$values[meta_ord$values[,2] > 0 ,2] * 100,2),'%]')
coordinates_df = cbind(coordinates_df, input_metadata)
coordinates_df$ID = rownames(coordinates_df)


coordinates_df$AMY1Groupfinal[coordinates_df$AMY1Groupfinal == 'Notset'& coordinates_df$AMY1CNfinal >8 ] = 'High'
coordinates_df$AMY1Groupfinal[coordinates_df$AMY1Groupfinal == 'Notset'& coordinates_df$AMY1CNfinal < 5 ] = 'Low'
coordinates_df$AMY1Groupfinal[coordinates_df$AMY1Groupfinal == 'Notset' ] = 'Medium'

pdf('./241204_beta_categorical.pdf', width = 5, height = 4)
ggplot(coordinates_df, aes(x = `PC_1 [14.04%]`, y = `PC_2 [10.45%]`, color = AMY1Groupfinal, label = ID )) + geom_point(size = 1)  + theme_classic()+ geom_text_repel(size = 2) + scale_colour_manual(values = c('orange','#034694', '#ADD8E6' )) + ggtitle('PcoA plot for Beta Diversity') + theme(plot.title = element_text(hjust = 0.5, size = 10), legend.title = element_text(size = 8), legend.text = element_text(size = 7), axis.title = element_text(size = 9), axis.text = element_text(size = 6))
dev.off()

pdf('./241204_beta_continu.pdf', width = 4.9, height = 4)
ggplot(coordinates_df, aes(x = `PC_1 [14.04%]`, y = `PC_2 [10.45%]`, color = AMY1CNfinal, label = ID)) + geom_point(size = 1) + scale_color_gradient(low = 'lightgray',high = 'red3') + theme_classic() + geom_text_repel(size = 2) + theme(plot.title = element_text(hjust = 0.5, size = 10), legend.title = element_text(size = 8), legend.text = element_text(size = 7), axis.title = element_text(size = 9), axis.text = element_text(size = 6))  + ggtitle('PcoA plot for Beta Diversity') 
dev.off()


# 2. With EC - Humann3 ####
## 2-1. Generate phyloseq data from metaphlan ####
### EC profile from Metaphlan ####
mphlanin = read.delim('./humann/merged/ecs_relab.tsv',check.names=FALSE,row.names = '# Gene Family') %>% as.data.frame() 
mphlanin = mphlanin *100
mphlanin = mphlanin[! rownames(mphlanin) %>% str_detect('[|]'),]
colnames(mphlanin)  = colnames(mphlanin) %>% word(1,sep = '_') %>% str_remove_all('TP01')
mphlanin= mphlanin[rowSums(mphlanin) > 0, ]

### Metatdata ####
input_metadata = readRDS('input_metadata.rds')

### Phyloseq object  ####
## source code from https://www.biostars.org/p/449688/
metaphlanToPhyloseq <- function(
    tax,
    metadat=NULL,
    simplenames=TRUE,
    roundtointeger=FALSE,
    split="|"){
  ### tax is a matrix or data.frame with the table of taxonomic abundances, rows are taxa, columns are samples
  ### metadat is an optional data.frame of specimen metadata, rows are samples, columns are variables
  ### if simplenames=TRUE, use only the most detailed level of taxa names in the final object
  ### if roundtointeger=TRUE, values will be rounded to the nearest integer
  xnames = rownames(tax)
  shortnames = gsub(paste0(".+\\", split), "", xnames)
  if(simplenames){
    rownames(tax) = shortnames
  }
  if(roundtointeger){
    tax = round(tax * 1e4)
  }
  x2 = strsplit(xnames, split=split, fixed=TRUE)
  taxmat = matrix(NA, ncol=max(sapply(x2, length)), nrow=length(x2))
  colnames(taxmat) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")[1:ncol(taxmat)]
  rownames(taxmat) = rownames(tax)
  for (i in 1:nrow(taxmat)){
    taxmat[i, 1:length(x2[[i]])] <- x2[[i]]
  }
  taxmat = gsub("[a-z]__", "", taxmat)
  taxmat = phyloseq::tax_table(taxmat)
  otutab = phyloseq::otu_table(tax, taxa_are_rows=TRUE)
  if(is.null(metadat)){
    res = phyloseq::phyloseq(taxmat, otutab)
  }else{
    res = phyloseq::phyloseq(taxmat, otutab, phyloseq::sample_data(metadat))
  }
  return(res)
}

phyloseqin= metaphlanToPhyloseq(mphlanin, metadat = input_metadata)

## 3. Calculate Beta diversity + Generate Plot ####
meta_ord <- ordinate(physeq = phyloseqin, method = "PCoA",  distance = "bray")
plot_ordination(physeq = phyloseqin, ordination = meta_ord)

coordinates_df = meta_ord$vectors %>% as.data.frame()  ## 72 49
colnames(coordinates_df) = paste0(colnames(coordinates_df) %>% str_replace_all('Axis[.]','PC_'),' [',round(meta_ord$values[meta_ord$values[,2] > 0 ,2] * 100,2),'%]')
coordinates_df = cbind(coordinates_df, input_metadata)
coordinates_df$ID = rownames(coordinates_df)


pdf('./241204_ec_beta_categorical.pdf', width = 5, height = 4)
ggplot(coordinates_df, aes(x = `PC_1 [24.73%]`, y = `PC_2 [21.46%]`, color = AMY1Groupfinal, label = ID )) + geom_point(size = 1)  + theme_classic()+ geom_text_repel(size = 2) + scale_colour_manual(values = c('orange','#034694', '#ADD8E6' )) + ggtitle('PcoA plot for Beta Diversity') + theme(plot.title = element_text(hjust = 0.5, size = 10), legend.title = element_text(size = 8), legend.text = element_text(size = 7), axis.title = element_text(size = 9), axis.text = element_text(size = 6))
dev.off()

pdf('./241204_ec_beta_continu.pdf', width = 4.9, height = 4)
ggplot(coordinates_df, aes(x = `PC_1 [24.73%]`, y = `PC_2 [21.46%]`, color = AMY1CNfinal, label = ID)) + geom_point(size = 1) + scale_color_gradient(low = 'lightgray',high = 'red3') + theme_classic() + geom_text_repel(size = 2) + theme(plot.title = element_text(hjust = 0.5, size = 10), legend.title = element_text(size = 8), legend.text = element_text(size = 7), axis.title = element_text(size = 9), axis.text = element_text(size = 6))  + ggtitle('PcoA plot for Beta Diversity') 
dev.off()
