
setwd('~/Documents/CU_project/Rotation2/TP01_meta/')
set.seed(1234)
library(Maaslin2)
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)


# 1. Metaphlan taxa check - see if all taxa go all the way down to species ####
# All goes down to species
input_data = read.csv('~/Documents/CU_project/Rotation2/TP01_meta/metaphlan/metaphlan_taxonomic_profiles.csv',check.names=FALSE, row.names = 'Taxonomy') %>% as.data.frame() 

taxa = rownames(input_data) %>% as.data.frame()
colnames(taxa) = 'clade_name'
taxa = separate_wider_delim(taxa, cols = clade_name, delim = "|", names = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), too_few = "align_start")

taxa = taxa %>% apply(2, function(column) word(column, 2, sep = "__")) %>% as.data.frame()

n = 1

for (h in colnames(taxa)[1:6]){
  print(h)
  
  for (i in unique(taxa[,h]) ){
    tmp = taxa[taxa[,h] == i,] %>% distinct
    last = apply(tmp, 2, function(col) tail(col[!is.na(col)], 1)) %>% unlist %>% names %>% tail(n = 1)
    
    taxa_low_tmp =  matrix(c(i,last), ncol = 2) %>% as.data.frame() ; colnames(taxa_low_tmp ) = c('ID','low')
    if (n == 1){
      taxa_low = taxa_low_tmp
    } else {
      taxa_low = rbind(taxa_low, taxa_low_tmp)
    }
    
    n = n+1
    
  }
  
  lowes = taxa_low$low %>% unique 
  print(lowes)
}


# 2. Input preparation ####

## Taxonomy profile data from Metaphlan - Row: sample / Col: Species ####
input_data = read.csv('~/Documents/CU_project/Rotation2/TP01_meta/metaphlan/metaphlan_taxonomic_profiles.csv',check.names=FALSE, row.names = 'Taxonomy') %>% t() %>% as.data.frame() 
rownames(input_data)  = rownames(input_data) %>% word(1,sep = '_') %>% str_remove_all('TP01')

## QC Metaphlan data - filter out those not bacteria ####
input_data= input_data[,colSums(input_data) > 0]
input_data= input_data[,colnames(input_data) %>% str_detect('k__Bacteria')] # 72 701


## Metadata - 59 samples + 13 of Diabetes ####
# metadata for 59 samples: metadata_collection_compliance_filtered_20220714_NAfilt_jy_241008.csv  
input_metadata_raw = read.csv('~/Documents/CU_project/Rotation2/TP01_meta/metadata_collection_compliance_filtered_20220714_NAfilt_jy_241008.csv', row.names = 'X')
input_metadata_raw  = input_metadata_raw %>% filter(timepoint == 'TP01') ; rownames(input_metadata_raw) = input_metadata_raw$participant_id; input_metadata_raw  = input_metadata_raw %>% dplyr::select(AMY1CNfinal,AMY1Groupfinal)


# metadata for 13 diabetes: All_AMY1CN_q_d_PCR.csv
input_metadata_dia = read.csv('~/Documents/CU_project/Rotation2/TP01_meta/All_AMY1CN_q_d_PCR.csv')
input_metadata_dia= input_metadata_dia[(input_metadata_dia$participant_id %>% startsWith('5')) & (input_metadata_dia$participant_id %in% rownames(input_data)),]
input_metadata_dia = input_metadata_dia %>% dplyr::select(participant_id,qdPCR_rounded)
rownames(input_metadata_dia) = input_metadata_dia$participant_id; input_metadata_dia$participant_id = NULL
colnames(input_metadata_dia) = 'AMY1CNfinal'
input_metadata_dia$AMY1Groupfinal = 'Notset'

# Merge two metadata
input_metadata = rbind(input_metadata_raw, input_metadata_dia) # dim:  72  2
input_metadata$AMY1Groupfinal[input_metadata$AMY1Groupfinal == 'Notset'& input_metadata$AMY1CNfinal >8 ] = 'High'
input_metadata$AMY1Groupfinal[input_metadata$AMY1Groupfinal == 'Notset'& input_metadata$AMY1CNfinal < 5 ] = 'Low'
input_metadata$AMY1Groupfinal[input_metadata$AMY1Groupfinal == 'Notset' ] = 'Medium'
saveRDS(input_metadata,'input_metadata.rds')

input_metadata = readRDS('input_metadata.rds')


# 3. Analysis #####

## 3-2. Using lowest level of metaphlan data: No association ####
input_data_low = input_data[,(colnames(input_data) %>% str_detect('s__')) ]
rowSums(input_data_low) # all 100
colSums(input_data_low) %>% summary 

maas_CN_contin = Maaslin2( input_data = input_data_low, input_metadata = input_metadata, output = "CN_contin_species", fixed_effects = c("AMY1CNfinal"),normalization  = 'NONE',min_prevalence=0.05, correction = 'BH' )

## 3-6. Using the lowest level of metaphlan data + LM for each + no MaAslin2 ####
### Lowest level data retrieval  ####
input_data_low = input_data[,(colnames(input_data) %>% str_detect('s__')) ]
rowSums(input_data_low) # all 100
colSums(input_data_low) %>% summary 

tmp = merge(input_data_low, input_metadata, by = "row.names") # dim: 72 441

### Run analysis ####
dir.create('CN_contin_LM')
res_coeff = as.data.frame(matrix(ncol = 4, nrow = 1)); colnames(res_coeff) = c("Estimate"  , "Std. Error", "t value"   , "Pr(>|t|)" )

for (i in colnames(input_data_low)){
  print(i)
  
  if (round(sum(tmp[,i]),2) != 0){
    res = glm( log(tmp[,i] + 1) ~ scale(tmp$AMY1CNfinal) ) %>% summary
    res_coeff_tmp = res$coefficients[2,]  %>% as.data.frame() %>% t %>% as.data.frame()
    rownames(res_coeff_tmp) = i %>% str_replace_all('[|]','.')
    res_coeff = rbind(res_coeff,res_coeff_tmp)
    
    
    if (res$coefficients[2,4] < 0.05){
      
      png(paste0('CN_contin_LM/signifi_scatter_',i,'.png'), res = 110, height = 30, width = 30,units = 'cm')
      plot(scale(tmp$AMY1CNfinal), tmp[,i],ylab = i, main=paste0('coeff: ', round(res$coefficients[2,1],2),'\n p-val: ', round(res$coefficients[2,4],2)))
      dev.off()
    }
  }
}


write.csv(res_coeff[-1,],'CN_contin_LM/res_coeff.csv')


## 3-7. Using lowest level of metaphlan data + categorical : No association ####
input_data_low = input_data[,(colnames(input_data) %>% str_detect('s__')) ]
input_metadata$AMY1Groupfinal[input_metadata$AMY1Groupfinal == 'Notset'& input_metadata$AMY1CNfinal >8 ] = 'High'
input_metadata$AMY1Groupfinal[input_metadata$AMY1Groupfinal == 'Notset'& input_metadata$AMY1CNfinal < 5 ] = 'Low'
input_metadata$AMY1Groupfinal[input_metadata$AMY1Groupfinal == 'Notset' ] = 'Medium'

idx = rownames(input_metadata)[input_metadata$AMY1Groupfinal != 'Medium']

maas_CN_contin = Maaslin2( input_data = input_data_low[idx,], input_metadata = input_metadata[idx,], output = "CN_cate_species", fixed_effects = c("AMY1Groupfinal"),normalization  = 'NONE',min_prevalence=0.05, correction = 'BH' )


# 4. Metaphlan heatmap ####
## Retrieve lowest level ####
input_data_low = input_data[,(colnames(input_data) %>% str_detect('s__')) ]
input_data_low = input_data_low[,colnames(input_data_low) %>% str_detect('Bacteria')] %>% t()

## Generate column annotation ####
input_metadata = readRDS('./input_metadata.rds')
colnames(input_data_low) = colnames(input_data_low) %>% word(1, sep = '_') %>% str_remove_all('TP01')
input_metadata = arrange(input_metadata, AMY1CNfinal)

idx = match(rownames(input_metadata),colnames(input_data_low))
colnames(input_data_low)[idx] == rownames(input_metadata) 
input_data_low = input_data_low[,idx]
colnames(input_data_low) == rownames(input_metadata) 
input_metadata %>% colnames()

colours <- list(
  "AMY1CNfinal" = colorRamp2(c(min(input_metadata$AMY1CNfinal, na.rm = TRUE), 
                               max(input_metadata$AMY1CNfinal, na.rm = TRUE)),
                             c('#FDDA0D', '#4B6F44')),  # Continuous scale
  "AMY1Groupfinal" = setNames(c('#034694', '#ADD8E6', 'orange'), 
                              input_metadata$AMY1Groupfinal %>% unique())
)

# Create column annotations
colAnn <- HeatmapAnnotation(
  df = input_metadata,
  which = 'col',
  col = colours,
  annotation_width = unit(c(1, 4), 'cm'),
  gap = unit(1, 'mm'),
  annotation_legend_param = list(
    AMY1CNfinal = list(direction = "horizontal", show_annotation_name = FALSE),
    AMY1Groupfinal = list(direction = "horizontal", show_annotation_name = FALSE)
  )
)


## With relative abundace percentage ####
pdf('241205_heatmap_metaphlan_relab.pdf', width = 10, height = 10)
ComplexHeatmap::Heatmap(input_data_low, name="EC Relative Abundance (%)",show_row_names = F, cluster_columns = F,  use_raster = F,column_title = '', top_annotation = colAnn, column_names_gp =  grid::gpar(fontsize = 7))
dev.off()

col_fun = colorRamp2(c(0, mean(input_data_low),max(input_data_low)), c( "blue4","lightgray", "red3"))
pdf('241205_heatmap_metaphlan_relab_max_mean.pdf', width = 10, height = 10)
ComplexHeatmap::Heatmap(input_data_low, name="EC Relative Abundance (%)",show_row_names = F, cluster_columns = F,  use_raster = F,column_title = '', top_annotation = colAnn, column_names_gp =  grid::gpar(fontsize = 7), col = col_fun)
dev.off()


