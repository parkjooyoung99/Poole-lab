
setwd('~/Documents/CU_project/Rotation2/TP01_meta/humann/merged/')
set.seed(1234)
library(Maaslin2)
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(rstatix)
library(Seurat)

# 1. Subset table to the EC_relab lowest level ####
input_data = read.delim('ecs_relab.tsv',check.names=FALSE,row.names = '# Gene Family') %>% as.data.frame() 
input_data_pct = input_data *100

input_large = input_data_pct[! rownames(input_data_pct) %>% str_detect('[|]'),]
colSums(input_large) # 100
input_low =  input_data_pct[ rownames(input_data_pct) %>% str_detect('[|]'),]
colSums(input_low) # 100


# 2. EC Relab distribution check ####
df = as.data.frame(matrix(ncol = 2)); colnames(df) = c('pct','ID')

for (i in colnames(input_low_pct)){
  tmp = as.data.frame(as.matrix(input_low_pct[,i])); tmp$ID = i
  colnames(tmp) =  c('pct','ID')
  df = rbind(df, tmp)
}

df = df[-1,]
df$log = log(df$pct+1)
ggplot(df, aes(x = log, col = ID) ) + geom_histogram() + theme(legend.position = "none")  +
  xlim(0.5, 0.85)
ggplot(df, aes(x = log, col = ID) ) + geom_histogram() + theme(legend.position = "none")  +
  xlim(0, 0.5)

# 3. EC Relab heatmap ####
## Generate column annotation ####
input_metadata = readRDS('./input_metadata.rds')
colnames(input_large) = colnames(input_large) %>% word(1, sep = '_') %>% str_remove_all('TP01')
input_metadata = arrange(input_metadata, AMY1CNfinal)

idx = match(rownames(input_metadata),colnames(input_large))
colnames(input_large)[idx] == rownames(input_metadata) 
input_large = input_large[,idx]
colnames(input_large) == rownames(input_metadata) 
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
pdf('241205_EC_heatmap_relab.pdf', width = 10, height = 10)
ComplexHeatmap::Heatmap(input_large, name="EC Relative Abundance (%)",show_row_names = F, cluster_columns = F,  use_raster = F,column_title = '', top_annotation = colAnn, column_names_gp =  grid::gpar(fontsize = 7))
dev.off()

col_fun = colorRamp2(c(0, 0.04,max(input_large)), c( "blue4","lightgray", "red3")) # 0.04: mean for every samples
pdf('241205_EC_heatmap_relab_max_mean.pdf', width = 10, height = 10)
ComplexHeatmap::Heatmap(input_large, name="EC Relative Abundance (%)",show_row_names = F, cluster_columns = F,  use_raster = F,column_title = '', top_annotation = colAnn, column_names_gp =  grid::gpar(fontsize = 7), col = col_fun)
dev.off()

## With log transformed ####
pdf('241205_EC_heatmap_relab_log.pdf', width = 10, height = 10)
ComplexHeatmap::Heatmap(log(input_large+1), name="EC Relative Abundance (%)",show_row_names = F, cluster_columns = F,  use_raster = F,column_title = '', top_annotation = colAnn, column_names_gp =  grid::gpar(fontsize = 7))
dev.off()

col_fun = colorRamp2(c(0, 0.04,max(log(input_large+1))), c( "blue4","lightgray", "red3")) # 0.04: mean for every samples
pdf('241205_EC_heatmap_relab_log_max_mean.pdf', width = 10, height = 10)
ComplexHeatmap::Heatmap(log(input_large+1), name="EC Relative Abundance (%)",show_row_names = F, cluster_columns = F,  use_raster = F,column_title = '', top_annotation = colAnn, column_names_gp =  grid::gpar(fontsize = 7), col = col_fun)
dev.off()


# 4. Analysis #####
## 4-1. Maaslin2 with EC relab ####
dir.create('EC_CN_contin')
colnames(input_large) = colnames(input_large) %>% word(1, sep = '_') %>% str_remove_all('TP01')
input_metadata = readRDS('../../input_metadata.rds')
input_metadata = input_metadata[colnames(input_large),]

maas_CN_contin = Maaslin2( input_data = t(input_large), input_metadata = input_metadata, output = "EC_CN_contin", fixed_effects = c("AMY1CNfinal"),normalization  = 'NONE',min_prevalence=0.05, correction = 'BH',transform = "LOG")

## 4-2. Spearman correlation with EC relab ####
dir.create('EC_CN_contin_spear')

input_large = t(input_large)
tmp = merge(input_large, input_metadata, by = "row.names") # dim: 72 441

res_coeff = as.data.frame(matrix(ncol = 2, nrow = 1)); colnames(res_coeff) = c("rho"   , "p-value" )

for (i in colnames(input_large)){
  print(i)
  
  if (round(sum(tmp[,i]),2) != 0){
    res = cor.test(  tmp[,i] , scale(tmp$AMY1CNfinal) , method = 'spearman')
    res_coeff_tmp = c(res$estimate, res$p.value) %>% as.data.frame() %>% t(); colnames(res_coeff_tmp) = c("rho"   , "p-value" ); rownames(res_coeff_tmp) = i
    res_coeff = rbind(res_coeff,res_coeff_tmp)
    
    
    if (res$p.value < 0.05){
      
      png(paste0('EC_CN_contin_spear/signifi_scatter_',i,'.png'), res = 110, height = 30, width = 30,units = 'cm')
      plot(scale(tmp$AMY1CNfinal), tmp[,i],ylab = i, main=paste0('rho: ', round(res$estimate,2),'\n p-val: ', round(res$p.value,2)))
      dev.off()
    }
  }
}


write.csv(res_coeff[-1,],'EC_CN_contin_spear/res_rho.csv')


### Draw scatter plot for CAZymes ####
input_metadata = readRDS('~/Documents/CU_project/Rotation2/TP01_meta/input_metadata.rds')
input_large = t(input_large)
tmp = merge(input_large, input_metadata, by = "row.names") # dim: 72 441
tmp$AMY1CNfinal_scale = tmp$AMY1CNfinal %>% scale

rho = read.csv('./EC_CN_contin_spear/res_rho.csv')
sig = c('2.4.1.291','2.4.1.305','3.2.1.39','2.4.1.290','2.4.1.303','4.2.1.3')

for (i in sig){
  print(i)
  
  pval = rho[rho$X == i,"p.value"] %>% round(2)
  rhoo = rho[rho$X == i,"rho"] %>% round(2)
  
  pdf(paste0('EC_CN_contin_spear/spear_',i,'.pdf'), width = 4.5*1.2,height = 4.1*1.2)
  print(ggplot(tmp, aes(x = AMY1CNfinal_scale , y = tmp[,i] ,color = AMY1Groupfinal)) + geom_point() + scale_color_manual( values = c( 'orange','#034694','#ADD8E6'))  + theme_classic() + NoLegend()  + ylab(i) + xlab('Scaled AMY1 CN')  )
  dev.off()
}

pdf('humann/merged/EC_CN_contin_spear/spear_2.4.1.291.pdf',width = 4.5*1.2,height = 5*1.2)
ggplot(tmp, aes(x = AMY1CNfinal_scale , y = tmp[,'2.4.1.291'] ,color = AMY1Groupfinal)) + geom_point() + scale_color_manual( values = c( 'orange','#034694','#ADD8E6'))  + theme_classic() + NoLegend()  + ylab('2.4.1.291') + xlab('Scaled AMY1 CN')
dev.off()

## 4-3. Using the lowest level of metaphlan data + LM for each + no MaAslin2 ####
### Lowest level data retrieval  ####
tmp = merge(t(input_large), input_metadata, by = "row.names") # dim: 72 2287

### Run analysis ####
dir.create('EC_CN_contin_LM')
res_coeff = as.data.frame(matrix(ncol = 4, nrow = 1)); colnames(res_coeff) = c("Estimate"  , "Std. Error", "t value"   , "Pr(>|t|)" )

for (i in rownames(input_large)){
  print(i)
  
  if (round(sum(tmp[,i]),2) != 0){
    res = glm( log(tmp[,i] + 1) ~ scale(tmp$AMY1CNfinal) ) %>% summary
    res_coeff_tmp = res$coefficients[2,]  %>% as.data.frame() %>% t %>% as.data.frame()
    rownames(res_coeff_tmp) = i %>% str_replace_all('[|]','.')
    res_coeff = rbind(res_coeff,res_coeff_tmp)
    
    
    if (res$coefficients[2,4] < 0.05){
      
      png(paste0('EC_CN_contin_LM/signifi_scatter_',i,'.png'), res = 110, height = 30, width = 30,units = 'cm')
      plot(scale(tmp$AMY1CNfinal), tmp[,i],ylab = i, main=paste0('coeff: ', round(res$coefficients[2,1],2),'\n p-val: ', round(res$coefficients[2,4],2)))
      dev.off()
    }
  }
}


write.csv(res_coeff[-1,],'EC_CN_contin_LM/res_coeff.csv')


## 4-4. Maaslin2 with EC relab - cate ####
dir.create('EC_CN_wilcox_cate')
colnames(input_large) = colnames(input_large) %>% word(1, sep = '_') %>% str_remove_all('TP01')
input_metadata = readRDS('../../input_metadata.rds')
input_metadata_sub = input_metadata %>% filter(AMY1Groupfinal != 'Medium')
input_large_sub = input_large[,rownames(input_metadata_sub)]


maas_CN_contin = Maaslin2( input_data = t(input_large_sub), input_metadata = input_metadata_sub, output = "EC_CN_cate", fixed_effects = c("AMY1Groupfinal"),normalization  = 'NONE',min_prevalence=0.05, correction = 'BH',transform = "LOG")

## 4-5. Wilcox test EC relab - cate ####
dir.create('EC_CN_wilcox_cate')
tmp = merge(t(input_large_sub), input_metadata_sub, by = "row.names") # dim:  49 2287
res_coeff = as.data.frame(matrix(ncol = 3, nrow = 1)); colnames(res_coeff) = c("Med_high"  , "Med_low", "p-val" )

for (i in colnames(tmp)){
  if (! i %in% c('AMY1CNfinal','AMY1Groupfinal','Row.names')){
    if (sum(tmp[,i])){
      formula <- as.formula(paste0("`", i, "` ~ AMY1Groupfinal"))
      res =  wilcox_test(formula = formula, data = tmp, p.adjust.method = 'bonferroni')
      
      if (res$p < 0.05){
        stat.test <- res %>%
          add_significance() %>%
          add_xy_position(x="condition",step.increase=1)
        stat.test$manual_position <-max(tmp[,i] + 0.01)
        stat.test$label <- stat.test$p.signif
        
        pdf(paste0('EC_CN_wilcox_cate/wilcox_sig_',i,'.pdf'),width = 4.5*1.2,height = 5*1.2)
        print(ggplot(tmp, aes(x = AMY1Groupfinal , y = tmp[,i])) + geom_violin(aes(fill = AMY1Groupfinal))+geom_boxplot( aes(fill = AMY1Groupfinal),alpha = 0.3, outlier.shape = NA)   + theme_classic() + NoLegend()  + ylab(i) + scale_fill_manual( values = c( 'orange','#034694')) + ggsignif::geom_signif(data=stat.test,aes(xmin=group1,xmax=group2,annotations=label,y_position=manual_position),manual=TRUE, inherit.aes=FALSE))
        dev.off()
        
        res_coeff_tmp = as.data.frame(matrix(ncol = 3, nrow = 1)); colnames(res_coeff_tmp) = c("Med_high"  , "Med_low", "p-val" )
        res_coeff_tmp[,3] = res$p
        idx = which(input_metadata_sub$AMY1Groupfinal != 'High')
        res_coeff_tmp[,2] = median(tmp[idx,i]) %>% round(5)
        idx = which(input_metadata_sub$AMY1Groupfinal == 'High')
        res_coeff_tmp[,1] = median(tmp[idx,i]) %>% round(5)
        rownames(res_coeff_tmp) = i
        
        res_coeff = rbind(res_coeff, res_coeff_tmp)
      }
    }
  }
}

write.csv(res_coeff[-1,],'EC_CN_wilcox_cate/res_wilcox.csv')

res_coeff = read.csv('EC_CN_wilcox_cate/res_wilcox.csv')
write.csv(res_coeff[res_coeff$Med_high < res_coeff$Med_low,],'EC_CN_wilcox_cate/res_low_wilcox.csv')
write.csv(res_coeff[res_coeff$Med_high > res_coeff$Med_low,],'EC_CN_wilcox_cate/res_high_wilcox.csv')

