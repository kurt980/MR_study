---
title: "0912_heritability_analysis"
author: "Chengyan Ji"
date: "2024-09-12"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

**This is for analyzing GCTA heritability results with updated instructions**
**09/12/2024: Started editing**
**09/13/2024: Added brain region plots (ggseg)**
**09/15/2024: 鉴于Mac上的很多Package安装起来很别扭，以下直接在我的windows跑的**
**09/26/2024: 当前画的是p-value筛选之后的**
**09/26/2024: 鉴于要做Presentation，现在用的是3月份的heritability画的图**
**10/30/2024: 修改了一下画图的代码，改了一下颜色**
**10/30/2024: 了解了一下改变大脑位置的代码，改了一下legend的代码**
**10/31/2024: 画boxplot在network(转移到1010_network_study了)**


## Read GREML Results

```{r}
library(dplyr)## Plot ROIs on brain regions
library(snpStats)
library(tibble)
library(dplyr)
library(AER)
library(readxl)
library(ggplot2)
library(ggseg)
library(ggsegGlasser)

### Load Utilities, data, tool and define directory
utilities_path <- "/Users/grc8mu/Desktop/DS/ABCD_data/ABCD_SCFC_utilities/" # "/Users/grc8mu/Desktop/DS/ABCD_data/ABCD_SCFC_utilities/"
bfile_path <- "/Users/grc8mu/Desktop/DS/ABCD_data/ABCD_bfile/"
cov_path <- "/Users/grc8mu/Desktop/DS/ABCD_data/Covariates/"
harm_path <- "/Users/grc8mu/Desktop/DS/ABCD_data/Harmonized_scfc/"
pheno_path <- "/Users/grc8mu/Desktop/DS/ABCD_data/Phenotype/"
tool_path <- "/Users/grc8mu/Desktop/DS/Tools/"

directory <- "" #"/Users/grc8mu/Desktop/DS/1010_NetworkStudy/"
load(paste0(utilities_path, "coupling_cbcl.RData"))
load(paste0(utilities_path, "scfc_coupling_matched.RData"))

h_2_values <- numeric(379)
p_values <- numeric(379)
SEs <- numeric(379)
Ns <- numeric(379)

base_path <- "" #"/Users/kurt/Documents/Data_Science_Project/0912_heritability/"
harms = c(2)

H2_harm_df <- read.csv(paste0(base_path, "greml_results/H2_harm_sd_df.csv"), header = TRUE) %>%
  select(roi_num, ROI, h2_harm2, p_value_harm2, SE_harm2, N_harm2) %>%
  rename(h2 = h2_harm2, p_value = p_value_harm2, SE = SE_harm2, N = N_harm2)

# Adjust p-values using BH and Bonferroni methods
H2_harm_df <- H2_harm_df %>%
  mutate(
    p_value_bh = p.adjust(p_value, method = "BH"),  # BH-adjusted p-values
    p_value_bonferroni = p.adjust(p_value, method = "bonferroni"),  # Bonferroni-adjusted p-values
    significant_raw = ifelse(p_value < 0.05, TRUE, FALSE),  # Significant based on raw p-values
    significant_BH = ifelse(p_value_bh < 0.05, TRUE, FALSE),  # Significant based on BH-adjusted p-values
    significant_bonferroni = ifelse(p_value_bonferroni < 0.05, TRUE, FALSE)  # Significant based on Bonferroni-adjusted p-values
  )

sum(H2_harm_df$significant_BH)

write.csv(H2_harm_df, "H2_pvalue_harm2.csv")
```

```{r}
for (harm in harms) {

  h_2_p_value <- data.frame(ROI = 1:nrow(H2_harm_df), h_2 = H2_harm_df[[paste0("h2_harm", harm)]], p_value = H2_harm_df[[paste0("p_value_harm", harm)]], SE = H2_harm_df[[paste0("SE_harm", harm)]], N = H2_harm_df[[paste0("N_harm", harm)]])
    
  #---------------------------------------plot distribution of H2 values (all)-------------------------------------
  library(ggplot2)
  p1 <- ggplot(h_2_p_value, aes(x = h_2)) +
    geom_histogram(bins = 50, fill = "skyblue", color = "black") +
    labs(title = "Distribution of h2 Values", x = "h2 Values", y = "Frequency") +
    theme_minimal()
  
  max_value <- max(h_2_p_value$h_2)
  mean_value <- mean(h_2_p_value$h_2)
  total_observations <- length(h_2_p_value$h_2)
  
  p1 <- p1 + annotate("text", x = Inf, y = Inf, label = sprintf("Max: %.3f\nMean: %.3f\nN: %d", max_value, mean_value, total_observations),
                    hjust = 1.1, vjust = 1.1, size = 5)
  print(p1)
  ggsave(paste0(base_path, "analyze_results/herit_harm", harm, "_H2_values_Hist.png"), plot=p1)
  
  #---------------------------------Plot Distribution of H2 values (Filtered)------------------------------------
  criterion <- 0.05 / 379
  
  filtered_df <- h_2_p_value[h_2_p_value$p_value < criterion, ]
  p2 <- ggplot(filtered_df, aes(x = h_2)) +
    geom_histogram(bins = 50, fill = "orange", color = "black") +
    labs(title = "Distribution of h_2 Values for only significant ROIs",
         x = "h_2 Value",
         y = "Frequency") +
    theme_minimal()
  
  max_value <- max(filtered_df$h_2)
  mean_value <- mean(filtered_df$h_2)
  total_observations <- length(filtered_df$h_2)
  
  p2 <- p2 + annotate("text", x = Inf, y = Inf, label = sprintf("Max: %.3f\nMean: %.3f\nN: %d", max_value, mean_value, total_observations),
                    hjust = 1.1, vjust = 1.1, size = 5)
  
  print(p2)
  ggsave(paste0(base_path, "analyze_results/herit_harm", harm, "_H2_values_ROIs_sig.png"), plot=p2)
  
  set.seed(123) # For reproducibility
  #-------------------------------------------Select the 5 rows with the smallest p-values---------------------------------------------------
  sampled_df <- filtered_df[order(filtered_df$p_value)[1:5], ]
  p3 <- ggplot(sampled_df, aes(x = as.factor(ROI), y = h_2)) +
    geom_point(aes(color = "h_2 Value"), size = 4) + # Plot h_2 values as points
    geom_errorbar(aes(ymin = h_2 - SE, ymax = h_2 + SE, color = "Standard Error"), width = 0.2) + # Add error bars for SE and use color to match legend
    scale_color_manual(name = "", values = c("h_2 Value" = "blue", "Standard Error" = "red")) + # Define custom colors for points and error bars
    labs(title = "ROIs with 5 Smallest P-values",
         x = "ROI",
         y = "h_2 Value") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    guides(color = guide_legend(override.aes = list(shape = 16, linetype = "solid"))) # Ensure the legend accurately represents both points and error bars
  
  print(p3)
  ggsave(paste0(base_path, "analyze_results/herit_harm", harm, "_Top5ROI_by_pvalue.png"), plot=p3)
  
  #---------------------------------------------------------5 Smallest P-values-----------------------------------------------------------
  # Select the 5 rows with the smallest p-values
  sampled_df <- filtered_df[order(filtered_df$h_2, decreasing = TRUE)[1:5], ]
  p4 <- ggplot(sampled_df, aes(x = as.factor(ROI), y = h_2)) +
    geom_point(aes(color = "h_2 Value"), size = 4) + # Plot h_2 values as points
    geom_errorbar(aes(ymin = h_2 - SE, ymax = h_2 + SE, color = "Standard Error"), width = 0.2) + # Add error bars for SE and use color to match legend
    scale_color_manual(name = "", values = c("h_2 Value" = "blue", "Standard Error" = "red")) + # Define custom colors for points and error bars
    labs(title = "ROIs with 5 Biggest H2 values",
         x = "ROI",
         y = "h_2 Value") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    guides(color = guide_legend(override.aes = list(shape = 16, linetype = "solid"))) # Ensure the legend accurately represents both points and error bars
  
  print(p4)
  ggsave(paste0(base_path, "analyze_results/herit_harm", harm, "_Top5ROI_by_h2.png"), plot=p4)
}
```


## Plot ROIs on brain regions
```{r}
library(RColorBrewer)
library(scales)
display.brewer.all()
color <- rev(brewer.pal(9, 'YlGn'))
show_col(color)
cl <- c("#1A1A1A","#252525","#67000D","#49006A","#542788","#006837","#41AB5D","#FED976","#FEB24C","#FD8D3C","#FC4E2A","#E31A1C")
show_col(cl)
# cl <- c("white","#73EEF3","#37a5f1","#FED976","#FEB24C", "#e037f1", "#FC4E2A", "#ff0004")
# cl <- c("#000000", "#000033", "#00008B", "#1E3A8A", "#1E90FF", "#87CEFA", "#B0E0E6", 
#             "#FF6347", "#FF8C00", "#FFA500", "#FFD700", "#fbff00")
cl <- c("#000000", "#00008B", "#1E90FF", "#FC4E2A", "#E31A1C", "#920707")
cl <- c("#ffff76", "#FEB24C", "#FC4E2A", "#E31A1C", "#920707")
show_col(cl)
display.brewer.pal(11,'RdGy')
```

```{r}
base_path <- ""
H2_harm_df_all <- read.csv(paste0(base_path, "greml_results/H2_harm_sd_df.csv"), header = T)
# names(H2_harm_df)[1] = 'ROI'
names(H2_harm_df_all)[2] = 'name'

########## Brain Visualization #########
library(ggseg)
library(ggplot2)
library(dplyr)
### glasser atlas ########
### The first 360 regions of interest (ROIs) are from the Glasser atlas
library(ggsegGlasser)

#############################################
## Visualize the sc-fc coupling values ######
#############################################
library(stringr)

for (harm in harms) {
  criterion <- 0.05
  H2_harm_df <- H2_harm_df_all %>% mutate(bh_pval_harm = p.adjust(pull(., paste0("p_value_harm", harm)), method = "BH")) %>%
    filter(!!sym(paste0("p_value_harm", harm)) < criterion)
    #filter(bh_pval_harm < criterion)
  
  mmp <- mmp_subcor[1:360] %>% str_remove(.,'_ROI') # mmp_subcor
roi_glasser <- mmp %>% str_split(.,"_",simplify = TRUE) %>% data.frame()
colnames(roi_glasser) <- c("hemi","region")
roi_glasser <- roi_glasser %>% 
  tibble::rownames_to_column(var = "ROI") %>%
  mutate(ROI = as.integer(ROI)) %>%
  #mutate(ROI = paste0("roi",ROI)) %>%
  mutate(label=ifelse(hemi=='L',paste0("lh_L_",region), paste0("rh_R_",region)))

  scfc_harm_df <- data.frame(ROI = H2_harm_df[["roi_num"]], h_2 = H2_harm_df[[paste0("h2_harm", harm)]], p_value = H2_harm_df[[paste0("p_value_harm", harm)]], SE = H2_harm_df[[paste0("SE_harm", harm)]], N = H2_harm_df[[paste0("N_harm", harm)]])
  
  #===========================================Plot brain regions for harm 0 1 and 2========================================
  #-------------------------------------------plot first 360 ROIs----------------------------------------------------------
  roi_glasser <- roi_glasser %>%
    inner_join(scfc_harm_df, by = "ROI")
  
  gldata <- tibble(
    label = roi_glasser$label,
    Estimate = roi_glasser$h_2,
  )
  ### cortical_position (i.e., how do we want the brain to be positioned)
  cortical_pos <- c("left lateral", "right lateral","left medial","right medial")
  ##### plot the coupling value in brain
  min_value <- min(gldata$Estimate, na.rm = TRUE)
  quantiles <- quantile(gldata$Estimate, probs = c(0.5), na.rm = TRUE)
  max_value <- max(gldata$Estimate, na.rm = TRUE)
  plot <- gldata %>%
    ggplot() +
    geom_brain(
      atlas = glasser, 
      position = position_brain(cortical_pos), # side ~ hemi
      # aes(fill = Estimate)) + scale_fill_gradientn(colours = cl, na.value = "#e8e8e8",
      #                                               breaks = c(0),
      #                                               labels = c(paste0("Min: ", round(min_value, 2)), 
      #                                                          paste0("Median: ", round(quantiles[1], 2)),
      #                                                          paste0("Max: ", round(max_value, 2)))
      #                                             ) + ggtitle("Glasser") + 
      aes(fill = Estimate)) + scale_fill_gradientn(colours = cl, #midpoint = 0.05,
                                                    na.value = "#d6d6d6", 
                                                    breaks=c(min_value, quantiles[1], max_value),
                                                    labels=c(round(min_value,2), round(quantiles[1],2), round(max_value,2)),
                                                    limits=c(min_value, max_value)) + ggtitle("Heritability (FDR < 0.05)") +
    theme_bw() + 
    theme(text = element_text(size = 10)) + 
    theme(panel.grid = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),panel.border = element_blank()) +
    theme(legend.position = 'bottom') + theme(legend.key.size = unit(1, 'cm')) + theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) + 
    theme()
  
  ggsave(paste0(base_path, "analyze_results/herit_cortical_harm",harm,"_filt005.png"), plot = plot, width = 10, height = 8, dpi = 300)
  
  #-------------------------------------------plot last 19 ROIs----------------------------------------------------------
  ### plot subcortical
  #### Get all data
  subcor <- mmp_subcor[1:379] %>% str_remove(.,'_ROI') %>% str_split(.,'_',simplify = T) %>% data.frame()
  colnames(subcor) <- c('region','hemi')
  subcor <- subcor %>%
    tibble::rownames_to_column(var = "ROI") %>%
    mutate(ROI = paste0("",ROI)) %>%
    slice((n() - 18):n())  # 或者使用 tail(subcor, 19)
  
  aseg$data$label
  subcor$label <- paste0(str_to_title(subcor[,'hemi']), '-',str_to_title(subcor[,'region']))
  subcor$label[19] <- 'Brain-Stem'
  subcor$label[subcor$region=='thalamus'] <- paste0(subcor$label[subcor$region=='thalamus'],'-','Proper')
  subcor$label[subcor$region=='diencephalon'] <- paste0(str_to_title(subcor$hemi[subcor$region=='diencephalon']),'-VentralDC')
  subcor$label[subcor$region=='cerebellum'] <- paste0(subcor$label[subcor$region=='cerebellum'],'-','Cortex')
  
  # Joining roi_glasser with ivw_df to align data
  subcor <- subcor %>%
    mutate(ROI = as.integer(ROI)) %>%
    inner_join(scfc_harm_df, by = "ROI")
  
  ### cortical_position (i.e., how do we want the brain to be positioned)
  cortical_pos <- c("left lateral", "right lateral","left medial","right medial")

  testdata <- tibble(
    label = subcor$label,
    h2 = subcor$h_2
  )
  min_value <- min(testdata$h2, na.rm = TRUE)
  quantiles <- quantile(testdata$h2, probs = c(0.5), na.rm = TRUE)
  max_value <- max(testdata$h2, na.rm = TRUE)
  plot <- testdata %>%
    #group_by(groups) %>%
    ggplot() +
    geom_brain(atlas = aseg, aes(fill = h2))+
    scale_fill_gradientn(#low = muted("lightgreen"),
                         # # mid = "white",
                         # high = muted("darkgreen"),midpoint = 0, 
                         colours = cl,
                         na.value = "#d6d6d6", 
                         breaks=c(min_value, quantiles[1], max_value),
                         labels=c(round(min_value,2), round(quantiles[1],2), round(max_value,2)),
                         limits=c(min_value, max_value)) + ggtitle("Heritability (FDR p < 0.05)") + 
    #ggtitle(paste(y_annot[i]))+
    theme_bw() +
    theme(panel.grid = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),panel.border = element_blank())+
    theme(text = element_text(size = 15)) +
    theme(legend.position = 'bottom') + theme(legend.key.size = unit(1, 'cm')) + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
  
  ggsave(paste0(base_path, "analyze_results/herit_subcortical_harm",harm,"_filt005.png"), plot = plot, width = 10, height = 8, dpi = 300)
}
```

## Plot Box plots

```{r}
# first read network data
network_df <- read_excel(paste0(utilities_path, "CAB-NP_v1.1_Labels-ReorderedbyNetworks.xlsx")) %>%
  rename(roi_num = INDEX)
```


