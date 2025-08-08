# 5.1.2.2 Inspect Diff expr phagocytosis results

library(decoupleR)
library(tidyverse)
library(pheatmap)
library(ggplotify)
library(patchwork)
library(ggrepel)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(qvalue)
library(msigdbr)

source("./functions.R")
options(future.globals.maxSize = (40000 * 1024 ^ 2)) #(~40000 Gb)

treatment_cols =  c(untreated = "#8D918B", IFN = "#3A5683", LPS = "#F8766D")

outdir="../../../data/results/phagocytosis/6.Diff_expr_limma_phagocytosis/"

# limma results

deg = list()
for(cont in c("untreated","IFN","LPS")){
  deg[[cont]] = read_csv(paste0(outdir,"/DiffExpr_log_fraction_phago_",
                                cont,".csv"))
  deg[[cont]]$contrast = cont
}

deg = do.call("rbind",deg)

deg = deg %>%
  dplyr::group_by(contrast) %>%
  dplyr::mutate(scaled_t = center_this(t)) %>%
  dplyr::mutate(scaled_t_flipped = -center_this(t)) %>%
  dplyr::mutate(normalised_t_higher_phago_dir = (scaled_t - min(scaled_t)) / (max(scaled_t) - min(scaled_t)),
                normalised_t_lower_phago_dir = (scaled_t_flipped - min(scaled_t_flipped)) / (max(scaled_t_flipped) - min(scaled_t_flipped))) %>%
  dplyr::mutate( diffexp = case_when(logFC > 0.5 & adj.P.Val < 0.05 ~ "UP",
                                     logFC < -0.5 & adj.P.Val < 0.05 ~ "DOWN",
                                     .default = "NO")) %>%
  ungroup()

table(deg$contrast, deg$adj.P.Val < 0.05)
table(deg$contrast, deg$diffexp)

# way fewer differentially expressed genes in LPS!
# LPS seems to decrease phagocytosis in macrophage, but lit says it increases phago in microglia?
# DOWN    NO    UP
# IFN        2357  8238  2252
# LPS          58 12764    37
# untreated  2750  7836  2601

ensembl = symbol_to_ensembl(deg$symbol)
ensembl = enframe(ensembl, name = "symbol", value = "gene")
deg = deg %>%
  dplyr::left_join(.,ensembl) %>%
  distinct() %>%
  dplyr::filter(!is.na(gene))

write_csv(deg,paste0(outdir,"/DiffExpr_phagocytosis_within_treatment_all.csv"))

# check LPS -  seems to have higher spread in gene expression between same donor in different pools, and 
# smaller range of scaled phagocytosis fraction
# so careful with those results


## volcano plots
plist = list()
for(treatment in unique(deg$contrast)){
  plist[[treatment]] = deg  %>% 
    dplyr::filter(contrast == treatment) %>%
  dplyr::filter(diffexp == "NO") %>%   # Filter the rows that meet the conditions
    dplyr::sample_frac(0.05) %>%         # Sample 5% of the filtered data
    dplyr::bind_rows(deg[deg$contrast==treatment,]  %>% filter(!(diffexp == "NO"))) %>%  # Combine with the rest of the data
    dplyr::mutate(pval=-log10(adj.P.Val)) %>% 
    
    dplyr::mutate(pval=case_when(pval=="Inf"~350,
                                 .default = pval)) %>% # substitute for very extreme pvalue, at the limit of the precision
    
    ggplot(aes(
      x = logFC,
      y = pval,
      col = diffexp,
      label = symbol
    )) +
    geom_point() +
    theme_minimal() +
    scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),breaks = c(0,10,50,100,200,400)) +
    ggrepel::geom_text_repel(color = "black") +
    scale_color_manual(values = c("DOWN"="#3A5683",
                                  "NO"= "#A7A9AE",
                                  "UP"="#74121D")) +
    geom_vline(xintercept = c(0.5, 0.5),linetype = "dashed", color = "grey") +
    geom_hline(yintercept = -log10(0.05),linetype = "dashed", color = "grey") +
    ylab("-log10(adj. p-val)") + 
    labs(fill = "Diff. expr.") +
    ggtitle(treatment)
}


png(paste0(outdir,"/volcanos_res_additive.png"),
    width = 15, height = 4, res = 400,units = "in", type = "cairo")
patchwork::wrap_plots(plist, guides = "collect")
dev.off()



#  t statistic - normalised as if for CELLECT

cellect_DEG = deg %>% 
  dplyr::select(gene, contrast, normalised_t_higher_phago_dir,normalised_t_lower_phago_dir) %>%
  group_by(contrast) %>%
  tidyr::pivot_wider(names_from = c(contrast),
                     values_from =c(normalised_t_higher_phago_dir,normalised_t_lower_phago_dir)) %>%
  tidyr::unnest(cols = c(normalised_t_higher_phago_dir_untreated, normalised_t_lower_phago_dir_untreated,
                         normalised_t_higher_phago_dir_IFN, normalised_t_lower_phago_dir_IFN,
                         normalised_t_higher_phago_dir_LPS,normalised_t_lower_phago_dir_LPS)) %>%
  replace_na(list(normalised_t_higher_phago_dir_untreated = 0, normalised_t_lower_phago_dir_untreated = 0,
                  normalised_t_higher_phago_dir_IFN = 0,normalised_t_lower_phago_dir_IFN = 0,
                  normalised_t_higher_phago_dir_LPS = 0, normalised_t_lower_phago_dir_LPS = 0)) %>%
  dplyr::left_join(deg[,c("symbol","gene")]) %>%
  distinct()

summary(cellect_DEG$normalised_t_higher_phago_dir_LPS)

###### plots of changing effect sizes of phagocytosis between treatments
# would probably need to run a proper interaction analysis

p1 = deg  %>%
  pivot_wider(id_cols = symbol, names_from = contrast, values_from = t, values_fill = 0) %>%
  # dplyr::mutate(correlation = case_when(IFN > 0 & LPS > 0 ~ "direct correlation",
  #                                       IFN < 0 & LPS < 0 ~ "direct correlation",
  #                                       IFN < 0 & LPS > 0 ~ "inverse correlation",
  #                                       IFN > 0 & LPS < 0 ~ "inverse correlation",
  #                                       .default = "no correlation")) %>%
  dplyr::group_by(symbol) %>%
  dplyr::mutate(foralpha = (abs(IFN - untreated) + max(abs(IFN),abs(untreated))) #to highlight the values far from 0
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(foralpha2 = (abs(foralpha) - min(abs(foralpha))) / (max(abs(foralpha)) - min(abs(foralpha)))) %>%
  ggplot(aes(
    x = IFN,
    y = untreated,
    label = symbol,
    alpha = foralpha2
  )) +
  geom_point() +
  theme_minimal() +
  # scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),breaks = c(0,10,50,100,200,400)) +
  ggrepel::geom_text_repel(color = "black") +
  
  geom_vline(xintercept = c(-0.25, 0.25),linetype = "dashed", color = "grey") +
  geom_hline(yintercept = c(-0.25, 0.25),linetype = "dashed", color = "grey") +
  xlab("t statistic for phagocytosis DEA in IFN") + 
  ylab("t statistic for phagocytosis DEA in untreated") + 
  labs(fill = "Diff. expr.") +
  ggtitle("phagocytosis DEA effects IFN vs LPS")

p2 = deg  %>%
  pivot_wider(id_cols = symbol, names_from = contrast, values_from = t, values_fill = 0) %>%
  
  dplyr::group_by(symbol) %>%
  dplyr::mutate(foralpha = (abs(LPS - untreated) + max(abs(LPS),abs(untreated))) #to highlight the values far from 0
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(foralpha2 = (abs(foralpha) - min(abs(foralpha))) / (max(abs(foralpha)) - min(abs(foralpha)))) %>%
  ggplot(aes(
    x = LPS,
    y = untreated,
    label = symbol,
    alpha = foralpha2
  )) +
  geom_point() +
  theme_minimal() +
  # scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),breaks = c(0,10,50,100,200,400)) +
  ggrepel::geom_text_repel(color = "black") +
  
  geom_vline(xintercept = c(-0.25, 0.25),linetype = "dashed", color = "grey") +
  geom_hline(yintercept = c(-0.25, 0.25),linetype = "dashed", color = "grey") +
  xlab("t statistic for phagocytosis DEA in LPS") + 
  ylab("t statistic for phagocytosis DEA in untreated") + 
  labs(fill = "Diff. expr.") +
  ggtitle("phagocytosis DEA effects IFN vs LPS")


p3 = deg  %>%
  pivot_wider(id_cols = symbol, names_from = contrast, values_from = t, values_fill = 0) %>%
  dplyr::mutate(correlation = case_when(IFN > 0 & LPS > 0 ~ "direct correlation",
                                        IFN < 0 & LPS < 0 ~ "direct correlation",
                                        IFN < 0 & LPS > 0 ~ "inverse correlation",
                                        IFN > 0 & LPS < 0 ~ "inverse correlation",
                                        .default = "no correlation")) %>%
  dplyr::group_by(symbol) %>%
  dplyr::mutate(foralpha = (abs(IFN - LPS) + max(abs(IFN),abs(LPS))) #to highlight the values far from 0
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(foralpha2 = (abs(foralpha) - min(abs(foralpha))) / (max(abs(foralpha)) - min(abs(foralpha)))) %>%
  ggplot(aes(
    x = IFN,
    y = LPS,
    label = symbol,
    alpha = foralpha2
  )) +
  geom_point() +
  theme_minimal() +
  # scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),breaks = c(0,10,50,100,200,400)) +
  ggrepel::geom_text_repel(color = "black") +
  
  geom_vline(xintercept = c(-0.25, 0.25),linetype = "dashed", color = "grey") +
  geom_hline(yintercept = c(-0.25, 0.25),linetype = "dashed", color = "grey") +
  xlab("t statistic for phagocytosis DEA in IFN") + 
  ylab("t statistic for phagocytosis DEA in LPS") + 
  labs(fill = "Diff. expr.") +
  ggtitle("phagocytosis DEA effects IFN vs LPS")

# right: higher risk effect in IFN left: lower risk in IFN
# up: higher risk effect in LPS down: lower risk LPS

png(paste0(outdir,"/DEA_PRS_t_statistic_correlation.png"),
    width = 10, height = 4,res = 400,units = "in" )
plot(p1 + p2 + p3)
dev.off()

png(paste0(outdir,"/DEA_PRS_t_statistic_correlation_IFN_LPS.png"),
    width = 10, height = 10,res = 400,units = "in" )
plot(p3)
dev.off()


############# pathway analysis and others ##########
# take top 10% values and check genes

top_10_percent_scores = cellect_DEG %>%
  summarise(across(
    starts_with("normalised"),                    # Apply to all columns (or use specific column names)
    ~ quantile(.x, probs = 0.9),     # Calculate the 90th percentile (top 10% quantile)
    .names = "top_10QS_{col}"       # Naming the new columns
  ))

deg_list = list()
for(treat in unique(deg$contrast)){
  deg_list[[treat]] = deg %>%
    dplyr::filter(contrast %in% treat) %>%
    dplyr::filter(normalised_t_higher_phago_dir >= unlist(top_10_percent_scores[,
                                                                               paste0("top_10QS_normalised_t_higher_phago_dir_",
                                                                                      treat)]) | normalised_t_lower_phago_dir >= unlist(top_10_percent_scores[,
                                                                                                                                                             paste0("top_10QS_normalised_t_lower_phago_dir_",treat)]) )
  
  
}


#  set of genes that have normalised t in top 10%
genes = do.call("rbind", deg_list) %>%
  dplyr::select(symbol) %>%
  distinct() %>%
  unlist()

# check t statistic - filtered for top 10% normalised t

deg_t_filtered =  deg %>%
  dplyr::filter(symbol %in% genes) %>%
  dplyr::select(symbol, scaled_t, contrast) %>%
  dplyr::arrange(symbol) %>%
  tidyr::pivot_wider(id_cols = c(symbol), names_from = contrast, values_from = scaled_t, values_fill =0) %>%
  # remove rows with NA
  tidyr::drop_na() %>%
  column_to_rownames(var = "symbol") %>% 
  as.matrix() # much fewer genes than in counts, as expected

# full list of t-statistic, out of the phagocytosis DEG
deg_t=  deg %>%
  dplyr::select(symbol, scaled_t, contrast) %>%
  dplyr::arrange(symbol) %>%
  tidyr::pivot_wider(id_cols = c(symbol), names_from = contrast, values_from = scaled_t, values_fill =0) %>%
  # remove rows with NA
  tidyr::drop_na() %>%
  column_to_rownames(var = "symbol") %>% 
  as.matrix() 

### check count data from limma - AveExpr (average log2 expression over all samples)
# 
counts =list()
for(treat in unique(deg$contrast)){
  
  counts[[treat]] = deg %>%
    dplyr::filter(contrast %in% treat) %>%
    dplyr::select(symbol, contrast, AveExpr) 
}

counts = do.call("rbind", counts)

# wide format, then matrix
counts = counts %>%
  dplyr::arrange(symbol) %>%
  tidyr::pivot_wider(id_cols = c(symbol), names_from = contrast, values_from = AveExpr) %>%
  # select genes that are in the top deg list
  dplyr::filter(symbol %in% genes) %>%
  # remove rows with NA
  tidyr::drop_na() %>%
  column_to_rownames(var = "symbol") %>% 
  as.matrix() 

######### TF activity inference ######
net = get_collectri(organism='human', split_complexes=FALSE)
net

# from counts
TF_scores = run_ulm(mat=counts, net=net, .source='source', .target='target',
                    .mor='mor', minsize = 5)


# check TF function in microglia
# https://www.jci.org/articles/view/90604

# from scaled t-statistic
# Run ulm
# top 10%
TF_scores_t_filtered = run_ulm(mat=deg_t_filtered, net=net, .source='source', .target='target',
                               .mor='mor', minsize = 5)
# using full list of genes
TF_scores_t = run_ulm(mat=deg_t, net=net, .source='source', .target='target',
                      .mor='mor', minsize = 5)
# I'm interested in the genes that are DEG across phagocytosis, so I'm more interested in the t-statistic
# for that DEG than in the counts themselves
# though both are related to the normalised t that goes into CELLECT

table(TF_scores_t_filtered$p_value < 0.05, TF_scores_t_filtered$condition) # 137 in total
table(TF_scores_t$p_value < 0.05, TF_scores_t$condition) # 220 in total
table(TF_scores$p_value < 0.05, TF_scores$condition) # 421 in total


############### plot heatmap ##########
n_tfs = 30
# Get top tfs with more variable means across clusters
variable_tfs = TF_scores_t_filtered %>%
  group_by(source) %>%
  summarise(std = sd(score)) %>%
  arrange(-abs(std)) 

sign_tfs =  TF_scores_t_filtered %>%
  dplyr::filter(p_value < 0.05) %>%
  dplyr::arrange(p_value) %>%
  distinct(source) %>%
  pull(source)

tfs = variable_tfs %>%
  dplyr::filter(source %in% sign_tfs) %>%
  head(n_tfs) %>%
  pull(source)

toplot = TF_scores_t_filtered %>%
  dplyr::filter(source %in% tfs) %>%
  # scale scores across treatments, per TF
  dplyr::group_by(source) %>%
  dplyr::mutate(scaled_score = scale_this(score)) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_wider(id_cols = 'condition', names_from = 'source',
                     values_from = 'scaled_score') %>%
  column_to_rownames('condition') %>%
  as.matrix() %>%
  t()

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("#3A5683", "white","#74121D"))(palette_length)

my_breaks = c(seq(min(toplot), 0, length.out=ceiling(palette_length/2) + 1),
              seq(0.05, max(toplot), length.out=floor(palette_length/2)))

# Plot - scaled to highlight differences between treatments (so red may not have a positive enrichment score!!! e.g. MYC)
pheatmap(toplot, border_color = NA, color=my_color, breaks = my_breaks) 

# now from top p-vals

top_sign_tfs = TF_scores_t_filtered %>%
  dplyr::filter(p_value < 0.05) %>%
  dplyr::arrange(p_value) %>%
  distinct(source) %>%
  head(n_tfs) %>%
  pull(source)  

toplot = TF_scores_t_filtered %>%
  dplyr::filter(source %in% top_sign_tfs) %>%
  # scale scores across treatments, per TF
  dplyr::group_by(source) %>%
  dplyr::mutate(scaled_score = scale_this(score)) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_wider(id_cols = 'condition', names_from = 'source',
                     values_from = 'scaled_score') %>%
  column_to_rownames('condition') %>%
  as.matrix() %>%
  t()

my_breaks = c(seq(min(toplot), 0, length.out=ceiling(palette_length/2) + 1),
              seq(0.05, max(toplot), length.out=floor(palette_length/2)))

# from top pvals - scaled to highlight differences between treatments (so red may not have a positive enrichment score!!!)

p_t_filtered = ggplotify::as.ggplot(pheatmap(toplot, border_color = NA, 
                                             color=my_color, breaks = my_breaks))  +
  ggtitle("Filtered t")
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11065001/


# now full t-stats from phagocytosis assay
top_sign_tfs = TF_scores_t %>%
  dplyr::filter(p_value < 0.05) %>%
  dplyr::arrange(p_value) %>%
  distinct(source) %>%
  head(n_tfs) %>%
  pull(source)  
toplot = TF_scores_t %>%
  dplyr::filter(source %in% top_sign_tfs) %>%
  # scale scores across treatments, per TF
  dplyr::group_by(source) %>%
  dplyr::mutate(scaled_score = scale_this(score)) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_wider(id_cols = 'condition', names_from = 'source',
                     values_from = 'scaled_score') %>%
  column_to_rownames('condition') %>%
  as.matrix() %>%
  t()

my_breaks = c(seq(min(toplot), 0, length.out=ceiling(palette_length/2) + 1),
              seq(0.05, max(toplot), length.out=floor(palette_length/2)))

# from full centered t statistic - scaled to highlight differences between treatments (so red may not have a positive enrichment score!!!)

p_t = ggplotify::as.ggplot(pheatmap(toplot, border_color = NA, color=my_color, breaks = my_breaks)) +
  ggtitle("Unfiltered t")
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11065001/

### counts


top_sign_tfs = TF_scores %>%
  dplyr::filter(p_value < 0.05) %>%
  dplyr::arrange(p_value) %>%
  distinct(source) %>%
  head(n_tfs) %>%
  pull(source)  

toplot = TF_scores %>%
  dplyr::filter(source %in% top_sign_tfs) %>%
  # scale scores across treatments, per TF
  dplyr::group_by(source) %>%
  dplyr::mutate(scaled_score = scale_this(score)) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_wider(id_cols = 'condition', names_from = 'source',
                     values_from = 'scaled_score') %>%
  column_to_rownames('condition') %>%
  as.matrix() %>%
  t()


my_breaks = c(seq(min(toplot), 0, length.out=ceiling(palette_length/2) + 1),
              seq(0.05, max(toplot), length.out=floor(palette_length/2)))

p_counts = ggplotify::as.ggplot(pheatmap(toplot, border_color = NA,
                                         color=my_color, breaks = my_breaks))  +
  ggtitle("Filtered counts")
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11065001/

# In the mic-resting GRN, APOE was regulated by CEBPA, CEBPB(x), CEBPD, IRF7, 
#and JUND TFs. In the mic-activated and mic-inflammatory GRNs, 
#APOE was regulated by CEBPB, CEBPD, and JUNB, with IRF7 also involved in the 
#mic-activated state. 

p_t + p_t_filtered + p_counts


########## not scaling
# p_t scores unscaled reflects which targets are upregulated (pos score) or downregulated
# (neg score) by phagocytosis
top_sign_tfs = TF_scores_t %>%
  dplyr::filter(p_value < 0.05) %>%
  dplyr::arrange(p_value ) %>%
  distinct(source) %>%
  head(n_tfs) %>%
  pull(source)  

toplot = TF_scores_t %>%
  dplyr::filter(source %in% top_sign_tfs) %>%
  tidyr::pivot_wider(id_cols = 'condition', names_from = 'source',
                     values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix() %>%
  t()

my_breaks = c(seq(min(toplot), 0, length.out=ceiling(palette_length/2) + 1),
              seq(0.05, max(toplot), length.out=floor(palette_length/2)))

p_t_unscaled = ggplotify::as.ggplot(pheatmap(toplot, border_color = NA, 
                                             color=my_color, breaks = my_breaks)) +
  ggtitle("Unfiltered t - unscaled")

pdf(paste0(outdir,"/TF_activity_plots_full_DEG_t_statistic_heatmap.pdf"), 
    height = 15, width = 10)
p_t + p_t_unscaled

dev.off()

# fixing heatmap

toplot = TF_scores_t %>%
  dplyr::filter(source %in% top_sign_tfs)  %>%
  dplyr::mutate(p_val_sig_plot =  c("***", "**", "*", "")[findInterval(p_value, c(0.001, 0.01, 0.05)) + 1]) %>%
  dplyr::arrange(p_value)

# hierarchical clustering based on score
toplot_wide = toplot %>%
  dplyr::select(source, condition, score) %>%
  tidyr::pivot_wider(names_from = condition, values_from = score, values_fill = 0) 

# Perform hierarchical clustering on 'source'
row_order = hclust(dist(toplot_wide[,-1]))$order
source_order = toplot_wide$source[row_order]

# Perform hierarchical clustering on 'condition'
column_order = hclust(dist(t(toplot_wide[,-1])))$order
condition_order = colnames(toplot_wide)[-1][column_order]

# Reorder the 'source' and 'condition' factors based on clustering
toplot$source = factor(toplot$source, levels = source_order)
toplot$condition = factor(toplot$condition, levels = condition_order)


## plot

p = ggplot(toplot, aes(x= condition , y=source, fill = score)) + 
  geom_raster() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), # lables vertical
        strip.text.y = element_blank()) +  #remove facet bar on y 
  scale_fill_gradient2(low="steelblue",mid = "white", high = "#74121D",name = "Enrichment score") +
  geom_text(aes(label = p_val_sig_plot), col = "white", size = 5) + 
  ggtitle("Transcription factor activity enrichment") +
  ylab("TF") +
  xlab("Treatment") +
  theme_minimal() + 
  theme(axis.ticks = element_blank(),
        panel.grid = element_blank()) +
  theme(legend.position="bottom")

p

##########
##########
##########
png(paste0(outdir,"/TF_activity_plots_full_DEG_t_statistic_unscaled_heatmap.png"), 
    height = 8, width = 3.5, res = 400, units = "in",)
plot(p)
dev.off()
# per treatment, which TF are overexpressed / underexpressed across phagocytosis scale
##########
##########
##########

p = list()
for(treat in c(unique(TF_scores_t$condition))){
  p[[treat]] = TF_scores_t %>%
    dplyr::filter(condition == treat) %>%
    dplyr::arrange(desc(abs(score))) %>%
    dplyr::filter(p_value < 0.05) %>%
    ggplot(aes(x = reorder(source,score), y = score)) + 
    geom_bar(aes(fill = score), stat = "identity") +
    scale_fill_gradient2(low = "#3A5683", high = "#74121D", 
                         mid = "whitesmoke", midpoint = 0) + 
    theme_minimal() +
    theme(axis.title = element_text(face = "bold", size = 12),
          axis.text.y = element_text(size =10),
          axis.text.x = element_text(size =10, face= "bold"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position = "none") +
    xlab(ifelse(
      treat == "untreated",
      "TFs",
      "")) +
    coord_flip() +
    ggtitle(treat)
  
}


pdf(paste0(outdir,"/TF_activity_plots_full_DEG_t_statistic_sign_scores_barplot.pdf"), 
    height = 12, width = 8)
(p$untreated + p$IFN + p$LPS) 
dev.off()


# write scores from t values

write_csv(TF_scores_t,paste0(outdir,"/TF_activity_scores_full_DEG_t_statistic.csv"))

############ plot volcano plots of target TF

interesting_TFs = c("PPARG","FOXO3","PGR","PAX7","SALL3")

p_list = list()
for(tf in interesting_TFs){
  
  
  targets = net %>%
    filter(source == tf) %>%
    arrange(target) %>%
    mutate(ID = target, color = "3") %>%
    column_to_rownames('ID')
  
  inter = sort(intersect(rownames(deg_t),rownames(targets)))
  targets = targets[inter, ]
  
  # join with full DEG info
  targets = targets %>%
    dplyr::rename(symbol = target) %>%
    dplyr::left_join(deg) %>%
    dplyr::rename(target = symbol)
  
  targets = targets %>%
    mutate(`regulatory direction` = case_when(mor > 0 & t > 0 & P.Value < 0.05 ~ 'activated',
                                              mor > 0 & t < 0  & P.Value < 0.05 ~ 'repressed',
                                              mor < 0 & t > 0  & P.Value < 0.05 ~ 'activated',
                                              mor < 0 & t < 0  & P.Value < 0.05 ~ 'repressed',
                                              .default = "NS" ))
  
  p_list[[tf]] = ggplot(targets, aes(x = logFC, y = -log10(P.Value), color = `regulatory direction`
  )) +
    geom_point(aes(alpha = scaled_t)) +
    scale_colour_manual(values = c("repressed" = "#3A5683", 
                                   "activated" = "#74121D",
                                   "NS"="#8D918B")) +
    geom_label_repel(aes(label = target), size=1) + 
    theme_minimal() +
    theme(legend.position = "bottom") +
    geom_vline(xintercept = 0, linetype = 'dotted') +
    geom_hline(yintercept = -log10(0.05), linetype = 'dotted') +
    ggtitle(tf) + 
    facet_grid(vars(contrast),scales = "free")
}

# blue means that the sign of multiplying the mor and t-value is negative, 
# meaning that these genes are “deactivating” the TF, 
# and red means that the sign is positive, meaning that these genes are “activating” the TF.


### check all APOE targets as well, and pull that info from sign values
APOE_TFs=  net %>%
  dplyr::filter(target == "APOE") %>%
  dplyr::arrange(target) %>%
  dplyr::mutate(TFs = source, color = "3") %>%
  column_to_rownames('source')

TF_scores_t = run_ulm(mat=deg_t, net=net, .source='source', .target='target',
                      .mor='mor', minsize = 5)

APOE_reg_scores = TF_scores_t %>%
  dplyr::filter(TF_scores_t$source %in% APOE_TFs$TFs) %>%
  dplyr::mutate(`regulatory direction` = case_when(p_value < 0.05 & score > 0 ~"targets activated",
                                                   p_value < 0.05 & score < 0 ~ "targets repressed",
                                                   .default = "NS"))

table(APOE_reg_scores$source,APOE_reg_scores$condition)

p = list()
for(treat in c(unique(APOE_reg_scores$condition))){
  p[[treat]] = APOE_reg_scores %>%
    dplyr::filter(condition == treat) %>%
    dplyr::arrange(desc(abs(score))) %>%
    ggplot(aes(x = reorder(source,score), y = score, fill =`regulatory direction` )) + 
    geom_bar( stat = "identity") +
    scale_fill_manual(values = c("targets activated" = "#74121D",
                                 "targets repressed" = "#3A5683", 
                                 "NS"="#8D918B")) +
    theme_minimal() +
    theme(axis.title = element_text(face = "bold", size = 12),
          axis.text.y = element_text(size =10),
          axis.text.x = element_text(size =10, face= "bold"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position = "none") +
    xlab(ifelse(
      treat == "untreated",
      "TFs",
      "")) +
    coord_flip() +
    ggtitle(treat)
  
}

##########
###########
##########
###### increased significance of regulators of AD-related genes in LPS? check for other AD candidates from jeremy
################

pdf(paste0(outdir,"/APOE_reg_TF_full_DEG_t_statistic_sign_scores_barplot.pdf"), 
    height = 12, width = 8)
(p$untreated + p$IFN + p$LPS) + patchwork::plot_annotation(title="TF scores of APOE regulators")
dev.off()

### checking phagocytosis candidates from Sam's screen
# pos score means KO phagocytose MORE - so that means downregulation in expression increases phagocytosis in our data (negative t statistic)
sam_phago = read_tsv("../../../../CRISPR/OTAR2065_phagocytosis_CRISPR/data/2024_07_sam_screen_results.tsv") %>%
  dplyr::arrange(FDR) %>%
  dplyr::rename(symbol = id)

# checking matching directionality in untreated
directionality = sam_phago %>%
  dplyr::left_join(deg) %>%
  dplyr::filter(contrast == "untreated") %>%
  dplyr::mutate(direction = case_when((t < 0 & Score > 0 & adj.P.Val < 0.05) | (t > 0 & Score < 0 & adj.P.Val < 0.05) ~ "same direction",
                                      (t < 0 & Score < 0 & adj.P.Val < 0.05) | (t > 0 & Score > 0 & adj.P.Val < 0.05)  ~ "opposite direction",
                                      .default = "NS DiffExpr"),
                direction_no_significance = case_when((t < 0 & Score > 0 ) | (t > 0 & Score < 0 ) ~ "same direction",
                                      (t < 0 & Score < 0 ) | (t > 0 & Score > 0 )  ~ "opposite direction"))
  

table(directionality$direction)
table(directionality$direction_no_significance)
# coin toss for chance of same direction

p = ggplot(directionality,aes(x = -Score*(-log10(FDR)), y = t, color = direction,label = symbol)) +
  geom_point(size = 5) + 
  geom_text_repel() + 
  scale_color_manual(values = c("same direction" = "#74121D",
                               "opposite direction" = "#3A5683", 
                               "NS DiffExpr"="#8D918B")) +
  theme_minimal()
# plotting to highlight those that are really significant in the CRISPR
pdf(paste0(outdir,"/sam_CRISPR_genes_diff_expr_directionality_dotplot.pdf"), 
    height = 4, width = 5)
plot(p)
dev.off()
# divide between TF and targets

net_sam_TF = net %>%
  dplyr::filter(source %in% sam_phago$symbol) # some

net_sam_targets = net %>%
  dplyr::filter(target %in% sam_phago$symbol) # some



tf_reg_scores = list()
for(gene in unique(net_sam_targets$target)){
  tf_subset = net_sam_targets %>%
    dplyr::filter(target == gene) %>%
    dplyr::arrange(target) %>%
    dplyr::mutate(TFs = source, color = "3") %>%
    column_to_rownames('source')
  
  tf_reg_scores[[gene]] = TF_scores_t %>%
    dplyr::filter(TF_scores_t$source %in% tf_subset$TFs) %>%
    dplyr::mutate(`regulatory direction` = case_when(p_value < 0.05 & score > 0 ~"targets activated",
                                                     p_value < 0.05 & score < 0 ~ "targets repressed",
                                                     .default = "NS"),
                  sam_gene = gene)
  
  
  
  
  
  
  
}

tf_reg_scores = do.call("rbind",tf_reg_scores) %>%
  dplyr::filter(`regulatory direction` !="NS")


p = tf_reg_scores %>%
  dplyr::group_by(source) %>%
  dplyr::filter(any(`regulatory direction` !="NS")) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(abs(score))) %>%
  ggplot(aes(x = reorder(source,score), y = score, fill =`regulatory direction` )) + 
  geom_bar( stat = "identity") +
  scale_fill_manual(values = c("targets activated" = "#74121D",
                               "targets repressed" = "#3A5683", 
                               "NS"="#8D918B")) +
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.y = element_text(size =10),
        axis.text.x = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  coord_flip() +
  facet_wrap(vars(condition))

pdf(paste0(outdir,"/all_Sam_targets_reg_TF_full_DEG_t_statistic_sign_scores_barplot.pdf"), 
    height = 12, width = 5)

plot(p)

dev.off()

write_csv(tf_reg_scores,paste0(outdir,"/all_Sam_targets_reg_TF_full_DEG_t_statistic_sign_scores.csv"))

### check Jeremy AD candidate genes

ad_candidate_genes = read_csv("../../../../resources/AD_PD_gene_sets_Andrew/Set1_Jeremy_candidates.csv") %>%
  # dplyr::filter(str_detect(disease,"AD")) %>%
  pull(gene_sym)


# divide between TF and targets

net_ad_TF = net %>%
  dplyr::filter(source %in% ad_candidate_genes) 

net_ad_targets = net %>%
  dplyr::filter(target %in% ad_candidate_genes) 



tf_reg_scores = net_ad_targets %>%
  dplyr::left_join(TF_scores_t)  %>%
  dplyr::mutate(dummy =paste(source,condition, score, sep = "_")) %>%
  dplyr::distinct(dummy, .keep_all = TRUE) %>%
  dplyr::select(-dummy) %>%
  tidyr::drop_na() %>%
  dplyr::mutate(`regulatory direction` = case_when(p_value < 0.05 & score > 0 ~"targets activated",
                                                   p_value < 0.05 & score < 0 ~ "targets repressed",
                                                   .default = "NS")) %>%
  distinct()


p = tf_reg_scores %>%
  dplyr::group_by(source) %>%
  dplyr::filter(any(`regulatory direction` !="NS")) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(abs(score))) %>%
  ggplot(aes(x = reorder(source,score), y = score, fill =`regulatory direction` )) + 
  geom_bar( stat = "identity") +
  scale_fill_manual(values = c("targets activated" = "#74121D",
                               "targets repressed" = "#3A5683", 
                               "NS"="#8D918B")) +
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.y = element_text(size =10),
        axis.text.x = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  coord_flip() +
  ggtitle("Scores of TF regulators of AD candidate genes") +
  facet_wrap(vars(condition))

pdf(paste0(outdir,"/AD_candidate_genes_targets_reg_TF_full_DEG_t_statistic_sign_scores_barplot.pdf"), 
    height = 12, width = 5)

plot(p)

dev.off()

# now check TFs among AD candidates

tf_reg_scores = TF_scores_t %>%
  dplyr::filter(TF_scores_t$source %in% net_ad_TF$source) %>%
  dplyr::mutate(`regulatory direction` = case_when(p_value < 0.05 & score > 0 ~"targets activated",
                                                   p_value < 0.05 & score < 0 ~ "targets repressed",
                                                   .default = "NS"))


p = list()
for(treat in c(unique(tf_reg_scores$condition))){
  p[[treat]] = tf_reg_scores %>%
    dplyr::filter(condition == treat) %>%
    dplyr::arrange(desc(abs(score))) %>%
    ggplot(aes(x = reorder(source,score), y = score, fill =`regulatory direction` )) + 
    geom_bar( stat = "identity") +
    scale_fill_manual(values = c("targets activated" = "#74121D",
                                 "targets repressed" = "#3A5683", 
                                 "NS"="#8D918B")) +
    theme_minimal() +
    theme(axis.title = element_text(face = "bold", size = 12),
          axis.text.y = element_text(size =10),
          axis.text.x = element_text(size =10, face= "bold"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position = "none") +
    xlab(ifelse(
      treat == "untreated",
      "TFs",
      "")) +
    coord_flip() +
    ggtitle(treat)
  
}



pp = (p$untreated + p$IFN + p$LPS) + patchwork::plot_annotation(title=paste0("TF AD candidates scores"))

pdf(paste0(outdir,"/AD_candidate_genes_reg_TF_full_DEG_t_statistic_sign_scores_barplot.pdf"), 
    height = 2, width = 3)
plot(pp)
dev.off()

### investigate plain candidate gene enrichment in DEA





## jeremy genes
jeremy_ad_candidates = read_table("../../../../resources/Jeremy_medrXiv_AD_candidate_genes.txt")
jeremy_set1_ad_candidates = read_csv("../../../../resources/AD_PD_gene_sets_Andrew/Set1_Jeremy_candidates.csv")
jeremy_set3_pd_candidates = read_csv("../../../../resources/AD_PD_gene_sets_Andrew/Set3_manually_curated.csv") %>%
  dplyr::filter(str_detect(disease,"PD")) %>%
  pull(gene_name)
# lead phagocytosis sam
sam_crispr = read_tsv("../../../../CRISPR/OTAR2065_phagocytosis_CRISPR/data/2024_07_sam_screen_results.tsv") %>%
  dplyr::arrange(FDR) %>%
  slice_head(n = 30) %>%
  pull(id)

# lead phagocytosis kampmann
kampmann_crispr_phago = read_csv("../../../../resources/CRISPRbrain_Kampmanm/kampmann_CRISPRi_iTF_microglia_phagocytosis.csv") %>%
  dplyr::mutate(FDR = qvalue(`P Value`)$qvalues) %>% # nothing seems to be significant
  dplyr::arrange(`P Value`) %>%
  slice_head(n = 30) %>%
  pull(Gene)

# lead immune activation kampmann
kampmann_crispr_immune = read_csv("../../../../resources/CRISPRbrain_Kampmanm/kampmann_CRISPRi_iTF_microglia_immune_activation.csv") %>%
  dplyr::mutate(FDR = qvalue(`P Value`)$qvalues) %>% # nothing seems to be significant
  dplyr::arrange(`P Value`) %>%
  slice_head(n = 30) %>%
  pull(Gene)
#### format my genes

# add entrez ids
deg$entrez = symbol_to_entrez(deg$symbol)

untreated_sorted =  deg %>%
  dplyr::filter(contrast == "untreated") %>%
  dplyr::arrange(desc(scaled_t)) %>%
  pull(t)

names(untreated_sorted) = deg %>%
  dplyr::filter(contrast == "untreated") %>%
  dplyr::arrange(desc(scaled_t)) %>%
  pull(symbol)
untreated_sorted = untreated_sorted[!is.na(names(untreated_sorted))]
# untreated_sorted = abs(untreated_sorted)

ifn_sorted =  deg %>%
  dplyr::filter(contrast == "IFN") %>%
  dplyr::arrange(desc(scaled_t)) %>%
  pull(t)

names(ifn_sorted) = deg %>%
  dplyr::filter(contrast == "IFN") %>%
  dplyr::arrange(desc(scaled_t)) %>%
  pull(symbol)
ifn_sorted = ifn_sorted[!is.na(names(ifn_sorted))]
# ifn_sorted = abs(ifn_sorted)

lps_sorted =  deg %>%
  dplyr::filter(contrast == "LPS") %>%
  dplyr::arrange(desc(scaled_t)) %>%
  pull(t)

names(lps_sorted) = deg %>%
  dplyr::filter(contrast == "LPS") %>%
  dplyr::arrange(desc(scaled_t)) %>%
  pull(symbol)
lps_sorted = lps_sorted[!is.na(names(lps_sorted))]
# lps_sorted = abs(lps_sorted)


pathways = list()
pathways[["jeremy_ad_candidates"]] = jeremy_ad_candidates$symbol
pathways[["set1_jeremy_ad_candidates"]] = jeremy_set1_ad_candidates$gene_sym # larger
pathways[["set3_jeremy_pd_candidates"]] = jeremy_set3_pd_candidates
pathways[["lead_phagocytosis_sam"]] = sam_crispr 
pathways[["lead_phagocytosis_kampmann"]] = kampmann_crispr_phago
pathways[["lead_immune_kampmann"]] = kampmann_crispr_immune 

### GSEA

ifn_res = list()
lps_res = list()
untreated_res = list()

untreated_res[["custom_genesets"]] = fgsea::fgsea(pathways=pathways, 
                                                  stats=untreated_sorted,
                                                  # scoreType="pos",
                                                  eps=0) %>%
  tidyr::as_tibble() %>%
  dplyr::arrange(desc(NES))

lps_res[["custom_genesets"]] = fgsea::fgsea(pathways=pathways, 
                                            stats=lps_sorted,
                                            # scoreType="pos",
                                            eps=0) %>%
  tidyr::as_tibble() %>%
  dplyr::arrange(desc(NES))


ifn_res[["custom_genesets"]] = fgsea::fgsea(pathways=pathways, 
                                            stats=ifn_sorted,
                                            # scoreType="pos",
                                            eps=0) %>%
  tidyr::as_tibble() %>%
  dplyr::arrange(desc(NES))


ifn_res[["custom_genesets"]]  = ifn_res[["custom_genesets"]] %>%
  dplyr::mutate(contrast = "IFN phagocytosis")
lps_res[["custom_genesets"]]  = lps_res[["custom_genesets"]] %>%
  dplyr::mutate(contrast = "LPS phagocytosis")
untreated_res[["custom_genesets"]] = untreated_res[["custom_genesets"]] %>%
  dplyr::mutate(contrast = "untreated phagocytosis")

res_df = rbind(ifn_res[["custom_genesets"]] , 
               lps_res[["custom_genesets"]],
               untreated_res[["custom_genesets"]]  )


p = res_df %>%
  dplyr::filter(pathway != "jeremy_ad_candidates") %>%
  dplyr::mutate(contrast = factor(contrast,levels = c("untreated phagocytosis","IFN phagocytosis","LPS phagocytosis")),
                pathway = case_when(pathway == "set3_jeremy_pd_candidates" ~ "PD candidates",
                                    pathway == "set1_jeremy_ad_candidates" ~ "AD candidates",
                                    pathway == "lead_phagocytosis_sam" ~ "Lead phagocytosis CRISPR",
                                    pathway == "lead_phagocytosis_kampmann" ~ "Lead phagocytosis CRISPR (Kampmann)",
                                    pathway == "lead_immune_kampmann" ~ "Lead immune CRISPR (Kampmann)"),
                minus_log10_pval = -log10(padj),
                p_val_sig_plot =  c("***", "**", "*", "")[findInterval(padj, c(0.001, 0.01, 0.05)) + 1]) %>%
  dplyr::mutate(pathway = factor(pathway,levels = c( "AD candidates", "PD candidates",
                                                     "Lead phagocytosis CRISPR", "Lead phagocytosis CRISPR (Kampmann)","Lead immune CRISPR (Kampmann)"))) %>%
  ggplot( aes(y=pathway, x=contrast, fill = NES)) + 
  geom_raster() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), # lables vertical
        strip.text.y = element_blank()) +  #remove facet bar on y 
  scale_fill_gradient2(low = "steelblue",mid = "white", high = "#74121D",name = "NES") +
  geom_text(aes(label = p_val_sig_plot), col = "white", size = 5) + 
  ggtitle("GSEA enrichment in candidate gene sets") +
  theme_minimal() + 
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12)) +
  # facet_grid(cols = vars(`gene list` )) +
  theme(legend.position="right")


png(paste0(outdir,"/GSEA_candidate_genes_tileplot.png"), 
    height = 4.5,width =10,units = "in", res = 500)
plot(p)
dev.off()
# top 30 most significant phago genes from kampmann are enriched in the untreated set (sorted)
# checking directionalities as with sam above in dotplt
directionality_kampmann = read_csv("../../../../resources/CRISPRbrain_Kampmanm/kampmann_CRISPRi_iTF_microglia_phagocytosis.csv") %>%
  dplyr::rename(symbol = Gene) %>%
  dplyr::mutate(t_kampmann = Phenotype * -log10(`P Value`)) %>%
  dplyr::arrange(desc(t_kampmann)) %>%
  dplyr::left_join(deg) %>%
  dplyr::filter(contrast == "untreated") %>%
  dplyr::mutate(direction = case_when((t < 0 & t_kampmann > 0 & adj.P.Val < 0.05) | (t > 0 & t_kampmann < 0 & adj.P.Val < 0.05) ~ "same direction",
                                      (t < 0 & t_kampmann < 0 & adj.P.Val < 0.05) | (t > 0 & t_kampmann > 0 & adj.P.Val < 0.05)  ~ "opposite direction",
                                      .default = "NS DiffExpr"),
                direction_no_significance = case_when((t < 0 & t_kampmann > 0 ) | (t > 0 & t_kampmann < 0 ) ~ "same direction",
                                                      (t < 0 & t_kampmann < 0 ) | (t > 0 & t_kampmann > 0 )  ~ "opposite direction"))


table(directionality_kampmann$direction)
table(directionality_kampmann$direction_no_significance)
# coin toss for chance of same direction, as with Sam

p = ggplot(directionality_kampmann,aes(x = t_kampmann, y = t, color = direction,label = symbol)) +
  geom_point(size = 5) + 
  geom_text_repel() + 
  scale_color_manual(values = c("same direction" = "#74121D",
                                "opposite direction" = "#3A5683", 
                                "NS DiffExpr"="#8D918B")) +
  theme_minimal()
# plotting to highlight those that are really significant in the CRISPR
pdf(paste0(outdir,"/kampmann_phago_CRISPR_genes_diff_expr_directionality_dotplot.pdf"), 
    height = 12, width = 8)
plot(p)
dev.off()
  
# do the genes that follow the same direction have anything in common?
#### do the phagocytosis genes that don't match directionality belong to a certain pathway?

directionality_kampmann %>%
  dplyr::filter(direction == "same direction") %>%
  pull(symbol) %>%
  paste(., collapse = " ") # mostly kinase GO, SPI1 TF

directionality_kampmann %>%
  dplyr::filter(direction == "opposite direction") %>%
  pull(symbol) %>%
  paste(., collapse = " ") # same, kinase GO, ZNF TF

directionality %>%
  dplyr::filter(direction == "same direction") %>%
  pull(symbol) %>%
  paste(., collapse = " ") # isopentenyl & mevalonate GO, lewy body disease

directionality %>%
  dplyr::filter(direction == "opposite direction") %>%
  pull(symbol) %>%
  paste(., collapse = " ") # apoptosis!!!!!, vesicle GO



############
### checking in gprofiler
# https://biit.cs.ut.ee/gprofiler/gost

  
  #########

# checking with directionality in the custom genesets
pathways[["untreated_pos"]] = deg %>%
  # dplyr::filter(contrast == "untreated" & adj.P.Val < 0.05) %>%
  dplyr::filter(symbol %in% unlist(deg_list$untreated[deg_list$untreated$scaled_t>0,"symbol"])) %>% # in top 10% significant genes
  pull(symbol)
pathways[["untreated_neg"]] = deg %>%
  # dplyr::filter(contrast == "untreated" & adj.P.Val < 0.05) %>%
  dplyr::filter(symbol %in% unlist(deg_list$untreated[deg_list$untreated$scaled_t<0,"symbol"])) %>% # in top 10% significant genes
  pull(symbol)
pathways[["LPS_pos"]] = deg %>%
  # dplyr::filter(contrast == "LPS" & adj.P.Val < 0.05) %>%
  dplyr::filter(symbol %in% unlist(deg_list$LPS[deg_list$LPS$scaled_t>0,"symbol"])) %>% # in top 10% significant genes
  pull(symbol)

pathways[["LPS_neg"]] = deg %>%
  # dplyr::filter(contrast == "LPS" & adj.P.Val < 0.05) %>%
  dplyr::filter(symbol %in% unlist(deg_list$LPS[deg_list$LPS$scaled_t<0,"symbol"])) %>% # in top 10% significant genes
  pull(symbol)

pathways[["IFN_pos"]] = deg %>%
  # dplyr::filter(contrast == "IFN" & adj.P.Val < 0.05) %>%
  dplyr::filter(symbol %in% unlist(deg_list$IFN[deg_list$IFN$scaled_t>0,"symbol"])) %>% # in top 10% significant genes
  pull(symbol)

pathways[["IFN_neg"]] = deg %>%
  # dplyr::filter(contrast == "IFN" & adj.P.Val < 0.05) %>%
  dplyr::filter(symbol %in% unlist(deg_list$IFN[deg_list$IFN$scaled_t<0,"symbol"]) )%>% # in top 10% significant genes
  pull(symbol)

jeremy_ad_candidates = read_table("../../../../resources/Jeremy_medrXiv_AD_candidate_genes.txt") %>%
  dplyr::arrange(desc(total_score)) %>%
  pull(total_score)
names(jeremy_ad_candidates) = read_table("../../../../resources/Jeremy_medrXiv_AD_candidate_genes.txt") %>%
  dplyr::arrange(desc(total_score)) %>%
  pull(symbol)
# lead phagocytosis sam
sam_crispr = read_tsv("../../../../CRISPR/OTAR2065_phagocytosis_CRISPR/data/2024_07_sam_screen_results.tsv") %>%
  dplyr::mutate(t = Score * -log10(FDR)) %>%
  dplyr::arrange(desc(t)) %>%
  pull(t)
names(sam_crispr) =  read_tsv("../../../../CRISPR/OTAR2065_phagocytosis_CRISPR/data/2024_07_sam_screen_results.tsv") %>%
  dplyr::mutate(t = Score * -log10(FDR)) %>%
  dplyr::arrange(desc(t)) %>%
  pull(id)

# lead phagocytosis kampmann
kampmann_crispr_phago = read_csv("../../../../resources/CRISPRbrain_Kampmanm/kampmann_CRISPRi_iTF_microglia_phagocytosis.csv") %>%
  dplyr::mutate(t = Phenotype * -log10(`P Value`)) %>%
  dplyr::arrange(desc(t)) %>%
  pull(t)

names(kampmann_crispr_phago) = read_csv("../../../../resources/CRISPRbrain_Kampmanm/kampmann_CRISPRi_iTF_microglia_phagocytosis.csv") %>%
  dplyr::mutate(t = Phenotype * -log10(`P Value`)) %>%
  dplyr::arrange(desc(t)) %>%
  pull(Gene)

# lead immune activation kampmann
kampmann_crispr_immune = read_csv("../../../../resources/CRISPRbrain_Kampmanm/kampmann_CRISPRi_iTF_microglia_immune_activation.csv") %>%
  dplyr::mutate(t = Phenotype * -log10(`P Value`)) %>%
  dplyr::arrange(desc(t)) %>%
  pull(t)

names(kampmann_crispr_immune) = read_csv("../../../../resources/CRISPRbrain_Kampmanm/kampmann_CRISPRi_iTF_microglia_immune_activation.csv") %>%
  dplyr::mutate(t = Phenotype * -log10(`P Value`)) %>%
  dplyr::arrange(desc(t)) %>%
  pull(Gene)

jeremy_ad_candidates_res = list()
sam_crispr_res = list()
kampmann_crispr_immune_res = list()
kampmann_crispr_phago_res = list()

jeremy_ad_candidates_res = fgsea::fgsea(pathways=pathways, 
                                        stats=jeremy_ad_candidates,
                                        scoreType="pos", # jeremy stats are all positive
                                        eps=0) %>%
  tidyr::as_tibble() %>%
  dplyr::arrange(desc(NES))

sam_crispr_res = fgsea::fgsea(pathways=pathways, 
                              stats=sam_crispr,
                              # scoreType="pos",
                              eps=0) %>%
  tidyr::as_tibble() %>%
  dplyr::arrange(desc(NES))


kampmann_crispr_immune_res= fgsea::fgsea(pathways=pathways, 
                                         stats=kampmann_crispr_immune,
                                         # scoreType="pos",
                                         eps=0) %>%
  tidyr::as_tibble() %>%
  dplyr::arrange(desc(NES))

kampmann_crispr_phago_res= fgsea::fgsea(pathways=pathways, 
                                        stats=kampmann_crispr_phago,
                                        # scoreType="pos",
                                        eps=0) %>%
  tidyr::as_tibble() %>%
  dplyr::arrange(desc(NES))


## nothing significant


# human hallmark pathways #########
# with directionality
deg$entrez = symbol_to_entrez(deg$symbol)

untreated_sorted =  deg %>%
  dplyr::filter(contrast == "untreated") %>%
  dplyr::arrange(desc(scaled_t)) %>%
  pull(scaled_t)

names(untreated_sorted) = deg %>%
  dplyr::filter(contrast == "untreated") %>%
  dplyr::arrange(desc(scaled_t)) %>%
  pull(entrez)
untreated_sorted = untreated_sorted[!is.na(names(untreated_sorted))]

ifn_sorted =  deg %>%
  dplyr::filter(contrast == "IFN") %>%
  dplyr::arrange(desc(scaled_t)) %>%
  pull(scaled_t)

names(ifn_sorted) = deg %>%
  dplyr::filter(contrast == "IFN") %>%
  dplyr::arrange(desc(scaled_t)) %>%
  pull(entrez)
ifn_sorted = ifn_sorted[!is.na(names(ifn_sorted))]

lps_sorted =  deg %>%
  dplyr::filter(contrast == "LPS") %>%
  dplyr::arrange(desc(scaled_t)) %>%
  pull(scaled_t)

names(lps_sorted) = deg %>%
  dplyr::filter(contrast == "LPS") %>%
  dplyr::arrange(desc(scaled_t)) %>%
  pull(entrez)
lps_sorted = lps_sorted[!is.na(names(lps_sorted))]


pathways = msigdbr("human", category="H")
pathways = split(as.character(pathways$entrez_gene), pathways$gs_name)

ifn_res[["gsea_msigdb_hallmark"]] = fgsea::fgsea(pathways=pathways, 
                                                 stats=ifn_sorted,
                                                 eps=0) %>%
  tidyr::as_tibble() %>%
  dplyr::arrange(desc(NES)) %>%
  dplyr::mutate(pathway = str_remove(pathway,"HALLMARK_")) %>%
  dplyr::mutate(pathway = str_replace_all(pathway,"_", " "))  %>%
  dplyr::mutate(pathway = str_to_sentence(pathway))


lps_res[["gsea_msigdb_hallmark"]] = fgsea::fgsea(pathways=pathways, 
                                                 stats=lps_sorted,
                                                 eps=0) %>%
  tidyr::as_tibble() %>%
  dplyr::arrange(desc(NES)) %>%
  dplyr::mutate(pathway = str_remove(pathway,"HALLMARK_")) %>%
  dplyr::mutate(pathway = str_replace(pathway,"_", " "))  %>%
  dplyr::mutate(pathway = str_to_sentence(pathway))

untreated_res[["gsea_msigdb_hallmark"]] = fgsea::fgsea(pathways=pathways, 
                                                       stats=untreated_sorted,
                                                       eps=0) %>%
  tidyr::as_tibble() %>%
  dplyr::arrange(desc(NES)) %>%
  dplyr::mutate(pathway = str_remove(pathway,"HALLMARK_")) %>%
  dplyr::mutate(pathway = str_replace(pathway,"_", " "))  %>%
  dplyr::mutate(pathway = str_to_sentence(pathway))

p2 = ifn_res[["gsea_msigdb_hallmark"]] %>%
  dplyr::filter(padj < 0.05) %>%
  ggplot(aes(reorder(pathway, NES), NES)) +
  geom_col() +
  coord_flip() +
  labs(x="Pathway", y="") + 
  theme_minimal() + 
  ggtitle("IFN phagocytosis")

p3 = lps_res[["gsea_msigdb_hallmark"]] %>%
  dplyr::filter(padj < 0.05) %>%
  ggplot(aes(reorder(pathway, NES), NES)) +
  geom_col() +
  coord_flip() +
  labs(x="", y="Normalized Enrichment Score") + 
  theme_minimal() + 
  ggtitle("LPS phagocytosis")

p1 = untreated_res[["gsea_msigdb_hallmark"]] %>%
  dplyr::filter(padj < 0.05) %>%
  ggplot(aes(reorder(pathway, NES), NES)) +
  geom_col() +
  coord_flip() +
  labs(x="", y="") + 
  theme_minimal() + 
  ggtitle("Untreated phagocytosis")

png(paste0(outdir,"/msigdb_hallmark_gsea_with_directionality.png"),
    width = 12, height = 5, res = 400,units = "in", type = "cairo")
(p1 + p2 + p3 ) + patchwork::plot_annotation(title = "DE genes: Hallmark pathways NES from GSEA")
dev.off()

# # which genes in these pathways
ifn_res[["gsea_msigdb_hallmark"]] = ifn_res[["gsea_msigdb_hallmark"]] %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::select("pathway", "leadingEdge","padj","NES","size") %>%
  tidyr::unnest(cols = c(leadingEdge)) %>%
  dplyr::rename(entrez = leadingEdge) %>%
  dplyr::inner_join(deg %>%
                      dplyr::filter(contrast == "IFN") %>%
                      dplyr::select(symbol,entrez)) %>%
  dplyr::group_by(pathway) %>%
  dplyr::reframe(padj = padj, NES = NES,size=size,
                 leading_edge = paste(symbol, collapse = ", ")) %>%
  distinct() %>%
  dplyr::arrange(padj) %>%
  dplyr::mutate(contrast = "IFN")

lps_res[["gsea_msigdb_hallmark"]] = lps_res[["gsea_msigdb_hallmark"]] %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::select("pathway", "leadingEdge","padj","NES","size") %>%
  tidyr::unnest(cols = c(leadingEdge)) %>%
  dplyr::rename(entrez = leadingEdge) %>%
  dplyr::inner_join(deg %>%
                      dplyr::filter(contrast == "LPS") %>%
                      dplyr::select(symbol,entrez)) %>% 
  dplyr::group_by(pathway) %>%
  dplyr::reframe(padj = padj, NES = NES,size=size,
                 leading_edge = paste(symbol, collapse = ", ")) %>%
  distinct() %>%
  dplyr::arrange(padj)  %>%
  dplyr::mutate(contrast = "LPS")

untreated_res[["gsea_msigdb_hallmark"]] = untreated_res[["gsea_msigdb_hallmark"]] %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::select("pathway", "leadingEdge","padj","NES","size") %>%
  tidyr::unnest(cols = c(leadingEdge)) %>%
  dplyr::rename(entrez = leadingEdge) %>%
  dplyr::inner_join(deg %>%
                      dplyr::filter(contrast == "untreated") %>%
                      dplyr::select(symbol,entrez)) %>%
  dplyr::group_by(pathway) %>%
  dplyr::reframe(padj = padj, NES = NES,size=size,
                 leading_edge = paste(symbol, collapse = ", ")) %>%
  distinct() %>%
  dplyr::arrange(padj) %>%
  dplyr::mutate(contrast = "Untreated")

# save results

rbind(ifn_res[["gsea_msigdb_hallmark"]] ,
      lps_res[["gsea_msigdb_hallmark"]] ,
      untreated_res[["gsea_msigdb_hallmark"]]) %>%
  write_csv(.,file=paste0(outdir,"/msigdb_halmmark_gsea_with_directionality_genelists.csv"))


### what are the genes of the IFN pathway targets of
df = rbind(ifn_res[["gsea_msigdb_hallmark"]] ,
           lps_res[["gsea_msigdb_hallmark"]] ,
           untreated_res[["gsea_msigdb_hallmark"]]) %>%
  dplyr::filter(pathway %in% c("Inflammatory response","Interferon gamma response","Mtorc1 signaling"))

net %>% 
  dplyr::filter(target %in% str_trim(str_split_1(df[1,]$leading_edge, pattern = ","))) %>%
  dplyr::group_by(source) %>%
  dplyr::summarise(n_targets = n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(n_targets))
# targets
targets = net %>% 
  dplyr::filter(target %in% str_trim(str_split_1(df[1,]$leading_edge, pattern = ","))) 
# sources
sources = net %>% 
  dplyr::filter(source %in% str_trim(str_split_1(df[1,]$leading_edge, pattern = ","))) 
#### progeny pathways

net =get_progeny(organism = 'human', top = 500)
net

progeny =run_mlm(mat=deg_t, net=net, .source='source', .target='target',
                 .mor='weight', minsize = 5)

toplot = progeny %>%
  dplyr::mutate(p_val_sig_plot =  c("***", "**", "*", "")[findInterval(p_value, c(0.001, 0.01, 0.05)) + 1]) %>%
  dplyr::arrange(p_value)

# hierarchical clustering based on score
toplot_wide = toplot %>%
  dplyr::select(source, condition, score) %>%
  tidyr::pivot_wider(names_from = condition, values_from = score, values_fill = 0) 

# Perform hierarchical clustering on 'source'
row_order = hclust(dist(toplot_wide[,-1]))$order
source_order = toplot_wide$source[row_order]

# Perform hierarchical clustering on 'condition'
column_order = hclust(dist(t(toplot_wide[,-1])))$order
condition_order = colnames(toplot_wide)[-1][column_order]

# Reorder the 'source' and 'condition' factors based on clustering
toplot$source = factor(toplot$source, levels = source_order)
toplot$condition = factor(toplot$condition, levels = condition_order)

p = ggplot(toplot, aes(x= condition , y=source, fill = score)) + 
  geom_raster() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), # lables vertical
        strip.text.y = element_blank()) +  #remove facet bar on y 
  scale_fill_gradient2(low="steelblue",mid = "white", high = "#74121D",name = "Enrichment score") +
  geom_text(aes(label = p_val_sig_plot), col = "white", size = 5) + 
  ggtitle("Pathway activity enrichment") +
  ylab("Pathway") +
  xlab("Treatment") +
  theme_minimal() + 
  theme(axis.ticks = element_blank(),
        panel.grid = element_blank()) +
  theme(legend.position="bottom")


png(paste0(outdir,"/progeny_pathway_enrichment_scaled_t.png"),
    width = 6, height = 5, res = 400,units = "in", type = "cairo")
plot(p)
dev.off()
# https://saezlab.github.io/decoupleR/articles/pw_bk.html

# EGFR: regulates growth, survival, migration, apoptosis, proliferation, and differentiation in mammalian cells
# JAK-STAT: involved in immunity, cell division, cell death, and tumor formation.
# MAPK: integrates external signals and promotes cell growth and proliferation.
# NFkB: regulates immune response, cytokine production and cell survival.
# PI3K: promotes growth and proliferation.
# TGFb: involved in development, homeostasis, and repair of most tissues.
# Trail: induces apoptosis.

pathway ='MAPK'

p_list = list()
for(treat in c("untreated","IFN","LPS")){
  df =net %>%
    filter(source == pathway) %>%
    arrange(target) %>%
    mutate(ID = target, color = "3") %>%
    column_to_rownames('target')
  inter =sort(intersect(rownames(deg_t),rownames(df)))
  df =df[inter, ]
  df['scaled_t'] =deg_t[inter,treat ]
  df =df %>%
    mutate(color = if_else(weight > 0 & scaled_t > 0, '2', color)) %>%
    mutate(color = if_else(weight > 0 & scaled_t < 0, '1', color)) %>%
    mutate(color = if_else(weight < 0 & scaled_t > 0, '1', color)) %>%
    mutate(color = if_else(weight < 0 & scaled_t < 0, '2', color))
  
  p_list[[treat]] = ggplot(df, aes(x = weight, y = scaled_t, color = color)) + geom_point() +
    scale_colour_manual(values = c("steelblue","#74121D","white")) +
    geom_label_repel(aes(label = ID)) + 
    theme_minimal() +
    theme(legend.position = "none") +
    geom_vline(xintercept = 0, linetype = 'dotted') +
    geom_hline(yintercept = 0, linetype = 'dotted') +
    ggtitle(treat)
}

plot(p_list$untreated | p_list$IFN | p_list$LPS) + patchwork::plot_annotation(title = pathway)






