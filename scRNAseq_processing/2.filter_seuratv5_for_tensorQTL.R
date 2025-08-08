
library(Seurat)
library(BPCells)
library(future)
library(vsn)
library(tidyverse)
library(SCpubr)
library(PCAtools)
# set this option when analyzing large datasets
options(future.globals.maxSize = 40000 * 1024^2) # 40Gb
# for new input_seurat assay slots
options(seurat.object.assay.version = "v5")

output_dir = "../../../OTAR2065_sc_eQTL/data/for_tensorQTL/"
dir.create(output_dir)
#biomartCacheClear()
# Load the seurat object

args = commandArgs(trailingOnly=TRUE)
message(length(args)," arguments provided")
if (length(args)<2) {
  stop("You need to detail 3 arguments: treatment and donor blacklist.n",
       call. = FALSE)
} else if (length(args) == 2) {
  treatment = args[1]
  donor_blacklist=args[2]
}

# # to test
# treatment="IFN"
# donor_blacklist="letw_5,lizq_3,zaie_1,romx_2,seru_7,qonc_2,iukl_1,curn_3,boqx_2,garx_2,sojd_3,yoch_6"

donor_blacklist = stringr::str_split(donor_blacklist,pattern = ",")[[1]]

# seurat data
input_seurat = readRDS(paste0("../../data/results/1.QC_v5/",
                              treatment,"_filtered_harmony/",
                              treatment,"_filtered_harmony.Rds"))

gc()

Layers(input_seurat[["RNA"]])
DefaultAssay(input_seurat) = "RNA"
Idents(input_seurat) = "cluster_full"

# protein-coding info
gene_info = read_csv("../../../resources/biomart/Homo_sapiens.GRCh38.111.genes.csv") %>%
  dplyr::rename(gene = gene_name) %>%
  dplyr::filter(gene_biotype == "protein_coding")

########
p1 =  SCpubr::do_BarPlot(input_seurat,
                         group.by = "Phase",
                         split.by = "cluster_full",
                         plot.title = "Proportion of cells per phase in each cluster",
                         position = "fill")

pdf(paste0(output_dir,"/",treatment,"_clusters_barplot_cycle_phase.pdf"),
    width = 9, height = 14)

plot(p1)
dev.off()


message("Naming proliferation categories")

input_seurat@meta.data$proliferation_status = "Not_proliferating"
if(treatment == "untreated"){
  input_seurat@meta.data[input_seurat@meta.data$cluster_full %in%  c(3,4) , "proliferation_status"] =  "Proliferating"
  
}
if(treatment == "LPS"){
  input_seurat@meta.data[input_seurat@meta.data$cluster_full %in%  c(4,5,7) , "proliferation_status"] =  "Proliferating"
  
}
if(treatment == "IFN"){  # much smaller proliferating clusters on IFN treated cells
  input_seurat@meta.data[input_seurat@meta.data$cluster_full %in%  c(4) , "proliferation_status"] =  "Proliferating"
  
}
message("Aggregating expression: subset to avoid memory issues")


pseudobulk = list()
metadata = list()
for(prolif in c("Not_proliferating","Proliferating")){
  subset_seurat = subset(x = input_seurat, subset = proliferation_status == prolif)
  pseudobulk[[prolif]] = AggregateExpression(subset_seurat, 
                                             assays = "RNA",
                                             return.seurat = FALSE,
                                             slot = "counts", # raw counts
                                             group.by = c("donor_id"))
  pseudobulk[[prolif]]  = as.data.frame(pseudobulk[[prolif]] )
  
  # remove genes that are not protein-coding
  pseudobulk[[prolif]] = pseudobulk[[prolif]] %>%
    tibble::rownames_to_column(var = "gene") %>%
    dplyr::filter(gene %in% gene_info$gene) %>%
    tibble::column_to_rownames(var = "gene")
    
    
  message("Adding metadata")
  
  metadata[[prolif]]  = subset_seurat@meta.data %>%
    dplyr::mutate(cols_names = paste(treatment,prolif,donor_id,sep = "_")) %>%
    dplyr::group_by(cols_names) %>%
    dplyr::reframe(count = n(),
                   proliferation_status = prolif,
                   donor_id=donor_id,
                   treatment = treatment) %>%
    dplyr::distinct() %>%
    dplyr::arrange(tolower(cols_names)) %>% # pseudobulk count columns are in alphabetical order (case insensitive)
    as.data.frame()
  
  colnames(pseudobulk[[prolif]] ) =metadata[[prolif]] $cols_names
  rm(subset_seurat)
  gc()
}
pseudobulk = Reduce(cbind, pseudobulk)
metadata = Reduce(rbind, metadata)
gene = rownames(pseudobulk)

rm(input_seurat)
gc()


message("Filtering by number of cells: remove donors with fewer than 100 cells per category")
message("Donors before filter in Not proliferating category: \n")
message(length(unique(metadata[metadata$proliferation_status=="Not_proliferating", "donor_id"])))
message("Donors before filter in Proliferating category: \n")
message(length(unique(metadata[metadata$proliferation_status=="Proliferating", "donor_id"])))

retain = metadata$count >=100
metadata = metadata[retain,]
pseudobulk = pseudobulk[,retain]

message("Donors after filter in Not proliferating category: \n")
message(length(unique(metadata[metadata$proliferation_status=="Not_proliferating", "donor_id"])))
message("Donors after filter in Proliferating category: \n")
message(length(unique(metadata[metadata$proliferation_status=="Proliferating", "donor_id"])))


message("Excluding donors from the blacklist")

retain = !(metadata$donor_id  %in% donor_blacklist)
metadata = metadata[retain,]
pseudobulk = pseudobulk[,retain]

message("Donors after filter in Not proliferating category: \n")
message(length(unique(metadata[metadata$proliferation_status=="Not_proliferating", 
                               "donor_id"])))
message("Donors after filter in Proliferating category: \n")
message(length(unique(metadata[metadata$proliferation_status=="Proliferating", 
                               "donor_id"])))


#log2 the summed counts applying size factors
message("Normalising to lib size, calculating log2(counts + 1)")

size_factors = colSums(pseudobulk)/mean(colSums(pseudobulk))
log_sum_counts = log2(t(t(pseudobulk)/size_factors) + 1)


# correlate size factors to number of cells
p = ggplot(data.frame(sf = size_factors,ncells=metadata$count),aes(x=ncells,y=sf)) +
  geom_point() + 
  theme_bw() +
  geom_smooth(method = "lm", se = FALSE)

pdf(file = paste0(output_dir,"/size_factors_vs_ncells_pseudobulk.pdf"),
    width = 5, height = 5)
plot(p)
dev.off()
#############
#log(ncells)
# To save expression per donor
message("Filtering again, now genes by expression levels")

# split tables according to proliferation status
metadata_list=list()
pseudobulk_list = list()
log_sum_counts_list = list()

write.table(metadata,paste0(output_dir,"/metadata.txt"))
write.table(pseudobulk,paste0(output_dir,"/pseudobulk_raw_counts.txt"))

pseudobulk_list[["Not_proliferating"]] = pseudobulk %>%
  as.data.frame() %>%
  dplyr::select(which(metadata$proliferation_status=="Not_proliferating"))
pseudobulk_list[["Proliferating"]] = pseudobulk %>%
  as.data.frame() %>%
  dplyr::select(which(metadata$proliferation_status=="Proliferating"))
log_sum_counts_list[["Not_proliferating"]] = log_sum_counts %>%
  as.data.frame() %>%
  dplyr::select(which(metadata$proliferation_status=="Not_proliferating"))
log_sum_counts_list[["Proliferating"]] = log_sum_counts %>%
  as.data.frame() %>%
  dplyr::select(which(metadata$proliferation_status=="Proliferating"))
metadata_list[["Not_proliferating"]] = metadata %>%
  as.data.frame() %>%
  dplyr::filter(proliferation_status == "Not_proliferating")
metadata_list[["Proliferating"]] = metadata %>%
  as.data.frame() %>%
  dplyr::filter(proliferation_status == "Proliferating")

cpm=list()
cpm[["Not_proliferating"]] =  apply(pseudobulk_list[["Not_proliferating"]],2, function(x) (x/sum(x))*1000000)
cpm[["Proliferating"]] =  apply(pseudobulk_list[["Proliferating"]],2, function(x) (x/sum(x))*1000000)

# To save list of genes per condition that pass the filters
gene_list = list()
# per proliferation cluster
# retain genes whose mean counts across donors are > 1CPM
# and retain genes that are present in at least 30% of the donors with at least 1CPM (rounding down)
filtered = cpm[["Not_proliferating"]] %>%
  as.data.frame() %>%
  dplyr::filter(rowMeans(.) >= 1 & (floor((rowSums(. > 0)/ncol(.))*100))>=30) 

log_sum_counts_list[["Not_proliferating"]] = log_sum_counts_list[["Not_proliferating"]][rownames(filtered),]

filtered = cpm[["Proliferating"]] %>%
  as.data.frame() %>%
  dplyr::filter(rowMeans(.) >= 1 & (floor((rowSums(. > 0)/ncol(.))*100))>=30) 
log_sum_counts_list[["Proliferating"]] = log_sum_counts_list[["Proliferating"]][rownames(filtered),]


message("Mean-SD plots")

# mean-sd plots
for(condition in names(log_sum_counts_list)){
  
  p1 =  vsn::meanSdPlot(as.matrix(log_sum_counts_list[[condition]]), plot = FALSE)$gg +theme_bw() +
    ggtitle("Libsize-scaled + log2 sum of raw counts (ranks)")
  
  
  p2 = vsn::meanSdPlot(as.matrix(log_sum_counts_list[[condition]]), plot = FALSE, ranks = FALSE)$gg +theme_bw() +
    ggtitle("Libsize-scaled + log2 sum of raw counts")
  
  pdf(file = paste0(output_dir,"/mean_sd_",condition,"_sum_pseudobulk.pdf"),
      width = 10, height = 7)
  plot(p1 + p2)
  dev.off()
}
# check that log_sum_counts_list has an expression more in line to what we'd expect from bulk
# with not too many genes on the left hand side of the curve (too low counts)
gc()

# write_rds(metadata_list,paste0(output_dir,"/metadata.rds"))
# write_rds(log_sum_counts_list,paste0(output_dir,"/log_sum_counts_list.rds"))
# metadata_list = read_rds(paste0(output_dir,"/metadata.rds"))
# log_sum_counts_list = read_rds(paste0(output_dir,"/log_sum_counts_list.rds"))
# save phenotype/expression table per condition
# bed format
#chr  start  end  gene_id  Sample1  Sample2

# look for gene id in biomart GRCh38 because that's what's being used in the RNA-seq file
# also in lifted over genotype files

message("Annotating Ensembl IDs for tensorQTL")
# ensembl v 111, GRCh38.p14 Downloaded from https://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz
# and filtered to genes
ensembl = read_csv("../../../resources/ENSEMBL_human_v111_2023/Homo_sapiens.GRCh38.111.genes.csv") %>%
  dplyr::select(seqname,start,end,gene_id,gene_name)

for(condition in names(log_sum_counts_list)){
  log_sum_counts_list[[condition]] = log_sum_counts_list[[condition]] %>%
    tibble::rownames_to_column(var = "gene_name") %>%
    dplyr::left_join(.,ensembl,by="gene_name") %>%
    dplyr::filter(!is.na(gene_id)) %>% # removing genes that are not found, around 900 from 12.5k
    dplyr::filter(seqname %in% c(1:22)) %>%  # Remove chromosomes not in 1:22
    dplyr::relocate(seqname,start,end,gene_id,gene_name) %>%
    dplyr::rename(`#chr` = seqname)
  message("There are ",nrow(log_sum_counts_list[[condition]])," genes left after filtering in ", condition)
}

message("Saving unscaled expression including gene name")

for(condition in names(log_sum_counts_list)){
  
  write.table(log_sum_counts_list[[condition]],paste0(output_dir,"/expr_sum_sizefactorsNorm_log2_",treatment, "_",condition ,".bed"),
              sep = "\t", quote = F, col.names = T, row.names = F)
  
}

# scaling for every gene across all samples
# so the vector for gene A [sample1, 2, ... sample x] has mean 0 sd 1, same for other genes
# that makes effect sizes comparable
message("Scaling expression and saving")

for(condition in names(log_sum_counts_list)){
  fixed_cols = log_sum_counts_list[[condition]] %>%
    dplyr::select("#chr","start","end","gene_id") %>%
     dplyr::mutate(end = start + 1) # fixing end for tensorQTL to be start + 1
          
  numeric_df = log_sum_counts_list[[condition]] %>%
    dplyr::select(!c(`#chr`,start,end,gene_id, gene_name)) %>%
    t() %>%
    as.data.frame() %>%
    mutate_all(as.numeric)
  scaled = base::scale(numeric_df) %>%
    t() %>%
    as.data.frame() %>%
    mutate_all(as.numeric) 
    
  scaled = cbind(fixed_cols,scaled)
  colnames(scaled) = c("#chr","start","end","gene_id",metadata_list[[condition]]$donor_id)    
  write.table(scaled,paste0(output_dir,"/expr_sum_sizefactorsNorm_log2_scaled_centered_",treatment, "_", condition,".bed"),
              sep = "\t", quote = F, col.names = T, row.names = F)
  
}

message("Saving metadata without PCs")

for(condition in names(metadata_list)){
  
  write.table(metadata_list[[condition]],paste0(output_dir,"/metadata_noPCs_",treatment, "_", condition,".txt"),
              sep = "\t", quote = F, col.names = T, row.names = F)
  
}

