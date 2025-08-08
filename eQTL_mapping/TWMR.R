# adapted from https://github.com/eleporcu/TWMR/blob/master/MR.R
#  from https://www.nature.com/articles/s41467-019-10936-0
# check also this: https://github.com/soreshkov/pyTWMR
# Explanation of how to prepare LD matrices and eQTL matrices is right before equation 3
# We ran a multivariable MR analysis for each ~16,000 gene, where we conditioned each gene's causal effect
# on the potential causal effects of all of its neighboring genes. Let us consider now one focal gene.
# We need to select instrument SNPs and exposure genes for the multivariate MR analysis that is destined
# to elucidate the focal gene’s multivariate causal effect on the outcome (GWAS). To this end, we first consider
# all the independent eQTLs for the focal gene with conditional P<1×10−3.
### In the reviewer report: We use GCTA only to select the independent SNPs and then we use the raw (univariate)
### effect estimates for such SNPs in the MR model – as the formula requires
### we can adapt this to select from the eQTL - pheno GWAS coloc the top colocalising SNPs

# Next, we include as exposures
# all the genes for which the selected SNPs are eQTLs. Finally, we extend the instruments to include all SNPs
# that are eQTLs for any of the exposure genes. Note that genes that do not share eQTLs with the focal gene
# do not alter the focal gene’s multivariate causal effect, hence do not need to be considered here.

### In the reviewer report: On average, for example for height, we included 8 SNPs and 2 genes in each MR
### model.
# To avoid numerical instability in our multiple regression model, we pruned SNPs that are in high LD (r2>0.1)
# also, in "Correlation between genes"
# To avoid numerical instability caused by near-colinearity in our multiple regression model and making choices
# between co-regulated genes, we removed one gene from each pair of genes with r2≥0.4.
library(tidyverse)
source("./functions.R")
options(scipen = 999) # prevents from showing very small numbers as 0



args = commandArgs(trailingOnly = TRUE)

# Checking arguments
if (length(args) < 4) {
  stop("4 arguments must be supplied", call. = FALSE)
}

focus_gene = as.character(args[1]) # gene to analyse as exposure for TWMR
treatment = as.character(args[2]) # treatment, for coloc, eQTL and GWAS files
phenotype = as.character(args[3]) # phenotype for GWAS
output_path = as.character(args[4]) # Output directory. Will create a treatment subdirectory

# to test
# focus_gene = "TREM2"
# treatment = "untreated"
# phenotype="phagocytosis"
# output_path="../../data/results/8.5.eQTL_MR/TWMR/output/"

### should probably separate input preparation (with LD filters and so on) from actual run, and save intermediate
# LD and effect files



# reading input and output files
read_input_files = function(gene, input_path) {
  # Read eQTL and GWAS effects file
  
  # read file containing eQTL (gene columns) and GWAS  (BETA_GWAS, last columns)
  #  effects per variant (first column, annoyingly called "genes")
  filecluster = read.table(
    paste0(input_path, gene, ".matrix"),
    header = TRUE,
    sep = " ",
    dec = "."
  )
  beta = as.matrix(filecluster[, 2:(ncol(filecluster) - 1)])
  
  # extracting beta for GWAS, here called "gamma" , G in paper
  gamma = as.matrix(filecluster[, ncol(filecluster)])
  
  # loading LD matrix for SNP effects on the tested gene, here called "C". Diagonal left to right is always 1 (LD with itself)
  LDmatrix = read.table(
    paste0(input_path, gene, ".ld"),
    header = FALSE,
    sep = " ",
    dec = "."
  )
  C = as.matrix(LDmatrix[, 1:ncol(LDmatrix)])
  
  return(
    list(
      beta = beta,
      gamma = gamma,
      C = C,
      LDmatrix = LDmatrix,
      eQTL_GWAS_effects = filecluster
    )
  )
}
# Main function to run TWMR analysis
run_TWMR_calculations = function(beta, gamma, C, Ngwas, N_eQTLs, n_for_LD) {
  # beta,gamma,c: output of read_input_files()
  # gene: focal gene for which we are testing effects on phenotype (provided GWAS)
  # Ngwas: number of participants in GWAS
  # N_eQTLs: number of participants in eQTL
  # n_for_LD: N used for LD calculation (I use 1000 Genomes, they used UK10K https://github.com/eleporcu/TWMR/issues/4)
  
  
  # aggregating  (sum) eQTL effects per gene for all variants, and removing those with overall effect on gene of 0
  # here is called beta but in their paper is E
  x = colSums(abs(beta))
  remove = which(x == 0)
  if (length(remove) > 0) {
    beta = beta[, -remove]
  }
  beta = as.matrix(beta)
  
  
  # Solve equation 2: Calculate alpha
  C_inv = solve(C)
  S = t(beta) %*% C_inv %*% beta
  H = (1 - 1 / sqrt(n_for_LD)) * S + (1 / sqrt(n_for_LD)) * diag(nrow(S))
  # here is the solution to equation 2 in the paper
  alpha = solve(H) %*% (t(beta) %*% C_inv %*% gamma)
  alpha = as.vector(alpha)
  
  # Solve equation 3: Calculate standard error (SE)
  GCG_inv = t(beta) %*% C_inv %*% beta
  GCG_inv = (1 - 1 / sqrt(n_for_LD)) * GCG_inv + (1 / sqrt(n_for_LD)) * diag(nrow(GCG_inv))
  GCG_inv = solve(GCG_inv)
  
  df_dg = GCG_inv %*% t(beta) %*% C_inv
  df_dG = (GCG_inv %x% (t(gamma) %*% C_inv %*% ((beta %*% GCG_inv %*% t(beta)) %*% C_inv + diag(nrow(beta))
  ))) +
    ((-t(gamma) %*% C_inv %*% beta %*% GCG_inv) %x% (GCG_inv %*% t(beta) %*% C_inv))
  J = cbind(df_dG, df_dg)
  
  SEs = c(rep(1 / sqrt(N_eQTLs), length(beta[1, ]) * length(beta[, 1])), rep(1 / sqrt(Ngwas), length(gamma[, 1])))
  R = diag(ncol(beta) + 1)
  Sigma = (SEs %*% t(SEs)) * (C %x% R)
  V = J %*% Sigma %*% t(J)
  se = sqrt(V[1, 1])
  
  # Calculate p-value
  N = nrow(beta)
  Ngene = ncol(beta)
  # causal effect Z-statistic
  Z = alpha[1] / se
  pval = 2 * pnorm(abs(Z), lower.tail = FALSE)
  
  # Prepare output
  
  # finally, for the focal gene "gene" detailed as input, we get
  # the causal effect of the gene expression on the outcome: alpha
  # its standard error: SE
  # its p-value: p
  # the number of SNPs used in the model: Nsnps
  # the number of genes used in the model: Ngene
  # Return results
  return(list(
    alpha = alpha[1],
    SE = se,
    P = pval,
    Nsnps = N,
    Ngene = Ngene
  ))
}



# Function to load and process coloc results
process_coloc_results = function(treatment, focus_gene, coloc_res) {
  coloc_res = coloc_res[[treatment]]

    if (sum(stringr::str_split_i(names(coloc_res),"_",i=2) %in% focus_gene) == 0) {
      lead_coloc_var = NA # there were no colocs
    }
    if (sum(stringr::str_split_i(names(coloc_res),"_",i=2) %in% focus_gene) > 0) {
      coloc_res = coloc_res[[names(coloc_res)[stringr::str_split_i(names(coloc_res),"_",i=2) %in% focus_gene]]]
      
      
      if (is_null(coloc_res)) {
        lead_coloc_var = NA # some might be NULL because it didn't fill the minimum number of snps needed for coloc
      }
      if (!is_null(coloc_res)) {
      
      if (coloc_res$summary["PP.H4.abf"] < 0.70) {
        lead_coloc_var = NA # there were no significant colocs
        
      } else{
        # selecting SNP with largest SNP.PP.H4
        lead_coloc_var = coloc_res$results %>%
          dplyr::select("snp", "SNP.PP.H4") %>%
          slice_max(SNP.PP.H4, with_ties = FALSE, n = 1) %>% # if there are ties, return first
          dplyr::pull(snp)
      }
    }
  }
  
  
  return(lead_coloc_var)
}


load_nominal_results = function(treatment) {
  nominal_res = readr::read_delim(
    paste0(
      "../../data/results/tensorqtl/best_results/sum_sizefactorsNorm_log2_scaled_centered_",
      treatment,
      "_Not_proliferating_common_500kb_window_tensorQTL_nominal.txt"
    )
  )
  
  
  # mapping ensembl to symbol
  ensembl_to_symbol = ensembl_to_gene_name(
    IDs = unique(nominal_res$phenotype_id),
    IDFrom = "GENEID",
    IDTo = "SYMBOL"
  )$map %>%
    dplyr::rename(phenotype_id = From, gene = To) %>%
    as_tibble() %>%
    dplyr::distinct(phenotype_id, .keep_all = TRUE)
  
  nominal_res = nominal_res %>%
    dplyr::left_join(ensembl_to_symbol)
  # adding info and filtering
  return(nominal_res)
  
}

#### very important function: how eQTL variants are selected at first
extract_relevant_variants_eQTL = function(focus_gene,
                                          nominal_res,
                                          coloc_effect_all_genes) {
  lead_effect_focus_gene = nominal_res %>%
    dplyr::select(variant_id, gene, pval_nominal, slope) %>%
    dplyr::filter(gene == focus_gene) %>%
    dplyr::slice_min(pval_nominal, na_rm = TRUE) %>% # selecting only min p-val variants, could also select all significant ones and filter post-hoc by LD:
    # dplyr::filter(pval_nominal<0.00001,na_rm = TRUE) %>%
    distinct()
  # check other eGenes for lead eQTL
  lead_effect_other_genes = nominal_res %>%
    dplyr::select(variant_id, gene, pval_nominal, slope) %>%
    dplyr::filter(variant_id %in% lead_effect_focus_gene$variant_id) %>% ### this may need an additional filter
    distinct() %>%
    dplyr::filter(gene != unique(lead_effect_focus_gene$gene)) # remove results for focus gene
  
  # adding together - doesn't seem to be the case here - search for case with multiple eQTL
  # table_eGenes_per_variant = nominal_res %>%
  #   dplyr::group_by(variant_id) %>%
  #   dplyr::distinct() %>%
  #   dplyr::filter(pval_nominal < 0.00001) %>%
  #   dplyr::summarise(n_eGenes= n()) %>%
  #   dplyr::ungroup() %>%
  #   dplyr::arrange(desc(n_eGenes))
  #
  # hist(table_eGenes_per_variant$n_eGenes)
  
  
  # loading other significant eQTLs for these nearby genes
  extended_eQTLs_other_genes = nominal_res %>%
    dplyr::select(variant_id, gene, pval_nominal, slope) %>%
    dplyr::filter(gene %in% lead_effect_other_genes$gene) %>%
    dplyr::group_by(gene) %>%
    dplyr::filter(rank(pval_nominal, ties.method = "first") == 1) %>%
    dplyr::filter(pval_nominal < 0.000001, na_rm = TRUE) %>% # take min breaking ties, must have p < 10-6
    distinct()
  
  
  
  # adding everything together
  if (is.null(dim(coloc_effect_all_genes))) {
    lead_effect_focus_gene_other_genes = rbind(lead_effect_focus_gene,
                                               lead_effect_other_genes,
                                               extended_eQTLs_other_genes) %>%
      dplyr::distinct() %>%
      dplyr::select(variant_id, gene, slope) %>%
      tidyr::pivot_wider(
        id_cols = c(variant_id),
        names_from = gene,
        values_from = slope,
        values_fill = 0
      ) # for each gene, variants with missing effect are 0
  } else{
    lead_effect_focus_gene_other_genes = rbind(
      coloc_effect_all_genes,
      lead_effect_focus_gene,
      lead_effect_other_genes,
      extended_eQTLs_other_genes
    ) %>%
      dplyr::distinct() %>%
      dplyr::select(variant_id, gene, slope) %>%
      tidyr::pivot_wider(
        id_cols = c(variant_id),
        names_from = gene,
        values_from = slope,
        values_fill = 0
      ) # for each gene, variants with missing effect are 0
    
    
    
  }
  
  return(
    list(
      "lead_effect_focus_gene_other_genes" = lead_effect_focus_gene_other_genes,
      "lead_effect_focus_gene" =  lead_effect_focus_gene
    )
  )
  
}

#### another very important function: how LD is processed and filtered
load_and_process_ld_file = function(coloc_effect_all_genes,
                                    lead_effect_focus_gene_other_genes,
                                    lead_effect_focus_gene) {
  # loading LD:
  # Loading chromosome LD file from TopMED (more comprehensive than 1000G)
  # consider subseting it first to tested variants so it doesn't take that long to load
  ld = readr::read_delim(
    paste0(
      "/lustre/scratch125/humgen/resources/TopLD/release_28032022/EUR/SNV/EUR_chr",
      stringr::str_split_i(
        lead_effect_focus_gene$variant_id[1],
        pattern = "_",
        i = 1
      ),
      "_no_filter_0.2_1000000_LD.csv.gz"
    )
  ) %>%
    dplyr::filter(
      SNP1 %in% stringr::str_split_i(
        lead_effect_focus_gene_other_genes$variant_id,
        pattern = "_",
        i = 2
      ) |
        SNP2 %in% stringr::str_split_i(
          lead_effect_focus_gene_other_genes$variant_id,
          pattern = "_",
          i = 2
        )
    ) %>%
    dplyr::distinct()
  
  if (nrow(ld) == 0) {
    ld_variants = NA
  } else{
    ld_filtered = ld %>%
      dplyr::mutate(
        SNP1 = paste(
          stringr::str_split_i(
            lead_effect_focus_gene$variant_id[1],
            pattern = "_",
            i = 1
          ),
          str_split_i(Uniq_ID_1, pattern = ":", i =
                        1),
          str_split_i(Uniq_ID_1, pattern = ":", i =
                        2),
          str_split_i(Uniq_ID_1, pattern = ":", i =
                        3),
          sep = "_"
        ),
        SNP2 = paste(
          stringr::str_split_i(
            lead_effect_focus_gene$variant_id[1],
            pattern = "_",
            i = 1
          ),
          str_split_i(Uniq_ID_2, pattern = ":", i =
                        1),
          str_split_i(Uniq_ID_2, pattern = ":", i =
                        2),
          str_split_i(Uniq_ID_2, pattern = ":", i =
                        3),
          sep = "_"
        )
      ) %>%
      dplyr::mutate(SNP1 = pmin(SNP1, SNP2), # Sort pair to maintain consistency
                    SNP2 = pmax(SNP1, SNP2)) %>%
      dplyr::select(SNP1, SNP2, R2) %>%
      dplyr::distinct()  # Remove any duplicates
    
    
    # check lead variants first: if there are two in high LD, keep only one
    
    if (length(lead_effect_focus_gene$variant_id) > 1) {
      ld_focus_gene = ld_filtered %>%
        dplyr::filter(
          SNP1 %in% lead_effect_focus_gene$variant_id[1]  &
            SNP2 %in% lead_effect_focus_gene$variant_id[2] |
            SNP2 %in% lead_effect_focus_gene$variant_id[1]  &
            SNP1 %in% lead_effect_focus_gene$variant_id[2]
        )
      
      ### here may throw errors:
      if (ld_focus_gene$R2 > 0.2) {
        lead_effect_focus_gene = lead_effect_focus_gene %>%
          dplyr::slice_max(abs(slope), with_ties = FALSE, n = 1)
      }
      
    }
    
    table(lead_effect_focus_gene_other_genes$variant_id %in% ld_filtered$SNP1)
    table(lead_effect_focus_gene_other_genes$variant_id %in% ld_filtered$SNP2)
    
    
    # check that lead coloc and lead eQTL variants are present in LD file
    if (is.null(dim(coloc_effect_all_genes))) {
      vars_to_check = lead_effect_focus_gene %>%
        dplyr::slice_min(pval_nominal)
    } else{
      vars_to_check = rbind(
        coloc_effect_all_genes,
        lead_effect_focus_gene %>% dplyr::slice_min(pval_nominal)
      )
    }
    
    
    if (sum(!vars_to_check$variant_id %in% ld_filtered$SNP1) > 0) {
      message(
        "Careful: some of the focus genes lead colocs or lead eQTL variants are not present in the original LD file. Those variants will be considered not to be in LD with the rest ... "
      )
      message("This might over/under estimate the aggregate SNP effects on the gene. ")
      print(
        lead_effect_focus_gene_other_genes %>%
          dplyr::filter(variant_id %in% vars_to_check$variant_id) %>%
          dplyr::mutate(
            in_ld_file = case_when(variant_id %in% ld_filtered$SNP1 ~ TRUE, .default = FALSE)
          )
      )
    }
    
    
    ##### IMPORTANT
    
    # For those variants from the selected eQTL results not present here, set R2 to 0
    # this may have unexpected consequences
    
    if (sum(!lead_effect_focus_gene_other_genes$variant_id %in% ld_filtered$SNP1) >
        0) {
      ld_filtered2 = ld_filtered %>%
        # Shuffling the second pair to avoid identity R2 (that should be 1)
        
        bind_rows(
          data.frame(
            SNP1 = lead_effect_focus_gene_other_genes[!lead_effect_focus_gene_other_genes$variant_id %in% ld_filtered$SNP1, "variant_id"]$variant_id,
            SNP2 = sample(lead_effect_focus_gene_other_genes[!lead_effect_focus_gene_other_genes$variant_id %in% ld_filtered$SNP1, "variant_id"]$variant_id),
            R2 = 0
          )
        )
    } else{
      ld_filtered2 = ld_filtered
    }
    
    
    ld_filtered2 = ld_filtered2 %>%
      dplyr::mutate(SNP1 = pmin(SNP1, SNP2), # Force consistent SNP1-SNP2 order
                    SNP2 = pmax(SNP1, SNP2)) %>%
      dplyr::distinct(SNP1, SNP2, .keep_all = TRUE)  # Keep only the first appearance of each pair
    
    # check I'm not losing lead variants
    lead_effect_focus_gene_other_genes$variant_id %in% ld_filtered2$SNP1
    
    ####### IMPORTANT: subset to only variants in low LD ########
    # may need to change to D' because of this ####
    select_ld_variants = ld_filtered2 %>%
      dplyr::filter(R2 < 0.2, SNP1 != SNP2) %>%  # Keep only pairs below 0.2 (and exclude self-pairs)
      dplyr::mutate(keep = pmin(SNP1, SNP2)) %>%  # Pick the lower-named variant to keep
      dplyr::distinct(keep) %>%  # Remove duplicate pairs
      dplyr::pull(keep)  # Extract only the variants to keep
    
    # if all the variants were in LD > 0.2, select the colocalising variant and lead eqtl variant
    
    if (purrr::is_empty(select_ld_variants)) {
      if (is.null(dim(coloc_effect_all_genes))) {
        select_ld_variants = unique(c(lead_effect_focus_gene_other_genes$variant_id))
        
      } else{
        select_ld_variants = unique(
          c(
            lead_effect_focus_gene_other_genes$variant_id,
            coloc_effect_all_genes$variant_id
          )
        )
        
      }
    }
    
    # ensure I'm not selecting the lead eQTL variant or the lead coloc variant
    
    lead_effect_focus_gene_other_genes$variant_id %in% select_ld_variants
    
    
    
    ###### the following filter is absurd - it keeps many variants in LD if I filter in one column,
    # filter on both with & or | is also stupid
    # trim directly on wide matrix ?
    # ld_filtered2 = ld_filtered2 %>%
    #   dplyr::filter(SNP1 %in% select_ld_variants) %>%
    #   dplyr::filter(SNP2 %in% select_ld_variants)
    
    
    ##############
    
    # Making wide table before filtering
    # Create a full list of unique SNPs (ensuring symmetry)
    unique_snps = unique(c(ld_filtered2$SNP1, ld_filtered2$SNP2))
    
    # Expand the data to include self-correlations (SNP1 = SNP2, R2 = 1)
    self_corr = tibble(SNP1 = unique_snps,
                       SNP2 = unique_snps,
                       R2 = 1)
    
    # Add self-correlations and make sure all SNP pairs exist
    ld_data = ld_filtered2 %>%
      dplyr::bind_rows(self_corr) %>%  # Add diagonal R2 = 1
      tidyr::complete(SNP1 = unique_snps, SNP2 = unique_snps) %>%
      dplyr::mutate(XY = map2(SNP1, SNP2, ~ paste0(sort(c(
        .x, .y
      )), collapse = "_"))) %>%
      dplyr::group_by(XY) %>%
      tidyr::fill(R2, .direction = "downup") %>%  # Fill in missing values within groups
      dplyr::ungroup() %>%
      dplyr::select(-XY)
    
    # Ensure all pairs exist
    ####### at this point I think the matrix is not symmetrical
    # need to find values where SNP1 = SNP2 and flipped and R2 > 0
    # because when the matrix is made wide it is not symmetrical on the side below the diagonal
    
    
    # check presence lead variants
    lead_effect_focus_gene_other_genes$variant_id %in% ld_data$SNP1
    
    
    # Fill missing R² values: Set self-correlation to 1, and all others to 0
    ld_data = ld_data %>%
      dplyr::mutate(R2 = replace_na(R2, 0))  # Fill missing R² values with 0, can have unexpected consequences
    
    # Convert to wide format
    ld_matrix = ld_data %>%
      tidyr::pivot_wider(
        id_cols = SNP1,
        names_from = SNP2,
        values_from = R2,
        values_fn = max
      ) %>% # taking max R2 in case of duplicates
      dplyr::rename(SNP = SNP1)
    
    ###### Important:
    # trim wide matrix to retain only set of variants in low LD
    # but
    # if there is only one variant remaining in select_ld_variants
    # take the pair in lowest LD with that variant
    # otherwise method won't run
    
    if (length(select_ld_variants) == 1) {
      to_subset = ld_filtered2 %>%
        dplyr::filter(SNP1 == select_ld_variants) %>%
        dplyr::slice_min(R2) %>%
        pull(SNP2)
      
      to_subset = unique(c(select_ld_variants, to_subset))
      
    } else{
      to_subset = select_ld_variants
    }
    
    ld_matrix = ld_matrix %>%
      dplyr::filter(SNP %in% to_subset) %>%
      dplyr::select(c("SNP", to_subset))
    
    
    # Ensure rows and columns have the same order
    ld_matrix = ld_matrix %>%
      dplyr::arrange(SNP)  # Sort rows
    
    snp_order = ld_matrix$SNP  # Save SNP order
    ld_matrix = ld_matrix %>%
      dplyr::select(SNP, all_of(snp_order)) %>%  # Ensure column order matches row order
      dplyr::rename(variant_id = SNP)
    
    #### retain variants that I selected before based on R2 values:
    
    
    # check presence lead variants
    lead_effect_focus_gene_other_genes$variant_id %in% ld_matrix$variant_id
    
    ### remove variants from effect file that are removed from the LD file
    lead_effect_focus_gene_other_genes = lead_effect_focus_gene_other_genes %>%
      dplyr::arrange(variant_id) %>%
      dplyr::filter(variant_id %in% ld_matrix$variant_id)
    
    
    # checking again I have the same number of variants
    ld_variants = ld_matrix %>%
      dplyr::filter(variant_id %in% lead_effect_focus_gene_other_genes$variant_id) %>%
      dplyr::select(c(
        "variant_id",
        lead_effect_focus_gene_other_genes$variant_id
      ))
  }
  # remove variants in high pairwise LD with all other variants
  
  return(ld_variants)
  
}

#Main
# input_path = "../../data/results/8.5.eQTL_MR/TWMR/input/"
# output_path = "../../data/results/8.5.eQTL_MR/TWMR/output/"
# dir.create(input_path)
dir.create(output_path)
# Example usage
# gene = "ENSG00000000419"
# Ngwas = 239087
# N_eQTLs = 32000
# n_for_LD = 3781

# Step 1: Read input files
# input_data = read_input_files(gene,input_path)

# input_data$eQTL_GWAS_effects
# GENES ENSG00000000419 ENSG00000101126    BETA_GWAS
# 1 rs7268202      0.00000000     -0.02890916  0.002336220
# 2 rs6013040     -0.02673929      0.03259156 -0.005177372
# 3 rs2426214     -0.05108933      0.00000000  0.004598900
# input_data$LDmatrix
# V1        V2       V3
# 1 1.0000000 0.0487694 0.151668
# 2 0.0487694 1.0000000 0.118453
# 3 0.1516680 0.1184530 1.000000
# Step 2: Perform TWMR calculations
# results= run_TWMR_calculations(input_data$beta, input_data$gamma, input_data$C, Ngwas, N_eQTLs, n_for_LD)

# Step 3: Write output file
# write_output_file(gene = gene, results = results, output_path = output_path)

# Print results
# print(results)

# load rsids for ldlink

# extract lead coloc variant
eqtl_coloc_subsets = readr::read_rds(
  file = paste0(
    "../../../OTAR2065_phenotypic_QTLs/data/results/",
    phenotype,
    "/3.coloc/myeQTL_myGWAS_subsets_for_coloc/",
    treatment,
    "_eQTL_subsets.rds"
  )
)
gwas_coloc_subsets = readr::read_rds(
  file = paste0(
    "../../../OTAR2065_phenotypic_QTLs/data/results/",
    phenotype,
    "/3.coloc/myeQTL_myGWAS_subsets_for_coloc/",
    treatment,
    "_my_GWAS_subsets.rds"
  )
)
nominal_res = load_nominal_results(treatment)

coloc_res = readr::read_rds(
  file = paste0(
    "../../../OTAR2065_phenotypic_QTLs/data/results/",
    phenotype,
    "/3.coloc/coloc_results/my_GWAS_my_eQTL/coloc_results_single_causal_variant_500kb.rds"
  )
)

lead_coloc_var = process_coloc_results(treatment, focus_gene, coloc_res = coloc_res)

if (anyNA(lead_coloc_var)) {
  lead_coloc_var_swapped = NA
  selected_gwas_coloc = NA
} else{
  lead_coloc_var_swapped = eqtl_coloc_subsets[[names(eqtl_coloc_subsets)[grepl(focus_gene, names(eqtl_coloc_subsets))]]] %>%
    dplyr::filter(variant_id == lead_coloc_var) %>%
    dplyr::pull(swapped_variant_id)
  selected_gwas_coloc = gwas_coloc_subsets[[names(gwas_coloc_subsets)[grepl(focus_gene, names(gwas_coloc_subsets))]]] %>%
    dplyr::filter(variant_id == lead_coloc_var) %>%
    dplyr::select(variant_id, beta)
  
}


# loading lead eQTL for this focus_gene, and other eGenes for lead and colocalised eQTL

if (anyNA(lead_coloc_var)) {
  coloc_effect_all_genes = NA
  
  
} else{
  coloc_effect_all_genes = nominal_res %>%
    dplyr::select(variant_id, gene, pval_nominal, slope) %>%
    dplyr::filter(variant_id %in% lead_coloc_var |
                    variant_id %in% lead_coloc_var_swapped)
  
  # check beta swap, change slope sign if true
  if (coloc_effect_all_genes$variant_id[1] == lead_coloc_var_swapped) {
    message("There is a ALT allee swap with respect to coloc lead variant: swapping slope...")
    coloc_effect_all_genes = coloc_effect_all_genes %>%
      dplyr::mutate(slope = -slope)
  }
  
}

# IMPORTANT step: check lead eQTLs for the eGene
effects = extract_relevant_variants_eQTL(
  focus_gene = focus_gene,
  nominal_res = nominal_res,
  coloc_effect_all_genes = coloc_effect_all_genes
)


if (nrow(effects$lead_effect_focus_gene) == 0) {
  message("No eQTL effects for gene ", focus_gene , " in treatment ", treatment)
  
  results = NA
  readr::write_rds(results, file = paste0(
    paste0(output_path, treatment, "/", phenotype, "/"),
    focus_gene,
    ".rds"
  ))
} else{
  
  # if at least two effects are identical, these will be in high LD for sure
  # remove those variants
  if (nrow(effects$lead_effect_focus_gene) != 1) {
    var_to_remove = c()
    
    for (n in 1:(nrow(effects$lead_effect_focus_gene)-1)) {
      if(effects$lead_effect_focus_gene[n,"slope"] == effects$lead_effect_focus_gene[n + 1,"slope"]){
        var_to_remove = append(var_to_remove,values = effects$lead_effect_focus_gene[n + 1,"variant_id"])
      }
    }
    var_to_remove = unlist(var_to_remove)
    
  }
  
  if(exists("var_to_remove")){
    if(!is_empty(var_to_remove)){
    effects$lead_effect_focus_gene = effects$lead_effect_focus_gene %>%
      dplyr::filter(!variant_id %in% var_to_remove)
    
    effects$lead_effect_focus_gene_other_genes = effects$lead_effect_focus_gene_other_genes %>%
      dplyr::filter(!variant_id %in% var_to_remove)
    }
  }
  ######## start of LD calculations ########
  
  
  ld_variants = load_and_process_ld_file(
    coloc_effect_all_genes,
    effects$lead_effect_focus_gene_other_genes,
    effects$lead_effect_focus_gene
  )
  ######## end of LD calculations ##########
  
  
  if (anyNA(ld_variants)) {
    message("No variants found in LD file for gene ", focus_gene , " in treatment ", treatment)
    
    results = NA
    readr::write_rds(results, file = paste0(
      paste0(output_path, treatment, "/", phenotype, "/"),
      focus_gene,
      ".rds"
    ))
  }else{
    # extract only variants in low LD from the beta matrix
    effects$lead_effect_focus_gene_other_genes = effects$lead_effect_focus_gene_other_genes %>%
      dplyr::filter(variant_id %in% ld_variants$variant_id)
    
    
    
    
    # extract GWAS beta for all variants
    
    pheno_gwas = readr::read_csv(
      file = paste0(
        "../../../OTAR2065_phenotypic_QTLs/data/results/",
        phenotype,
        "/2.check_association_results/lm_1pct_filtered_deflated/all_res_pqtl_",
        phenotype,
        ".csv"
      )
    )
    
    filtered_pheno_gwas = pheno_gwas %>%
      dplyr::filter(snp %in% ld_variants$variant_id) %>%
      dplyr::select(contains(treatment), snp) %>%
      dplyr::select(contains("coef"), snp) %>%
      dplyr::rename(variant_id = snp,
                    BETA_GWAS = paste0("coef_", treatment)) %>%
      dplyr::mutate(BETA_GWAS = tidyr::replace_na(BETA_GWAS, 0))#  NA variant effects are 0
    
    
    ######### putting together beta table and ld table and testing ########
    
    
    Ngwas = 261
    N_eQTLs = 261
    n_for_LD = 13160 # info from TopMED readme
    
    
    input_data = list(
      "eQTL_GWAS_effects" = effects$lead_effect_focus_gene_other_genes %>%
        dplyr::left_join(filtered_pheno_gwas) %>%
        dplyr::rename(GENES = variant_id) %>%
        dplyr::mutate(BETA_GWAS = tidyr::replace_na(BETA_GWAS, 0)) %>% #  NA variant effects are 0
        as.data.frame(),
      "LDmatrix" = ld_variants %>%
        dplyr::select(-variant_id) %>%
        as.data.frame()
    )
    
    ####### IMPORTANT ########
    # the focus gene needs to be in the second column of eQTL_GWAS_effects (immediately after "GENES" )
    # otherwise the results won't be correct
    # the authors of TWMR shouldn't have made the code so that the selection is positional IMHO
    # so should potentially fix that
    
    column_order = c("GENES",focus_gene,setdiff(colnames(effects$lead_effect_focus_gene_other_genes),c(focus_gene,"variant_id")),"BETA_GWAS")
    
    input_data$eQTL_GWAS_effects = input_data$eQTL_GWAS_effects %>%
      dplyr::select(all_of(column_order)) %>%
      as.data.frame()
    
    if(which(colnames(input_data$eQTL_GWAS_effects) %in% focus_gene) !=2){
      message("ERROR: the focus gene is not in second position in the effect matrix for TWMR. Exiting ...")
      break()
    }
    
    ###################
    
    input_data[["beta"]] = as.matrix(input_data$eQTL_GWAS_effects[, 2:(ncol(input_data$eQTL_GWAS_effects) - 1)])
    
    # extracting beta for GWAS, here called "gamma" , G in paper
    input_data[["gamma"]] = as.matrix(input_data$eQTL_GWAS_effects[, ncol(input_data$eQTL_GWAS_effects)])
    
    # loading LD matrix for SNP effects on the tested gene, here called "C". Diagonal left to right is always 1 (LD with itself)
    input_data[["C"]] = as.matrix(input_data$LDmatrix[, 1:ncol(input_data$LDmatrix)])
    
    # if there is only one variant, this won't work
    # eg. focus_gene = "SDSL"
    # treatment = "LPS"
    
    # Step 2: Perform TWMR calculations
    results = run_TWMR_calculations(input_data$beta,
                                    input_data$gamma,
                                    input_data$C,
                                    Ngwas,
                                    N_eQTLs,
                                    n_for_LD)
    
    # what if I use only the top coloc variant
    
    # input_data[["beta"]] = as.matrix(input_data$eQTL_GWAS_effects[, 2:(ncol(input_data$eQTL_GWAS_effects) - 1)])[1,]
    #
    # # extracting beta for GWAS, here called "gamma" , G in paper
    # input_data[["gamma"]] = as.matrix(input_data$eQTL_GWAS_effects[, ncol(input_data$eQTL_GWAS_effects)])[1,]
    #
    # # loading LD matrix for SNP effects on the tested gene, here called "C". Diagonal left to right is always 1 (LD with itself)
    # input_data[["C"]]= as.matrix(input_data$LDmatrix[, 1:ncol(input_data$LDmatrix)])[1,1]
    # results= run_TWMR_calculations(input_data$beta, input_data$gamma, input_data$C, Ngwas, N_eQTLs, n_for_LD) # does not work with one variant
    
    # Step 3: Write output file
    dir.create(paste0(output_path, treatment, "/", phenotype, "/"),
               recursive = TRUE)
    
    results$LD_matrix = input_data$LDmatrix
    results$eQTL_GWAS_effects = input_data$eQTL_GWAS_effects
    readr::write_rds(results, file = paste0(
      paste0(output_path, treatment, "/", phenotype, "/"),
      focus_gene,
      ".rds"
    ))
    
    
    scales::scientific(results$P) # "3.6e-14" for LRRK2, 4.3 e-14 when removing SNPs in high LD first, alpha barely changes
    # 5.6e-12 when trimming ld code
    # 9.96e-01 for SLC2A13, near LRRK2 and with many eQTLs in the matrix from that gene
    # TREM2 untreated p-val 0.99
    
    
    #### check distribution of TREM2 p-vals with LD and betas
    #
    # input_data_list <- map(1:30, function(i) {
    #   factor <- 1 + (i - 1) * 0.1  # Increasing factor starting at 1
    #
    #   beta_modified <- input_data$beta
    #   beta_modified[, 2] <- beta_modified[, 2] * factor  # Modify column for "TREM2"
    #
    #   list(
    #     beta = beta_modified,
    #     gamma =input_data$gamma,
    #     C =input_data$C # Identity matrix for C
    #   )
    # })
    # # Use pmap to iterate over each set of inputs
    # # results <- pmap(input_data_list, ~ run_TWMR_calculations(..1, ..2, ..3, Ngwas, N_eQTLs, n_for_LD))
    # results2= run_TWMR_calculations(input_data_list[[4]]$beta, input_data_list[[4]]$gamma, input_data_list[[4]]$C, Ngwas, N_eQTLs, n_for_LD)
    # results3= run_TWMR_calculations(input_data_list[[8]]$beta, input_data_list[[8]]$gamma, input_data_list[[8]]$C, Ngwas, N_eQTLs, n_for_LD)
    #
    #
    #
    # ### checking correlations of genes - the authors remove genes
    # # that show  r2 ≥ 0.4. The correlation r2 was estimated as Pearson’s correlation between the Z scores of the shared, independent eQTLs.
    # cor(input_data$beta[,3],input_data$beta[,2],)^2
    #
    # cor_to_remove = cor(input_data$beta, use = "pairwise.complete.obs") %>%
    #     as.data.frame() %>%
    #     rownames_to_column(var = "Var1") %>%
    #   dplyr::select(TREM2,Var1)
    #
    #
    # print(cor_long)
    
    # Print results
    print(results)
    
  }
}

