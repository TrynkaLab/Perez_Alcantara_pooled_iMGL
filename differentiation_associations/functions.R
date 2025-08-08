# functions

load_sc_donors = function( pools = paste0("pool",2:17) ,
                           directory =  "../../../OTAR2065_scRNA_data_processing/data/"){
  lines_in_pools = read.table(paste0(directory,"lines_in_pools.txt"), 
                              header = TRUE)
  #lines_in_pools has three columns: Pool, N_lines, Lines (last one has lines separated by ;)
  sc_ncells = list()
  for (pool in pools) {
    message("Working on pool ", pool,"...")
    files = list.dirs(
      paste0(
        directory,"cellranger/",
        pool
      ),
      recursive = FALSE,
      full.names = FALSE
    )
    for (file in files) {
      message("... file ", file,"...")
      
      sc_ncells[[paste0(pool,".",file)]]  = read.table(paste0("../../../OTAR2065_scRNA_data_processing/data/cellranger/",
                                                              pool,"/",file,"/vireoOutput/donor_ids.tsv"),
                                                       header = TRUE)
      
      # if doublets, put doublet in donor_id
      # if unassigned, check which probability is higher
      sc_ncells[[paste0(pool,".",file)]] = sc_ncells[[paste0(pool,".",file)]] %>%
        dplyr::mutate(new_donor_id = case_when(donor_id == "doublet" ~ best_doublet,
                                               donor_id == "unassigned" & prob_max < prob_doublet ~ best_doublet,
                                               donor_id == "unassigned" & prob_max >= prob_doublet ~ best_singlet,
                                               .default = donor_id)) %>%
        tidyr::separate_rows(new_donor_id, sep = ",") %>%
        distinct()  %>%
        dplyr::count(new_donor_id) %>%
        dplyr::rename(Line = new_donor_id, Frequency = n)
      
      
      
      
      lines = unlist(strsplit(lines_in_pools[lines_in_pools$Pool == pool,"Lines"],split = ";"))
      sc_ncells[[paste0(pool, ".", file)]] =  sc_ncells[[paste0(pool, ".", file)]] %>%
        dplyr::rows_insert(tibble(Line = lines[!lines %in% sc_ncells[[paste0(pool, ".", file)]]$Line])) %>%
        tidyr::replace_na(list(Frequency = 0))
      
      sc_ncells[[paste0(pool,".",file)]]  = sc_ncells[[paste0(pool,".",file)]] %>%
        dplyr::reframe(Line = Line,
                       prop = Frequency / sum(Frequency),
                       pool = pool, 
                       sample = file, 
                       stage = "microglia")
      
      
    }
  }
  
  sc_ncells = do.call(rbind,Map(cbind,set = names(sc_ncells),sc_ncells))
  return(sc_ncells)
}

alluvial_plot = function(long_df){
  p = ggplot(long_df, aes(y = prop,x = sample,
                          alluvium=Line,stratum = Line,
                          fill = Line, colour = Line,
                          label = Line)) +
    scale_x_discrete(expand = c(.1, .1)) +
    scale_fill_viridis_d(option = "D")+
    scale_color_viridis_d(option = "D")+
    geom_flow(alpha = .6) +
    geom_stratum(alpha = .0) +
    # geom_text(colour = "black", position = position_stack(vjust = 0.5),size = 3) +
    ggfittext::geom_fit_text(stat = "stratum",colour = "black", width = 1/4, min.size = 3) +
    theme_bw() + 
    theme(legend.position = "none",axis.text= element_text(face="bold"),
          axis.title = element_text(face="bold"), plot.title =element_text(face="bold") ) +
    labs(title = paste0("Changing donor proportions for ",unique(long_df$pool))) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  plot(p)
}
# Find the closest value in Y to a given value in X
find_closest_value = function(x, y) {
  y[which.min(abs(y - x))]
}

# scaled microglia vs preMAC proliferations
prolif_microglia_premac = function(data_with_error){
  efficiency = list()
  for(p in paste0("pool",2:17)){
    if(p == "pool3"){
      efficiency[[p]] = data_with_error %>%
        dplyr::filter(stage!="iPSC" & pool=="pool3") %>%
        tidyr::pivot_wider(names_from = sample,values_from = prop_adjusted_mean,id_cols = Line) %>%
        dplyr::reframe(line=Line,
                       sc_untreated = log1p(h3_diff2_untreated) / log1p(d35_preMAC),
                       sc_IFN =  log1p(h3_diff2_IFNg) / log1p(d35_preMAC),
                       sc_LPS =  log1p(h3_diff2_LPS) / log1p(d35_preMAC),
                       phago_untreated =  log1p(`phago_s3D-1`) / log1p(d57_preMAC),
                       phago_IFN =  log1p(`phago_s3E-1`) / log1p(d57_preMAC),
                       phago_LPS =  log1p(`phago_s3F-1`) / log1p(d57_preMAC)) %>%
        # Puigdevall et al prolif. rate: not much different distribution from raw props.
        
        tidyr::pivot_longer(cols =-line,names_to = "sample",values_to = "scaled_proportion") %>%
        dplyr::mutate(pool="pool3")
      # don't have right preMAC for migration
      # don't have right preMACs for pool2
    }
    if(p == "pool4"){
      efficiency[[p]] = data_with_error %>%
        dplyr::filter(stage!="iPSC" & pool=="pool4") %>%
        tidyr::pivot_wider(names_from = sample,values_from = prop,id_cols = Line) %>%
        dplyr::reframe(line=Line,
                       sc_untreated = log1p(h4_diff2_untreated) / log1p(d35_preMAC),
                       sc_IFN =  log1p(h4_diff2_IFNg) / log1p(d35_preMAC),
                       sc_LPS =  log1p(h4_diff2_LPS) / log1p(d35_preMAC),
                       phago_untreated =  log1p(`phago_s4D-1`) / log1p(d57_preMAC),
                       phago_IFN =  log1p(`phago_s4E-1`) / log1p(d57_preMAC),
                       phago_LPS =  log1p(`phago_s4F-1`) / log1p(d57_preMAC)) %>%
        tidyr::pivot_longer(cols =-line,names_to = "sample",values_to = "scaled_proportion") %>%
        dplyr::mutate(pool="pool4")
      # don't have right preMAC for migration
      
    }
    if(p == "pool5"){
      efficiency[[p]] = data_with_error %>%
        dplyr::filter(stage!="iPSC" & pool=="pool5") %>%
        tidyr::pivot_wider(names_from = sample,values_from = prop_adjusted_mean,id_cols = Line) %>%
        dplyr::reframe(line=Line,
                       sc_untreated = log1p(P5_diff2_untreated) / log1p(d35_preMAC),
                       sc_IFN =  log1p(P5_diff2_IFN) / log1p(d35_preMAC),
                       sc_LPS =  log1p(P5_diff2_LPS) / log1p(d35_preMAC)) %>%
        tidyr::pivot_longer(cols =-line,names_to = "sample",values_to = "scaled_proportion") %>%
        dplyr::mutate(pool="pool5")
      # don't have seeded cells for migration
      # don't have right preMAC for phagocytosis
      
    }
    if(p == "pool6"){
      efficiency[[p]] = data_with_error %>%
        dplyr::filter(stage!="iPSC" & pool=="pool6") %>%
        tidyr::pivot_wider(names_from = sample,values_from = prop_adjusted_mean,id_cols = Line) %>%
        dplyr::reframe(line=Line,
                       sc_untreated = log1p(P6_diff2_untreated) / log1p(d35_preMAC),
                       sc_IFN =  log1p(P6_diff2_IFN) / log1p(d35_preMAC),
                       sc_LPS =  log1p(P6_diff2_LPS) / log1p(d35_preMAC)) %>%
        tidyr::pivot_longer(cols =-line,names_to = "sample",values_to = "scaled_proportion") %>%
        dplyr::mutate(pool="pool6")
      # don't have seeded cells for migration
      # don't have right preMAC for phagocytosis
      
    }
    if(p == "pool7"){
      efficiency[[p]] = data_with_error %>%
        dplyr::filter(stage!="iPSC" & pool=="pool7") %>%
        tidyr::pivot_wider(names_from = sample,values_from = prop_adjusted_mean,id_cols = Line) %>%
        dplyr::reframe(line=Line,
                       sc_untreated = log1p(P7_diff2_untreated) / log1p(D36_preMAC),
                       sc_IFN =  log1p(P7_diff2_IFN) / log1p(D36_preMAC),
                       sc_LPS =  log1p(P7_diff2_LPS) / log1p(D36_preMAC),
                       migr_untreated =  log1p(Mig17_Untr) / log1p(D47_preMAC),
                       migr_LPS =  log1p(Mig17_LPS) / log1p(D47_preMAC),
                       phago_untreated =  log1p(`phago_s7A-1`) / log1p(D54_preMAC),
                       phago_IFN =  log1p(`phago_s7B-1`) / log1p(D54_preMAC),
                       phago_LPS =  log1p(`phago_s7C-1`) / log1p(D54_preMAC)) %>%
        tidyr::pivot_longer(cols =-line,names_to = "sample",values_to = "scaled_proportion") %>%
        dplyr::mutate(pool="pool7")
    }
    if(p == "pool8"){
      efficiency[[p]] = data_with_error %>%
        dplyr::filter(stage!="iPSC" & pool=="pool8") %>%
        tidyr::pivot_wider(names_from = sample,values_from = prop_adjusted_mean,id_cols = Line) %>%
        dplyr::reframe(line=Line,
                       sc_untreated = log1p(P8_diff2_untreated) / log1p(D36_preMAC),
                       sc_IFN =  log1p(P8_diff2_IFN) / log1p(D36_preMAC),
                       sc_LPS =  log1p(P8_diff2_LPS) / log1p(D36_preMAC),
                       migr_untreated =  log1p(Mig16_Untr) / log1p(D40_preMAC),
                       migr_LPS =  log1p(Mig16_LPS) / log1p(D40_preMAC),
                       phago_untreated =  log1p(`phago_s8A-1`) / log1p(D54_preMAC),
                       phago_IFN =  log1p(`phago_s8A-1`) / log1p(D54_preMAC),
                       phago_LPS =  log1p(`phago_s8A-1`) / log1p(D54_preMAC)) %>%
        tidyr::pivot_longer(cols =-line,names_to = "sample",values_to = "scaled_proportion") %>%
        dplyr::mutate(pool="pool8")
    }
    if(p == "pool9"){
      efficiency[[p]] = data_with_error %>%
        dplyr::filter(stage!="iPSC" & pool=="pool9") %>%
        tidyr::pivot_wider(names_from = sample,values_from = prop_adjusted_mean,id_cols = Line) %>%
        dplyr::reframe(line=Line,
                       sc_untreated_A = log1p(P9_Diff2_Untreated_A) / log1p(D35_PreMac),
                       sc_untreated_B = log1p(P9_Diff2_Untreated_B) / log1p(D35_PreMac),
                       sc_untreated_C = log1p(P9_Diff2_Untreated_B) / log1p(D35_PreMac),
                       sc_untreated_D = log1p(P9_Diff2_Untreated_B) / log1p(D35_PreMac),
                       sc_IFN_A =  log1p(P9_Diff2_IFN_A) / log1p(D35_PreMac),
                       sc_IFN_B =  log1p(P9_Diff2_IFN_B) / log1p(D35_PreMac),
                       sc_IFN_C =  log1p(P9_Diff2_IFN_C) / log1p(D35_PreMac),
                       sc_LPS_A =  log1p(P9_Diff2_LPS_A) / log1p(D35_PreMac),
                       sc_LPS_B =  log1p(P9_Diff2_LPS_B) / log1p(D35_PreMac),
                       sc_LPS_C =  log1p(P9_Diff2_LPS_C) / log1p(D35_PreMac),
                       migr_untreated =  log1p(Untr_cells) / log1p(D39_PreMac),
                       migr_LPS =  log1p(LPS_cells) / log1p(D39_PreMac),
                       phago_untreated =  log1p(`phago_s9A-1`) / log1p(D50_PreMac),
                       phago_LPS =  log1p(`phago_s9C-1`) / log1p(D50_PreMac)) %>%
        tidyr::pivot_longer(cols =-line,names_to = "sample",values_to = "scaled_proportion") %>%
        dplyr::mutate(pool="pool9")
    }
    if(p == "pool10"){
      efficiency[[p]] = data_with_error %>%
        dplyr::filter(stage!="iPSC" & pool=="pool10") %>%
        tidyr::pivot_wider(names_from = sample,values_from = prop_adjusted_mean,id_cols = Line) %>%
        dplyr::reframe(line=Line,
                       sc_untreated_A = log1p(P10_Diff2_Untreated_A) / log1p(D36_PreMac),
                       sc_untreated_B = log1p(P10_Diff2_Untreated_B) / log1p(D36_PreMac),
                       sc_untreated_C = log1p(P10_Diff2_Untreated_B) / log1p(D36_PreMac),
                       sc_untreated_D = log1p(P10_Diff2_Untreated_B) / log1p(D36_PreMac),
                       sc_untreated_old_A = log1p(P10_Diff5_Untreated_A) / log1p(D54_PreMac),
                       sc_untreated_old_B = log1p(P10_Diff5_Untreated_B) / log1p(D54_PreMac),
                       sc_IFN_A =  log1p(P10_Diff2_IFN_A) / log1p(D36_PreMac),
                       sc_IFN_B =  log1p(P10_Diff2_IFN_B) / log1p(D36_PreMac),
                       sc_IFN_C =  log1p(P10_Diff2_IFN_C) / log1p(D36_PreMac),
                       sc_IFN_old = log1p(P10_Diff5_IFN) / log1p(D54_PreMac),
                       sc_LPS_A =  log1p(P10_Diff2_LPS_A) / log1p(D36_PreMac),
                       sc_LPS_B =  log1p(P10_Diff2_LPS_B) / log1p(D36_PreMac),
                       sc_LPS_C =  log1p(P10_Diff2_LPS_C) / log1p(D36_PreMac),
                       sc_LPS_old = log1p(P10_Diff5_LPS) / log1p(D54_PreMac),
                       migr_untreated =  log1p(Untr_cell) / log1p(D43_PreMac),
                       migr_LPS =  log1p(LPS_cells) / log1p(D43_PreMac),
                       phago_untreated =  log1p(`phago_s10A-1`) / log1p(D50_PreMac),
                       phago_IFN =  log1p(`phago_s10B-1`) / log1p(D50_PreMac),
                       phago_LPS =  log1p(`phago_s10C-1`) / log1p(D50_PreMac)) %>%
        tidyr::pivot_longer(cols =-line,names_to = "sample",values_to = "scaled_proportion") %>%
        dplyr::mutate(pool="pool10")
    }
    if(p == "pool11"){
      efficiency[[p]] = data_with_error %>%
        dplyr::filter(stage!="iPSC" & pool=="pool11") %>%
        tidyr::pivot_wider(names_from = sample,values_from = prop_adjusted_mean,id_cols = Line) %>%
        dplyr::reframe(line=Line,
                       sc_untreated_A = log1p(P11_Diff2_Untreated_A) / log1p(day35_preMAC),
                       sc_untreated_B = log1p(P11_Diff2_Untreated_B) / log1p(day35_preMAC),
                       sc_IFN =  log1p(P11_Diff2_IFN) / log1p(day35_preMAC),
                       sc_LPS =  log1p(P11_Diff2_LPS) / log1p(day35_preMAC),
                       migr_untreated =  log1p(pool11_UNTR) / log1p(day39_preMAC),
                       migr_LPS =  log1p(pool11_LPS) / log1p(day39_preMAC),
                       phago_untreated =  log1p(pool11_A1) / log1p(day49_preMAC),
                       phago_IFN =  log1p(pool11_B1) / log1p(day49_preMAC),
                       phago_LPS =  log1p(pool11_C1) / log1p(day49_preMAC)) %>%
        tidyr::pivot_longer(cols =-line,names_to = "sample",values_to = "scaled_proportion") %>%
        dplyr::mutate(pool="pool11")
    }
    if(p == "pool13"){
      efficiency[[p]] = data_with_error %>%
        dplyr::filter(stage!="iPSC" & pool=="pool13") %>%
        tidyr::pivot_wider(names_from = sample,values_from = prop_adjusted_mean,id_cols = Line) %>%
        dplyr::reframe(line=Line,
                       sc_untreated_A = log1p(P13_Diff2_Untreated_A) / log1p(day36_preMAC),
                       sc_untreated_B = log1p(P13_Diff2_Untreated_B) / log1p(day36_preMAC),
                       sc_IFN =  log1p(P13_Diff2_IFN) / log1p(day36_preMAC),
                       sc_LPS =  log1p(P13_Diff2_LPS) / log1p(day36_preMAC),
                       migr_untreated =  log1p(pool13_UNTR) / log1p(day43_preMAC),
                       migr_IFN =  log1p(pool13_IFN) / log1p(day43_preMAC),
                       migr_LPS =  log1p(pool13_LPS) / log1p(day43_preMAC),
                       phago_untreated =  log1p(pool13_A1) / log1p(day50_preMAC),
                       phago_IFN =  log1p(pool13_B1) / log1p(day50_preMAC),
                       phago_LPS =  log1p(pool13_C1) / log1p(day50_preMAC)) %>%
        tidyr::pivot_longer(cols =-line,names_to = "sample",values_to = "scaled_proportion") %>%
        dplyr::mutate(pool="pool13")
    }
    if(p == "pool14"){
      efficiency[[p]] = data_with_error %>%
        dplyr::filter(stage!="iPSC" & pool=="pool14") %>%
        tidyr::pivot_wider(names_from = sample,values_from = prop_adjusted_mean,id_cols = Line) %>%
        dplyr::reframe(line=Line,
                       sc_untreated_A = log1p(P14_Diff2_Untreated_A) / log1p(D36_PreMac),
                       sc_untreated_B = log1p(P14_Diff2_Untreated_B) / log1p(D36_PreMac),
                       sc_IFN_A =  log1p(P14_Diff2_IFN_A) / log1p(D36_PreMac),
                       sc_IFN_B =  log1p(P14_Diff2_IFN_B) / log1p(D36_PreMac),
                       sc_LPS_A =  log1p(P14_Diff2_LPS_A) / log1p(D36_PreMac),
                       sc_LPS_B =  log1p(P14_Diff2_LPS_B) / log1p(D36_PreMac),
                       migr_untreated =  log1p(Mig24_Untr) / log1p(D39_PreMac),
                       migr_LPS =  log1p(Mig24_LPS) / log1p(D39_PreMac)) %>%
        tidyr::pivot_longer(cols =-line,names_to = "sample",values_to = "scaled_proportion") %>%
        dplyr::mutate(pool="pool14")
    }
    if(p == "pool15"){
      efficiency[[p]] = data_with_error %>%
        dplyr::filter(stage!="iPSC" & pool=="pool15") %>%
        tidyr::pivot_wider(names_from = sample,values_from = prop_adjusted_mean,id_cols = Line) %>%
        dplyr::reframe(line=Line,
                       sc_untreated_A = log1p(P15_Diff2_Untreated_A) / log1p(D36_PreMac),
                       sc_untreated_B = log1p(P15_Diff2_Untreated_B) / log1p(D36_PreMac),
                       sc_IFN_A =  log1p(P15_Diff2_IFN_A) / log1p(D36_PreMac),
                       sc_IFN_B =  log1p(P15_Diff2_IFN_B) / log1p(D36_PreMac),
                       sc_LPS_A =  log1p(P15_Diff2_LPS_A) / log1p(D36_PreMac),
                       sc_LPS_B =  log1p(P15_Diff2_LPS_B) / log1p(D36_PreMac),
                       migr_untreated =  log1p(Mig23_Untr) / log1p(D39_PreMac),
                       migr_IFN =  log1p(Mig23_IFN) / log1p(D39_PreMac),
                       migr_LPS =  log1p(Mig23_LPS) / log1p(D39_PreMac)) %>%
        tidyr::pivot_longer(cols =-line,names_to = "sample",values_to = "scaled_proportion") %>%
        dplyr::mutate(pool="pool15")
    }
    
    
  }
  efficiency = do.call("rbind",efficiency)
  return(efficiency)
}

fix_burden_matrices <- function(df) {
  df %>%
    tibble::as_tibble(rownames = "gene") %>%
    dplyr::rename_with(~str_extract(., "(?<=-)\\w+"), .cols = contains("H"))
}

# subset column names to shared lines, including "gene"
subset_burden_shared_lines = function(x,y) { 
  x %>%
    dplyr::select(c("gene",sort(lubridate::intersect(colnames(.),y)))) 
}

filter_rare_variant_genes = function(tb,fraction) { 
  message("The fraction corresponds to ",floor((ncol(tb)-1)*fraction), " lines")
  tb %>%
    tibble::column_to_rownames(var="gene") %>%
    dplyr::filter(rowSums(. > 0) >= floor(ncol(.)*fraction)) %>%
    dplyr::filter(rowSums(. == 0) >= floor(ncol(.)*fraction)) %>%
    tibble::rownames_to_column(var="gene")
    
}

rare_variant_to_long = function(tbl_list) { 
  # returns a tibble
  for(i in names(tbl_list)){
    tbl_list[[i]] = tbl_list[[i]] %>%
      tidyr::pivot_longer(.,cols = !gene,names_to = "line",values_to = "rare_burden") %>%
      dplyr::mutate(rare_mutation_type = i)
  }
  return(do.call("rbind",tbl_list))
  
}

ggcheck_the_model <- function(m1){
  gg1 <- ggcheck_the_qq(m1)
  gg2 <- ggcheck_the_spreadlevel(m1)
  cowplot::plot_grid(gg1, gg2, nrow = 1)
}

ggcheck_the_qq = function(m1,
                          line = "robust",
                          n_boot = 200){
  n <- nobs(m1)
  m1_res <- residuals(m1)
  #sigma_m1_res <- sigma(m1)
  
  normal_qq <- ppoints(n) %>%
    qnorm()
  sample_qq <- m1_res[order(m1_res)]
  
  # mean + sd
  parametric_slope <- sd(sample_qq)
  parametric_intercept <- mean(sample_qq)
  
  # quartiles
  m1_quartiles <- quantile(m1_res, c(.25, .75))
  qnorm_quartiles <- qnorm( c(.25, .75))
  m1_diff <- m1_quartiles[2] - m1_quartiles[1]
  qnorm_diff <- qnorm_quartiles[2] - qnorm_quartiles[1] # = 1.349
  quartile_slope <- m1_diff/qnorm_diff
  quartile_intercept <- median(m1_quartiles) # median of quartiles not quantiles
  
  # robust uses MASS:rlm (default arguments?)
  qq_observed <- data.table(normal_qq = normal_qq,
                            sample_qq = sample_qq)
  m2 <- rlm(sample_qq ~ normal_qq, data = qq_observed)
  robust_intercept <- coef(m2)[1]
  robust_slope <- coef(m2)[2]
  
  # re-sample ribbon
  set.seed(1)
  resample_qq_model <- numeric(n_boot*n)
  Y <- simulate(m1, n_boot)
  fd <- model.frame(m1) %>%
    data.table
  inc <- 1:n
  for(sim_i in 1:n_boot){
    # parametric bound
    fd[, (1) := Y[,sim_i]]
    m1_class <- class(m1)[1]
    if(m1_class == "lm"){
      ff <- lm(formula(m1), data = fd) 
    }
    if(m1_class == "lmerModLmerTest" | m1_class == "lmerMod"){
      ff <- lmer(formula(m1), data = fd)
    }
    y_res <- residuals(ff)
    resample_qq <- y_res[order(y_res)]
    resample_qq_model[inc] <- resample_qq
    inc <- inc + n
    
    # robust bound
    qq_resampled <- data.table(normal_qq = normal_qq,
                               resample_qq = resample_qq)
    m2_resample <- rlm(resample_qq ~ normal_qq, data = qq_resampled)
    
  }
  
  qq_sim <- data.table(normal_qq = normal_qq,
                       resample_qq_model = resample_qq_model)
  
  qq_ci_model <- qq_sim[, .(median = median(resample_qq_model),
                            lower = quantile(resample_qq_model, 0.025),
                            upper = quantile(resample_qq_model, 0.975)),
                        by = normal_qq]
  m2_boot <- rlm(median ~ normal_qq, data = qq_ci_model)
  robust_intercept_boot <- coef(m2_boot)[1]
  robust_slope_boot <- coef(m2_boot)[2]
  
  ggplot(data = qq_observed,
         aes(x = normal_qq, y = sample_qq)) +
    
    # ribbon
    geom_ribbon(data = qq_ci_model,
                aes(ymin = lower,
                    ymax = upper,
                    y = median,
                    fill = "band"),
                fill = "gray",
                alpha = 0.6) +
    # draw points
    geom_point() +
    
    # robust
    geom_abline(aes(intercept = robust_intercept,
                    slope = robust_slope,
                    color = "robust"),
                show.legend = FALSE,
                size = 0.75) +

    xlab("Normal Quantiles") +
    ylab("Sample Quantiles") +
    
    scale_color_manual(values = pal_okabe_ito[c(1:2,5:6)]) +
    theme_minimal_grid() +
    NULL
  
}

# manual scaling (and centering to the mean)

scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}


## ----ggcheck_the_spreadlevel----------------------------------------------------
ggcheck_the_spreadlevel <- function(m1,
                                    n_boot = 200){
  n <- nobs(m1)
  m1_res <- residuals(m1)
  m1_scaled <- m1_res/sd(m1_res)
  m1_root <- sqrt(abs(m1_scaled))
  m1_fitted <- fitted(m1)
  
  m2 <- lm(m1_root ~ m1_fitted)
  m2_intercept <- coef(m2)[1]
  m2_slope <- coef(m2)[2]
  
  plot_data <- data.table(
    m1_res = sqrt(abs(m1_scaled)),
    fitted = m1_fitted
  )
  
  ggplot(data = plot_data,
         aes(x = fitted, y = m1_res)) +

    geom_point() +
    
    geom_smooth(method = lm) +
   
  xlab("Fitted") +
    ylab("root abs-scaled-residual") +
    
    
    scale_color_manual(values = pal_okabe_ito[c(1:2,5:6)]) +
    theme_minimal_grid() +
    NULL
}


extract_SKATO_results <- function(dlist) {
  dlist = dlist[[1]] ### careful here
  
  df = data.frame("p_val"= dlist$p.value, 
                  "gene_name" = dlist$gene_name,
                  "gene_id"=dlist$gene_id)
  
  if(!is.null(dlist$resampling_pval)){
    
    # calculating proportion of permuted pval < actual pval
    # a.k.a permuted p-val
    df$resampling_pval = dlist$resampling_pval
  }else{
    df$resampling_pval = NA
  }
  return(df)
}


extract_burden_results <- function(dlist) {

  df = data.frame("p_val"= dlist$coefficients["rare_burden","Pr(>|t|)"], 
                  "coef"= dlist$coefficients["rare_burden","Estimate"],
                  "gene_name" = dlist$gene)
  
  return(df)
}


## qqplot
gg_qqplot = function(pvals) {
  df = data.frame(
    observed = -log10(sort(pvals)),
    expected = -log10(ppoints( length(pvals)))
  )
  ggplot(df) +
    geom_point(aes(expected, observed), shape = 1, size = 3, col = "grey40") +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    theme_bw() + 
    xlab( expression(paste("Expected -log"[10], plain(P)))) +
    ylab(expression(paste("Observed -log"[10], plain(P))))
}

