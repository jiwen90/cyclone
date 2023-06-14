
save_heatmap_png <- function(x, filename) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   png(filename, width=2300, height=2000, units="px", res=500)
   ComplexHeatmap::draw(x)
   dev.off()
}

save_heatmap_pdf <- function(x, filename) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=7, height=6)
   ComplexHeatmap::draw(x)
   dev.off()
}

paired_colset <- RColorBrewer::brewer.pal(12,"Paired")
# remove 6th and 8th
paired_colz <- paired_colset[-c(6,8, 12)]

selected_colz <- c("#438EC0", "#63A8A0", "#7D54A5", "#A6CEE3", "#3BA432", "#f16689")

auto_colset <- function(n=2, colz="paired"){
  if (colz == "paired") {
     ret <- grDevices::colorRampPalette(paired_colz)(n)
  } else {
     ret <- grDevices::colorRampPalette(selected_colz)(n)
  }
  return(ret)
}

cycling_cutoff <- 0.12
fc_filter_plotting <- 0.02

fix_pd1 <- function(x) {
  # replace aPD1 with aPD-1 in vector
  x <- gsub("aPD1", "aPD-1", x)
  return(x)
}

temp3 <- function(x)
{
  if (grepl("0", x, fixed=TRUE))
  {
    return (NA)
  }
  return (x)
}

compute_trb_clones_nt <- function(m)
{
  m$TRB_chain_cdr3_nt[is.na(m$TRB_chain_cdr3_nt)] <- 0
  m <- m %>% mutate(TRB_nt_dataset_clonotype_id = lapply(TRB_chain_cdr3_nt, temp3))
  return(m)
}

compute_trb_clones_aa <- function(m)
{
  m$TRB_chain_cdr3[is.na(m$TRB_chain_cdr3)] <- 0
  m <- m %>% mutate(TRB_aa_dataset_clonotype_id = lapply(TRB_chain_cdr3, temp3))
  return(m)
}

cyclone_prepare_trajectories <- function(seurat, ambiguity = 'combine', calculate_nt_clones = TRUE)
{
  mdata <- seurat@meta.data[, !colnames(seurat@meta.data) %in% 'clone_id'] %>% tibble::rownames_to_column('barcode')
  
  message("Computing clonotypes")
  
  if (calculate_nt_clones){
    message("Clonotypes are calculated using TRB nt sequences")
    clonotypes <-
      compute_trb_clones_nt(mdata) %>%
      filter(!is.na(TRB_nt_dataset_clonotype_id)) %>%
      mutate(clone_id = paste0(patient_alias, TRB_nt_dataset_clonotype_id))
  }
  else {
    message("Clonotypes are calculated using TRB aa sequences")
    clonotypes <-
      compute_trb_clones_aa(mdata) %>%
      filter(!is.na(TRB_aa_dataset_clonotype_id)) %>%
      mutate(clone_id = paste0(patient_alias, TRB_aa_dataset_clonotype_id))
  }
  message(nrow(clonotypes) - length(unique(clonotypes$barcode)), " cells with ambigiuous clonotypes identified")
  
  if (ambiguity == 'combine')
  {
    message("Combining ambiguous clonotypes")
    
    # concatenate multiple clonotypes per barcode
    clonotypes_unique <-
      clonotypes%>%
      select(barcode, clone_id) %>%
      distinct %>%
      arrange(barcode, clone_id) %>%
      group_by(barcode) %>%
      summarize(clone_id = paste0(unique(clone_id), collapse = '|')) %>%
      ungroup
    
    clonotypes <- left_join(clonotypes %>% select(!matches('^TR[A,B]|TCR')) %>% select(-clone_id) %>% distinct(), clonotypes_unique)
  } else {
    # throws all (other option could be to randomly remove duplicated entry)
    ambi_barcodes <- clonotypes %>% mutate(duplicated = duplicated(barcode)) %>% filter(duplicated) %>% pull(barcode)
    clonotypes <- clonotypes %>% filter(!barcode %in% ambi_barcodes) %>% mutate(clone_id = TRB_nt_dataset_clonotype_id)
  }
  
  seurat@meta.data <- left_join(mdata, select(clonotypes, barcode, clone_id), by = 'barcode') %>% tibble::column_to_rownames('barcode')
  
  return(seurat)
}

### Function 2 - cyclone_prepare_trajectories
# This function calculate clone trajectories and clusters based on "clone_id" and "timepoint" features in the meta data table.

cyclone_find_trajectories <- function(seurat,
                                      n_clust = 6,
                                      min_cells = 5,
                                      scale = TRUE,
                                      ambiguity = 'combine',
                                      clustering_method = 'hclust',
                                      correlation_method = 'pearson',
                                      timepoint_column = 'timepoint',
                                      downsample = FALSE,
                                      min_t_cells = 1, #1e3,
                                      drop_missing_timepoints = TRUE,
                                      classify_vdj = FALSE,
                                      delta = FALSE,
                                      reg = 0.00,
                                      min_fc = 2,
                                      min_c = 0.25,
                                      custom_order = NULL){
  stopifnot(!is.null(n_clust))
  
  if(is.null(seurat@meta.data$clone_id)) {
    seurat <- cyclone_prepare_trajectories(seurat, ambiguity = ambiguity)
  }
  # Delete previous cluster assignment from the object's meta data table
  mdata <- seurat@meta.data[, !colnames(seurat@meta.data) %in% 'k'] %>%
    tibble::rownames_to_column('barcode') %>%
    filter(!is.na(clone_id))
  
  mdata <- mdata %>% dplyr::rename('timepoint' = !!sym(timepoint_column))
  
  # Remove samples in which the number of T-cells is below the threshold
  mdata_f <- mdata %>% dplyr::count(patient_alias, timepoint) %>%
             filter(n >= min_t_cells) %>% select(-n) %>% left_join(., mdata)
  
  # Optional - Subset each sample to a set number of T-cells 
  if(downsample) {
    mdata_f <- mdata_f %>% group_by(patient_alias, timepoint) %>% 
      sample_n(min_t_cells) %>% ungroup
    message("Downsampling T-cells to ", min_t_cells, " per patient and timepoint")
  }
  
  # Compute relative abundance statistics per clone
  clonestats <- .compute_clonestats(mdata_f, classify_vdj)
  
  clonestats_f <-
    clonestats %>%
    arrange(clone_id, timepoint) %>%
    group_by(clone_id) %>%
    mutate(delta_fc = (perc + reg) / (lag(perc, order_by = timepoint, n = 1) + reg)) %>%
    mutate(delta_c = perc - lag(perc, order_by = timepoint, n = 1)) %>%
    ungroup %>%
    mutate(delta_fc = ifelse(is.na(delta_fc), 1, delta_fc)) %>% # baseline should be 1
    mutate(delta_c = ifelse(is.na(delta_c), 0, delta_c)) # baseline should be 0

  clonestats_f <- clonestats_f %>%
    group_by(clone_id) %>%
    mutate(min_perc = min(perc), max_perc = max(perc)) %>%
    ungroup

  ntimepoints <- length(unique(as.character(clonestats_f$timepoint)))

  # save fc and c to csv
  clonestats_f <- clonestats_f %>%
    group_by(clone_id) %>%
    mutate(timepoint_max = timepoint[which.max(perc)]) %>%
    mutate(timepoint_min = timepoint[which.min(perc)]) %>%
    mutate(delta_fc_max = max(delta_fc)) %>%
    mutate(delta_fc_min = min(delta_fc)) %>%
    mutate(delta_c_max = max(delta_c)) %>%
    mutate(delta_c_min = min(delta_c)) %>%
    ungroup
    
  clonestats_f %>%
    write_csv(paste0("fold change vs abs change_", ntimepoints, "_", clustering_method, ".csv"))

  clonestats_f <- clonestats_f %>%
    group_by(clone_id) %>%
    filter(any((delta_fc > min_fc & delta_c > min_c & timepoint_max != 0) | # increasing case
                  (delta_fc < 1 / min_fc & delta_c < -min_c & timepoint_max == 0))) %>% # decreasing case
    ungroup
  
  clonestats_f %>%
    write_csv(paste0("fold change vs abs change_", ntimepoints, "_", clustering_method, "_filtered.csv"))
  
  # Create a matrix of clone-level relative abundance by time point
  mat <- clonestats_f %>%
    filter(n_max >= min_cells) %>%
    select(clone_id, timepoint, perc) %>%
    tidyr::spread(timepoint, perc) %>%
    column_to_rownames('clone_id') %>%
    as.matrix() %>%
    t()
  
  ## ONLY USE IF PATIENTS WITH MISSING SAMPLES WERE REMOVED
  if (TRUE) {
    mat[is.na(mat)] <- 0
  }

  # Optional - remove missing time points
  drop_missing_timepoints <- TRUE
  if(drop_missing_timepoints) {
    mat <- mat[, colSums(is.na(mat)) == 0]
  }
  
  message("Working with ", ncol(mat), " clones")
  
  if (scale) {
    message("Scaling Clone Count Matrix")
    mat  <- (mat - rowMeans(mat, na.rm = TRUE)) / apply(mat, 1, sd, na.rm = TRUE)
  }
  
  # Compute correlation matrix between clone relative abundances
  #mat_cor <- suppressWarnings(cor(mat, method = correlation_method, use = "pairwise.complete.obs"))
  mat_cor <- cor(mat, method = correlation_method, use = "pairwise.complete.obs")

  # Compute clone clusters based on the correlation matrix
  message("Clustering Clones")
  if (clustering_method == 'hcpc')
  {
    res_pca  <- FactoMineR::PCA(mat_cor, graph=FALSE)
    res_hcpc <- FactoMineR::HCPC(res_pca, nb.clust = n_clust, min = 3, max = NULL, graph = FALSE, order = TRUE)
    
    df_k <- tibble(clone_id = rownames(res_hcpc$data.clust), k = paste0('Traj ', res_hcpc$data.clust$clust))
  }
  
  if (clustering_method == 'hclust')
  {
    #d    <- tgstat::tgs_dist(t(mat_cor))\
    d    <- dist(mat_cor, method = 'euclidean')
    hc   <- hclust(d, method = 'ward.D2')
    ct   <- cutree(hc, k = n_clust)
    if (!is.null(custom_order) & length(custom_order) == n_clust & all(custom_order %in% unique(ct))) {
      # replace cluster numbers with custom order
      ct <- factor(ct, levels = unique(ct), labels = custom_order)
    }
    df_k <- tibble(clone_id = names(ct), k = paste0('Traj ', ct))
  }
  
  if (clustering_method == 'pam')
  {
    res  <- cluster::pam(mat_cor, k = n_clust)
    df_k <- tibble(clone_id = names(res$clustering), k = paste0('Traj ', res$clustering))
  }

  if (clustering_method == 'absolute')
  {
    mat <- data.frame(t(mat))
    # move first column (Baseline) to last column
    mat <- mat[, c(2:ncol(mat), 1)]
    # assign cluster by timepoint with max frequency
    df_k <- mat %>%
      mutate(k = paste0('Traj ', apply(., 1, which.max))) %>%
      tibble::rownames_to_column('clone_id') %>%
      select(clone_id, k) %>%
      distinct(clone_id, .keep_all = TRUE)
  }
  
  # Add trajectory information back to the Seurat object
  ## Add cluster assignment column to the meta data table
  mdata_cl <- left_join(mdata, df_k) %>% 
    select(barcode, clone_id, k) %>% 
    mutate(k = ifelse(is.na(k), 'Other', k))
  seurat@meta.data <- left_join(mdata, mdata_cl) %>% column_to_rownames('barcode')
  
  ## Add clone statistics and correlation table to the object 
  seurat@misc[['mat_cor']]    <- mat_cor
  seurat@misc[['clonestats']] <- left_join(clonestats, df_k, by = 'clone_id') %>% tidyr::replace_na(list(k = 'Other'))
  
  return(seurat)
}

### Function 3 - .compute_clonestats
# A helper function for "cyclone_find_trajectories", calculates clone relative abundance and size statistics on a clone-timepoint level

.compute_clonestats <- function(seurat, classify_vdj)
{
  if(any(class(seurat) %in% 'Seurat'))
  {
    clonotypes <- seurat@meta.data
  } else {
    clonotypes <- seurat
  }
  
  clonotypes <- data.table::as.data.table(clonotypes)
  clonestats <- clonotypes[which(!is.na(clone_id)), ][, .(n = .N),
                                                      by = intersect(colnames(clonotypes),
                                                                     c('patient_alias', 'timepoint', 'treatment', 'clone_id', 'response', 'TRB_chain_cdr3'))]
  
  # Extract levels of the timepoint feature
  if(class(clonestats$timepoint) != 'factor')
  {
    levels <- levels(factor(as.character(clonotypes$timepoint)))
  } else {
    levels <- levels(clonotypes$timepoint)
  }
  # Set the data frame to a long format and fill with 0 counts where needed
  clonestats <- clonestats %>%
    tidyr::spread(timepoint, n) %>%
    tidyr::pivot_longer(!one_of('patient_alias', 'treatment', 'clone_id', 'response', 'TRB_chain_cdr3')) %>%
    group_by(patient_alias, treatment, name) %>%
    mutate(all_na = all(is.na(value)),
           value = ifelse(!all_na & is.na(value), 0, value),
           name = factor(name, levels = levels)) %>%
    ungroup %>%
    select(-all_na) %>%
    dplyr::rename(timepoint = name, n = value)
  # Compute % of total T-cells per timepoint and patient
  clonestats <- clonestats %>% group_by(patient_alias, treatment, timepoint) %>%
    mutate(n_tcells_tp = sum(n)) %>% ungroup %>%
    mutate(perc = n / n_tcells_tp * 100)
  
  # Compute per patient the minimum, maximum and range of clone sizes 
  clonestats <- clonestats %>% group_by(clone_id, patient_alias, treatment) %>% 
    mutate(n_max = max(n, na.rm = T),
           n_min = min(n, na.rm = T),
           range = n_max - n_min) %>% ungroup
  
  # Compute clone rank
  clonestats <- clonestats %>% group_by(clone_id, patient_alias, treatment) %>% 
    mutate(range_perc = max(perc, na.rm = TRUE) - min(perc, na.rm = TRUE)) %>% 
    group_by(patient_alias, treatment) %>% ungroup
  
  # Add trajectory column in case it existed in the Seurat object
  if(any(colnames(clonotypes) %in% "k")) {
    clonestats <- left_join(clonestats, distinct(select(clonotypes, clone_id, k)))
  }

  # specificity
  if (classify_vdj) {
    vdj <- read.csv("vdjdb_trb.tsv", sep="\t", header=T)
    vdj <- vdj %>% select(CDR3, Epitope.gene, Epitope.species) %>% 
                  group_by(CDR3) %>% 
                  filter(row_number() == 1) %>% # prevent duplicates
                  ungroup()
    clonestats <- merge(clonestats, vdj, by.x="TRB_chain_cdr3", by.y="CDR3", all.x=T)
  }

  if(any(class(seurat) %in% 'Seurat')) {
    seurat@misc$clonestats <- clonestats
    return(seurat)
  } else {
    return(clonestats)
  }
}

## Plotting functions
### Plotting function 1 - Plot correlation heatmap of clone dynamics 
plot_trajectories_cor <- function(seurat) {
  mat_cor <- seurat@misc$mat_cor
  # number of each k
  clone_order <- seurat@meta.data %>% 
    filter(k != 'Other', !is.na(k)) %>%
    select(clone_id, k) %>%
    distinct %>% 
    arrange(k)
  
  # number of each k
  k_counts <- clone_order %>%
    group_by(k) %>%
    summarise(n = n())
  
  annotation <- data.frame(Trajectory = clone_order$k)
  rownames(annotation) <- clone_order$clone_id

  colz <- auto_colset(length(unique(clone_order$k)), colz="selected")
  names(colz) <- unique(clone_order$k)
  ann_colors = list(
    Trajectory = colz
  )

  ha <- ComplexHeatmap::HeatmapAnnotation(Trajectory = clone_order$k, 
                                          col = ann_colors, which="row", 
                                          show_annotation_name=c(Trajectory=FALSE))

  ht <- ComplexHeatmap::Heatmap(mat_cor[clone_order$clone_id, clone_order$clone_id],
                                col = rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")),
                                left_annotation=ha, name="Correlation",
                                cluster_columns=FALSE, cluster_rows=FALSE,
                                show_row_names = FALSE, show_column_names = FALSE)
  return (ht)
}

### Plotting function 2- Plot clone trajectories by cluster across time points
plot_trajectories_area <- function(seurat,
                                   facet_by = 'k',
                                   group_by = 'clone_id',
                                   col_by = 'k',
                                   ncols = 4,
                                   summary_function = sum,
                                   ymax = NULL,
                                   split_treatment = FALSE,
                                   normalize = FALSE,
                                   exclude = FALSE) {
  ggdata <- seurat@misc$clonestats %>% 
    filter(k != 'Other') %>%
    group_by(k) %>%
    mutate(rank_cl = base::rank(n_max, ties.method = 'first')) %>%
    ungroup() %>%
    arrange(k, rank_cl) %>%
    mutate(clone_id = forcats::fct_inorder(factor(clone_id))) %>% 
    filter(!is.na(timepoint))
  
  # if exclude trajectory that goes down
  if (exclude) {
    last_traj <- ggdata %>% pull(k) %>% unique %>% sort %>% tail(1)
    ggdata <- ggdata %>% filter(k != last_traj)
  }
  
  if (col_by != "k") {
    # order clone_id by treatment, then by clone size
    ggdata <- ggdata %>%
      group_by(!!sym(col_by), clone_id) %>%
      mutate(rank_cl = base::rank(n_max, ties.method = 'first')) %>%
      ungroup() %>%
      arrange(!!sym(col_by), rank_cl) %>%
      mutate(clone_id = forcats::fct_inorder(factor(clone_id))) %>% 
      filter(!is.na(timepoint))
  }
  
  if(group_by != 'clone_id') {
    ggdata <- ggdata %>%
      group_by(timepoint, patient_alias, treatment, k) %>%
      summarize(perc = summary_function(perc)) %>%
      ungroup
  }

  # get number of unique patients per treatment
  if (normalize) {
    ggdata <- ggdata %>%
        group_by(treatment) %>%
        mutate(npatients = n_distinct(patient_alias)) %>%
        ungroup() %>%
        mutate(perc = perc / npatients * min(npatients))
  }
  
  ggdata <- ggdata %>% filter(!is.na(perc))
  if(group_by == 'patient_alias') {
    p <- ggplot(ggdata, aes(timepoint, perc,
                            group = !!sym(group_by),
                            col = factor(!!sym(group_by)))) +
      geom_line() +
      theme_bw() + 
      geom_point() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      facet_wrap(vars(!!sym(facet_by)), scales = 'free', ncol = ncols) +
      {if(is.numeric(ymax)) coord_cartesian(ylim = (c(0, ymax)))}
  } else {
    colz <- auto_colset(length(unique(ggdata[[col_by]])), colz=ifelse(col_by == "k", "selected", "paired"))
    p <- ggplot(ggdata, aes(timepoint, perc,
                            group = !!sym(group_by),
                            fill = factor(!!sym(col_by)))) +
      geom_area(color = 'black', lwd = 0.05) +
      theme_bw() + 
      labs(fill=col_by) +
      scale_fill_manual(values=colz) +
      {if (split_treatment) {
        facet_grid(vars(treatment), vars(!!sym(facet_by)), scales = 'free', labeller=labeller(treatment = fix_pd1))#, space = 'free') #if want box to be proportional
      } else {
        facet_wrap(vars(!!sym(facet_by)), scales = 'fixed', ncol = ncols)
      }} +
      {if(is.numeric(ymax)) coord_cartesian(ylim = (c(0, ymax)))}
  }
  if(col_by == 'clone_id'){
    p <- p + theme(legend.position = 'none')
  }
  
  p <- p + labs(x = 'Weeks', y = 'Clone Abundance (%)')
  return(p)
}

### Plotting function 3 - Plot clone trajectories by cluster and cell type across time points
plot_trajectories_ctype <- function(seurat,
                                    normalization = 'none',
                                    delta = 'delta',
                                    trajectories = NULL,
                                    ctype_exclude = NULL,
                                    top_n = Inf,
                                    regularization = 5,
                                    ncols = 2,
                                    facet_by = "k",
                                    group_by = 'celltype',
                                    norm_cell_n = FALSE,
                                    show_labels = TRUE,
                                    label_top_n = 3,
                                    exclude = FALSE)
{
  ggdata <- seurat@meta.data %>% filter(!is.na(clone_id))
  
  # Set timepoint column as a factor
  if(class(ggdata$timepoint) != 'factor') {
    ggdata$timepoint <- as.factor(ggdata$timepoint)
  }
  
  # Normalize time-points to same # t-cells
  if(norm_cell_n) {
    set.seed(19)
    n_cells <- ggdata %>% count(timepoint) %>% pull(n) %>% min
    ggdata  <- ggdata %>% group_by(timepoint) %>% sample_n(n_cells, replace = FALSE)
    message('Normalizing to ', n_cells, ' number of cells per time point')
  }

  # exclude trajectory that goes down (if used absolute clustering)
  if (exclude) {
    last_traj <- ggdata %>% filter(k != "Other") %>% pull(k) %>% unique %>% sort %>% tail(1)
    ggdata <- ggdata %>% filter(k != last_traj)
  }

  ggdata2 <- ggdata %>% filter(k != 'Other', !is.na(timepoint)) %>%
    group_by(timepoint, k, celltype, Phase) %>%
    summarise(phase_counts=n()) %>%
    ungroup() %>%
    group_by(timepoint, k, celltype) %>%
    summarise(cycling = sum(phase_counts[Phase == 'G2M']) / sum(phase_counts) * 100) %>%
    ungroup()

  ggdata <- ggdata %>% filter(!is.na(timepoint)) %>%
    group_by(timepoint, patient_alias) %>%
    mutate(count=n()) %>%
    ungroup() %>%
    filter(k != 'Other') %>%
    group_by(timepoint, k, clone_id, celltype) %>%
    summarise(freq = n() / count) %>%
    ungroup() %>%
    distinct() %>%
    tidyr::complete(timepoint, k, celltype, fill = list(freq = 0)) %>%
    group_by(timepoint, k, celltype) %>%
    summarise(n = sum(freq) * 100)

  ggdata <- ggdata %>% left_join(ggdata2, by = c('timepoint', 'k', 'celltype'))

  y_label_suffix <- '%'
  
  # Format y axis values either as delta (default) or FC
  if(delta == 'fc') {
    ggdata <- ggdata %>%
      group_by(k, celltype) %>%
      mutate(delta = log2((n + regularization) / (n[timepoint == levels(timepoint)[1]] + regularization))) %>%
      ungroup
    
    y_label <- 'log2FC'
  }
  
  if(delta == 'delta') {
    ggdata <- ggdata %>%
      group_by(k, celltype) %>%
      mutate(delta = n - n[timepoint == levels(timepoint)[1]]) %>%
      ungroup
    
    y_label <- paste0('Delta (', y_label_suffix, ')')
  }

  # Optional - remove specific cell types
  if(!is.null(ctype_exclude)) {
    ggdata <- ggdata %>% filter(!celltype %in% ctype_exclude)
  }
  
  ggdata <- ggdata %>% 
    mutate(celltype = forcats::fct_lump_n(celltype,
                                           n = top_n,
                                           w = abs(delta),
                                           other_level = 'other')) %>%
    filter(celltype != 'other')
  # Optional - remove specific trajectories
  if(!is.null(trajectories)) {
    ggdata <- ggdata %>% filter(k %in% trajectories)
  }
  # Plotting
  p <- ggdata %>%
    ggplot(aes(timepoint, delta, group = !!sym(group_by), col = !!sym(group_by))) +
    geom_point() +
    geom_point(aes(size = cycling), alpha=0.3) +
    geom_line() +
    scale_size_binned_area(n.breaks = 8) +
    facet_wrap(vars(!!sym(facet_by)), scales = 'free', ncol = ncols) +
    theme_bw() + 
    labs(x = "Weeks", y = y_label)
  
  if(show_labels) {
    ggdata_text_repel <- ggdata %>%
      group_by(celltype, k) %>%
      mutate(keep = abs(delta) == max(abs(delta))) %>%
      filter(keep) %>% 
      group_by(celltype, k) %>%
      sample_n(1)

    # choose labels that exceed threshold
    ggdata_text_repel <- ggdata_text_repel %>% 
      group_by(k) %>%
      filter(delta > fc_filter_plotting) %>%
      ungroup()
    p <- p + ggrepel::geom_text_repel(data = ggdata_text_repel,
                                      aes(label = ""),
                                      nudge_x = 0.35,
                                      size = 4,
                                      show_guide=F)
  }
  
  return(p)
}

### Plotting function 3 - Plot clone trajectories by cluster and cell type across time points
plot_trajectories_ctype_treatment <- function(seurat,
                                    normalization = 'none',
                                    delta = 'delta',
                                    trajectories = NULL,
                                    ctype_exclude = NULL,
                                    top_n = Inf,
                                    regularization = 5,
                                    ncols = 2,
                                    facet_by = "k",
                                    group_by = 'celltype',
                                    norm_cell_n = FALSE,
                                    show_labels = TRUE,
                                    exclude = FALSE,
                                    selected = NULL,
                                    label_top_n = 3,
                                    filter_magnitude = FALSE)
{
  ggdata <- seurat@meta.data %>% filter(!is.na(clone_id))
  
  # Set timepoint column as a factor
  if(class(ggdata$timepoint) != 'factor') {
    ggdata$timepoint <- as.factor(ggdata$timepoint)
  }
  
  # Normalize time-points to same # t-cells
  if(norm_cell_n) {
    set.seed(19)
    n_cells <- ggdata %>% count(timepoint, treatment) %>% pull(n) %>% min
    ggdata  <- ggdata %>% group_by(timepoint, treatment) %>% sample_n(n_cells, replace = FALSE)
    message('Normalizing to ', n_cells, ' number of cells per time point')
  }

  if (exclude) {
    last_traj <- ggdata %>% filter(k != "Other") %>% pull(k) %>% unique %>% sort %>% tail(1)
    ggdata <- ggdata %>% filter(k != last_traj)
  }

  ggdata2 <- ggdata %>% filter(k != 'Other', !is.na(timepoint)) %>%
    group_by(timepoint, treatment, k, celltype, Phase, patient_alias) %>%
    summarise(phase_counts=n()) %>%
    ungroup() %>%
    group_by(timepoint, treatment, k, celltype, patient_alias) %>%
    summarise(cycling = sum(phase_counts[Phase == 'G2M']) / sum(phase_counts)) %>%
    ungroup() %>%
    group_by(timepoint, treatment, k, celltype) %>%
    summarise(cycling = mean(cycling)) # average patient

  ggdata <- ggdata %>% filter(!is.na(timepoint)) %>%
    group_by(timepoint, treatment, patient_alias) %>%
    mutate(count=n()) %>%
    ungroup() %>%
    filter(k != 'Other') %>%
    group_by(timepoint, treatment, k, clone_id, celltype, patient_alias) %>%
    summarise(freq = n() / count) %>%
    ungroup() %>%
    distinct() %>%
    tidyr::complete(timepoint, treatment, k, celltype, patient_alias, fill = list(freq = 0)) %>%
    group_by(timepoint, treatment, k, celltype, patient_alias) %>%
    summarise(n = sum(freq) * 100) %>%
    ungroup() %>%
    group_by(timepoint, treatment, k, celltype) %>%
    summarise(n = mean(n)) # average patient

  ggdata <- ggdata %>% left_join(ggdata2, by = c('timepoint', 'treatment', 'k', 'celltype'))

  y_label_suffix <- '%'
  
  # Format y axis values either as delta (default) or FC
  if(delta == 'fc') {
    ggdata <- ggdata %>%
      group_by(k, treatment, celltype) %>%
      mutate(delta = log2((n + regularization) / (n[timepoint == levels(timepoint)[1]] + regularization))) %>%
      ungroup
    
    y_label <- 'log2FC'
  } else if(delta == 'delta') {
    ggdata <- ggdata %>%
      group_by(k, treatment, celltype) %>%
      mutate(delta = n - n[timepoint == levels(timepoint)[1]]) %>%
      ungroup
    
    y_label <- paste0('Delta (', y_label_suffix, ')')
  } else {
    ggdata <- ggdata %>% mutate(delta = n)
    y_label <- "Absolute %"
  }



  # Optional - remove specific cell types
  if(!is.null(ctype_exclude)) {
    ggdata <- ggdata %>% filter(!celltype %in% ctype_exclude)
  }
  
  ggdata <- ggdata %>% 
    mutate(celltype = forcats::fct_lump_n(celltype,
                                           n = top_n,
                                           w = abs(delta),
                                           other_level = 'other')) %>%
    filter(celltype != 'other')
  # Optional - specific trajectories
  if(!is.null(trajectories)) {
    ggdata <- ggdata %>% filter(k %in% trajectories)
  }
  # Plotting
  p <- ggdata %>% 
    ggplot(aes(timepoint, delta, group = !!sym(group_by), col = !!sym(group_by))) +
    geom_point(size=1.5) +
    geom_point(aes(size = cycling), alpha=0.3) +
    geom_line(size=0.5) +
    facet_grid(vars(treatment), vars(!!sym(facet_by)), scales = 'free', space = 'fixed', labeller=labeller(treatment = fix_pd1)) +
    theme_bw() + 
    labs(x = "Weeks", y = y_label) +
    scale_radius(range=c(0.5, 8), breaks=c(0, 0.25, 0.50, 0.75, 1.00), labels=c(0, 25, 50, 75, 100)) +
    guides(size=guide_legend(title="% cycling"))
  
  if(show_labels) {
    ggdata_text_repel <- ggdata %>%
      group_by(celltype, k) %>%
      mutate(keep = abs(delta) == max(abs(delta))) %>%
      filter(keep) %>% 
      group_by(celltype, k) %>%
      sample_n(1)
    # choose labels that exceed threshold
    ggdata_text_repel <- ggdata_text_repel %>% 
      group_by(k) %>%
      filter(delta > fc_filter_plotting) %>%
      ungroup()
    
    p <- p + ggrepel::geom_text_repel(data = ggdata_text_repel,
                                      aes(label = !!sym(group_by)),
                                      nudge_x = 0.35,
                                      size = 4,
                                      show.legend=F)
  }
  
  return(p)
}

### Plotting function 3 - Plot clone trajectories by cluster and cell type across time points
plot_trajectories_ctype_treatment_patient <- function(seurat,
                                    normalization = 'none',
                                    delta = 'delta',
                                    trajectories = NULL,
                                    ctypes = NULL,
                                    top_n = Inf,
                                    regularization = 5,
                                    ncols = 2,
                                    facet_by = "k",
                                    group_by = 'celltype',
                                    norm_cell_n = FALSE,
                                    show_labels = TRUE,
                                    exclude = FALSE,
                                    selected = NULL,
                                    treatments = NULL,
                                    label_top_n = 3,
                                    filter_magnitude = FALSE)
{
  ggdata <- seurat@meta.data %>% filter(!is.na(clone_id))
  
  # Set timepoint column as a factor
  if(class(ggdata$timepoint) != 'factor') {
    ggdata$timepoint <- as.factor(ggdata$timepoint)
  }
  
  # Normalize time-points to same # t-cells
  if(norm_cell_n) {
    set.seed(19)
    n_cells <- ggdata %>% count(timepoint, treatment) %>% pull(n) %>% min
    ggdata  <- ggdata %>% group_by(timepoint, treatment) %>% sample_n(n_cells, replace = FALSE)
    message('Normalizing to ', n_cells, ' number of cells per time point')
  }

  if (exclude) {
    last_traj <- ggdata %>% filter(k != "Other") %>% pull(k) %>% unique %>% sort %>% tail(1)
    ggdata <- ggdata %>% filter(k != last_traj)
  }

  ggdata <- ggdata %>% filter(!is.na(timepoint)) %>%
    group_by(timepoint, treatment, patient_alias) %>%
    mutate(count=n()) %>%
    ungroup() %>%
    filter(k != 'Other') %>%
    group_by(timepoint, treatment, k, clone_id, celltype, patient_alias) %>%
    summarise(freq = n() / count) %>%
    ungroup() %>%
    distinct()
  
  if (!is.null(selected)) {
    ggdata <- ggdata %>% filter(k %in% selected)
  }

  if (!is.null(treatments)) {
    ggdata <- ggdata %>% filter(treatment %in% treatments)
  }

  ggdata <- ggdata %>%
    tidyr::complete(timepoint, treatment, k, celltype, patient_alias, fill = list(freq = 0)) %>%
    group_by(timepoint, treatment, k, celltype, patient_alias) %>%
    summarise(n = sum(freq) * 100)

  y_label_suffix <- '%'
  
  # Format y axis values either as delta (default) or FC
  if(delta == 'fc') {
    ggdata <- ggdata %>%
      group_by(k, treatment, celltype) %>%
      mutate(delta = log2((n + regularization) / (n[timepoint == levels(timepoint)[1]] + regularization))) %>%
      ungroup
    
    y_label <- 'log2FC'
  } else if(delta == 'delta') {
    ggdata <- ggdata %>%
      group_by(k, treatment, celltype) %>%
      mutate(delta = n - n[timepoint == levels(timepoint)[1]]) %>%
      ungroup
    
    y_label <- paste0('Delta (', y_label_suffix, ')')
  } else {
    ggdata <- ggdata %>% mutate(delta = n)
    y_label <- "Absolute %"
  }

  # Optional - specific cell types
  if(!is.null(ctypes)) {
    ggdata <- ggdata %>% filter(celltype %in% ctypes)
  }
  
  ggdata <- ggdata %>% 
    mutate(celltype = forcats::fct_lump_n(celltype,
                                           n = top_n,
                                           w = abs(delta),
                                           other_level = 'other')) %>%
    filter(celltype != 'other')

  # Plotting
  p <- ggdata %>% 
    ggplot(aes(timepoint, delta, group = patient_alias, col = celltype)) +
    geom_point(size=1.5) +
    geom_line(size=0.5) +
    facet_grid(vars(celltype), vars(k), scales = 'fixed', space = 'free', labeller=labeller(treatment = fix_pd1)) +
    theme_bw() + 
    labs(x = "Weeks", y = y_label)
  
  return(p)
}

### Plotting function 3 - Plot clone trajectories by cluster and cell type across time points
plot_trajectories_ctype_comp <- function(seurat,
                                    normalization = 'none',
                                    delta = 'delta',
                                    trajectories = NULL,
                                    ctype_exclude = NULL,
                                    top_n = Inf,
                                    regularization = 5,
                                    ncols = 2,
                                    facet_by = "k",
                                    group_by = 'celltype',
                                    norm_cell_n = FALSE,
                                    show_labels = TRUE,
                                    label_top_n = 3,
                                    exclude = FALSE)
{
  ggdata <- seurat@meta.data %>% filter(!is.na(clone_id))
  
  # Set timepoint column as a factor
  if(class(ggdata$timepoint) != 'factor') {
    ggdata$timepoint <- as.factor(ggdata$timepoint)
  }
  
  # Normalize time-points to same # t-cells
  if(norm_cell_n) {
    set.seed(19)
    n_cells <- ggdata %>% count(timepoint) %>% pull(n) %>% min
    ggdata  <- ggdata %>% group_by(timepoint) %>% sample_n(n_cells, replace = FALSE)
    message('Normalizing to ', n_cells, ' number of cells per time point')
  }

  # ggdata <- ggdata %>% filter(k != 'Other', !is.na(timepoint))  %>%
  #   count(timepoint, k, celltype) %>%
  #   ungroup() %>%
  #   tidyr::complete(timepoint, k, celltype, fill = list(n = 0))
  # if exclude trajectory that goes down
  if (exclude) {
    last_traj <- ggdata %>% filter(k != "Other") %>% pull(k) %>% unique %>% sort %>% tail(1)
    ggdata <- ggdata %>% filter(k != last_traj)
  }
  ggdata <- ggdata %>% filter(k != 'Other', !is.na(timepoint)) %>%
    group_by(timepoint, k, celltype) %>%
    summarise(n=n(), cycling_count=sum(cycling1 > cycling_cutoff), cycling = cycling_count / n) %>%
    ungroup() %>%
    tidyr::complete(timepoint, k, celltype, fill = list(n = 0))

  ggdata$cycling <- scales::rescale(ggdata$cycling, to=c(1, 8), from=range(0, 0.1))
  
  # Format abundance values either as absolute counts (default) or relative abundance
  if(normalization != 'absolute')
  {
    ggdata <- ggdata %>% 
      group_by(timepoint, k) %>%
      mutate(n = n / sum(n) * 100) %>%
      ungroup %>% 
      filter(!is.na(n)) # if no cells at timepoint + trajectory will return NaN
    
    y_label_suffix <- '%'
  } else {
    y_label_suffix <- '# cells'
  }
  
  # Format y axis values either as delta (default) or FC
  if(delta == 'fc') {
    ggdata <- ggdata %>%
      group_by(k, celltype) %>%
      mutate(delta = log2((n + regularization) / (n[timepoint == levels(timepoint)[1]] + regularization))) %>%
      ungroup
    
    y_label <- 'log2FC'
  }
  
  if(delta == 'delta') {
    ggdata <- ggdata %>%
      group_by(k, celltype) %>%
      mutate(delta = n - n[timepoint == levels(timepoint)[1]]) %>%
      ungroup
    
    y_label <- paste0('Delta (', y_label_suffix, ')')
  }

  # Optional - remove specific cell types
  if(!is.null(ctype_exclude)) {
    ggdata <- ggdata %>% filter(!celltype %in% ctype_exclude)
  }
  
  ggdata <- ggdata %>% 
    mutate(celltype = forcats::fct_lump_n(celltype,
                                           n = top_n,
                                           w = abs(delta),
                                           other_level = 'other')) %>%
    filter(celltype != 'other')
  # Optional - remove specific trajectories
  if(!is.null(trajectories)) {
    ggdata <- ggdata %>% filter(k %in% trajectories)
  }
  # Plotting
  p <- ggdata %>% mutate(timepoint = factor(as.integer(timepoint))) %>%
    ggplot(aes(timepoint, delta, group = !!sym(group_by), col = !!sym(group_by))) +
    geom_point() +
    #geom_point(aes(size = cycling), alpha=0.3) +
    geom_line() +
    scale_size_binned_area(n.breaks = 8) +
    facet_wrap(vars(!!sym(facet_by)), scales = 'free', ncol = ncols) +
    theme_bw() + 
    labs(y = y_label)
  
  if(show_labels) {
    ggdata_text_repel <- ggdata %>% mutate(timepoint = factor(as.integer(timepoint))) %>%
      group_by(celltype, k) %>%
      mutate(keep = abs(delta) == max(abs(delta))) %>%
      filter(keep) %>% 
      group_by(celltype, k) %>%
      sample_n(1)

    ggdata_text_repel <- ggdata_text_repel %>% 
      group_by(k) %>%
      filter(abs(delta) > fc_filter_plotting) %>%
      ungroup()
    
    p <- p + ggrepel::geom_text_repel(data = ggdata_text_repel,
                                      aes(label = !!sym(group_by)),
                                      nudge_x = 0.35,
                                      size = 4)
  }
  
  return(p)
}

### Plotting function 3 - Plot clone trajectories by cluster and cell type across time points
plot_trajectories_ctype_treatment_comp <- function(seurat,
                                    normalization = 'none',
                                    delta = 'delta',
                                    trajectories = NULL,
                                    ctype_exclude = NULL,
                                    top_n = Inf,
                                    regularization = 5,
                                    ncols = 2,
                                    facet_by = "k",
                                    group_by = 'celltype',
                                    norm_cell_n = FALSE,
                                    show_labels = TRUE,
                                    label_top_n = 3,
                                    exclude = FALSE)
{
  ggdata <- seurat@meta.data %>% filter(!is.na(clone_id))
  
  # Set timepoint column as a factor
  if(class(ggdata$timepoint) != 'factor') {
    ggdata$timepoint <- as.factor(ggdata$timepoint)
  }
  
  # Normalize time-points to same # t-cells
  if(norm_cell_n) {
    set.seed(19)
    n_cells <- ggdata %>% count(timepoint, treatment) %>% pull(n) %>% min
    ggdata  <- ggdata %>% group_by(timepoint, treatment) %>% sample_n(n_cells, replace = FALSE)
    message('Normalizing to ', n_cells, ' number of cells per time point')
  }

  if (exclude) {
    last_traj <- ggdata %>% filter(k != "Other") %>% pull(k) %>% unique %>% sort %>% tail(1)
    ggdata <- ggdata %>% filter(k != last_traj)
  }

  ggdata <- ggdata %>% filter(k != 'Other', !is.na(timepoint)) %>%
    group_by(timepoint, treatment, k, celltype) %>%
    summarise(n=n(), cycling_count=sum(cycling1 > cycling_cutoff), cycling = cycling_count / n) %>%
    ungroup() %>%
    tidyr::complete(timepoint, treatment, k, celltype, fill = list(n = 0))

  ggdata$cycling <- scales::rescale(ggdata$cycling, to=c(1, 8), from=range(0, 0.1))
  
  # Format abundance values either as absolute counts (default) or relative abundance
  if(normalization != 'absolute')
  {
    ggdata <- ggdata %>% 
      group_by(timepoint, treatment, k) %>%
      mutate(n = n / sum(n) * 100) %>%
      ungroup %>% 
      filter(!is.na(n)) # if no cells at timepoint + trajectory will return NaN
    
    y_label_suffix <- '%'
  } else {
    y_label_suffix <- '# cells'
  }
  
  # Format y axis values either as delta (default) or FC
  if(delta == 'fc') {
    ggdata <- ggdata %>%
      group_by(k, treatment, celltype) %>%
      mutate(delta = log2((n + regularization) / (n[timepoint == levels(timepoint)[1]] + regularization))) %>%
      ungroup
    
    y_label <- 'log2FC'
  }
  
  if(delta == 'delta') {
    ggdata <- ggdata %>%
      group_by(k, treatment, celltype) %>%
      mutate(delta = n - n[timepoint == levels(timepoint)[1]]) %>%
      ungroup
    
    y_label <- paste0('Delta (', y_label_suffix, ')')
  }

  # Optional - remove specific cell types
  if(!is.null(ctype_exclude)) {
    ggdata <- ggdata %>% filter(!celltype %in% ctype_exclude)
  }
  
  ggdata <- ggdata %>% 
    mutate(celltype = forcats::fct_lump_n(celltype,
                                           n = top_n,
                                           w = abs(delta),
                                           other_level = 'other')) %>%
    filter(celltype != 'other')
  # Optional - remove specific trajectories
  if(!is.null(trajectories)) {
    ggdata <- ggdata %>% filter(k %in% trajectories)
  }
  # Plotting
  p <- ggdata %>% mutate(timepoint = factor(as.integer(timepoint))) %>%
    ggplot(aes(timepoint, delta, group = !!sym(group_by), col = !!sym(group_by))) +
    geom_point() +
    #geom_point(aes(size = cycling), alpha=0.3) +
    geom_line() +
    scale_size_binned_area(n.breaks = 8) +
    #facet_wrap(vars(!!sym(facet_by)), scales = 'free', ncol = ncols) +
    facet_grid(vars(treatment), vars(!!sym(facet_by)), scales = 'free', space = 'free', labeller=labeller(treatment = fix_pd1)) +
    theme_bw() + 
    labs(y = y_label)
  
  if(show_labels) {
    ggdata_text_repel <- ggdata %>% mutate(timepoint = factor(as.integer(timepoint))) %>%
      group_by(celltype, treatment, k) %>%
      mutate(keep = abs(delta) == max(abs(delta))) %>%
      filter(keep) %>% 
      group_by(celltype, treatment, k) %>%
      sample_n(1)
    
    ggdata_text_repel <- ggdata_text_repel %>% 
      group_by(k) %>%
      filter(abs(delta) > fc_filter_plotting) %>%
      ungroup()
    
    p <- p + ggrepel::geom_text_repel(data = ggdata_text_repel,
                                      aes(label = !!sym(group_by)),
                                      nudge_x = 0.35,
                                      size = 4)
  }
  
  return(p)
}

