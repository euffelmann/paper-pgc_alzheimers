miami_f <- function(
    sst_1,
    sst_2,
    label_1,
    label_2,
    loci_1 = NULL,
    loci_2 = NULL,
    annotate_loci = F,
    max_log_p_value = 20
    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    ### input:
    ##    sst_1 & sst_2:      sumstats in GWAS catalog format
    ##    label_1 & label_2:  labels for plot facets
    ##    loci_1 & loci_2:    locus index variants + genes
    ##    max_log_p_value:    y-axis limit    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    
    ) { 
  
  sst <- NULL
  ## add gene labels
  if (annotate_loci) {
    
    sst[[1]] <- sst_1 %>%
      left_join(loci_1, by = c("variant_id" = "index_variant_id"))
    
    sst[[2]] <- sst_2 %>%
      left_join(loci_2, by = c("variant_id" = "index_variant_id"))
    
  } else {
    
    sst[[1]] <- sst_1
    sst[[2]] <- sst_2
    
  }

  ## remove 99% of SNPs with p_value < 0.01 to speed up plotting
  for (i in c(1,2)) {
    sst[[i]] <- sst[[i]] %>%
      mutate(
        p_value = as.numeric(p_value),
        keep = ifelse(p_value > 0.01, runif(n()) < 0.01, TRUE)
      ) %>%
      filter(keep)
  }
  
  ## add labels
  label_1 <- paste0(label_1,"~(N[eff]==", round(max(sst[[1]]$neff)), ")")
  label_2 <- paste0(label_2,"~(N[eff]==", round(max(sst[[2]]$neff)), ")")
  sst[[1]] <- sst[[1]] %>%
    mutate(label = label_1)
  sst[[2]] <- sst[[2]] %>%
    mutate(label = label_2)
  
  sst_combined <- rbind(sst[[1]], sst[[2]]) %>%
    mutate(
      log10_p = case_when(
        label == label_1 ~ -log10(p_value),
        label == label_2 ~  log10(p_value)
      )) %>%
    filter(abs(log10_p) < max_log_p_value) %>%
    arrange(chromosome, base_pair_location)
  
  
  temp <- sst_combined %>%
    # Compute chromosome size
    group_by(chromosome) %>% 
    summarise(chr_len=max(base_pair_location)) %>% 
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
    dplyr::select(-chr_len) %>%
    # Add this info to the initial dataset
    left_join(sst_combined, ., by=c("chromosome"="chromosome")) %>%
    # Add a cumulative position of each gene and mark genes that will be highlighted
    arrange(chromosome, base_pair_location) %>%
    mutate(
      BPcum = base_pair_location + tot
    )
  
  ## the dummy object is needed to fix the scale limits to the same range for both manhattan plots
  dummy <- data.frame(BPcum = rep(range(temp$BPcum), 2),
                      log10_p = c(0, round(max_log_p_value, digits = -1) + 2, 0, round(-max_log_p_value, digits = -1) - 2), # manually change to desired scale-limits
                      label = c(label_1, label_1, label_2, label_2), stringsAsFactors = F)
  
  # prepare X axis
  axisdf <- temp %>% group_by(chromosome) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  axis_labels_df <- axisdf %>% 
    mutate(log10_p = 0, label = label_2)  # assign to only one facet
  
  ## make the plot
  p <- temp %>%
    ggplot(aes(x=BPcum, y=log10_p)) +
    # Show all points
    geom_point(aes(color=as.factor(chromosome)), alpha=1, size=1.3, show.legend = FALSE) +
    scale_color_manual(values = rep(c("#3690c0", "#bdbdbd"), 23 )) +
    # custom X axis:
    scale_x_continuous(label = axisdf$chromosome, breaks= axisdf$center ) +
    facet_wrap(~label, scales = "free_y", ncol = 1, strip.position = "right", labeller = label_parsed) +
    scale_y_continuous(
      expression(-log[10](italic(P))),
      expand = c(0.015, 0),
      breaks = seq(
        round(-max_log_p_value, digits = -1),
        round(max_log_p_value, digits = -1),
        5
      ),
      labels = abs(seq(
        round(-max_log_p_value, digits = -1),
        round(max_log_p_value, digits = -1),
        5
      ))
    ) +
    geom_hline(
      data = . %>% group_by(label) %>% slice(1),
      aes(yintercept = ifelse(log10_p > 0,-log10(5e-8), log10(5e-8))),
      color = "#D82148",
      linetype = "dashed"
    ) +
    ## this is needed to fix the scale limits to the same range for both manhattan plots
    geom_blank(data = dummy) +
    geom_text(
      data = axis_labels_df,
      aes(x = center, y = log10_p, label = chromosome),
      inherit.aes = FALSE,
      vjust = -0.5,
      size = 3
    ) +
    # Custom the theme:
    coord_cartesian(clip = 'off') +
    theme_classic() +
    theme( 
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      strip.background = element_blank(),
      legend.position="top",
      legend.text = element_text(size=8),
      legend.title = element_text(size=8),
      legend.key.size = unit(0.3,"line")
    ) +
    xlab("Chromsome")

  ## add gene labels to index variants
  if (annotate_loci) {
    
    p <- p +
      shadowtext::geom_shadowtext(
        aes(label = ifelse(p_value < 5e-8 & !is.na(symbol) & log10_p < 0, symbol, "")),
        angle = 90,
        hjust = 1.3,
        size = 3,
        colour = "#bd0026",
        bg.colour = "white",
        bg.r = 0.2
      ) +
      shadowtext::geom_shadowtext(
        aes(label = ifelse(p_value < 5e-8 & !is.na(symbol) & log10_p > 0, symbol, "")),
        angle = 90,
        hjust = -0.2,
        size = 3,
        colour = "#bd0026",
        bg.colour = "white",
        bg.r = 0.2
      )
  }
  
  
  return(p)
    
}

  