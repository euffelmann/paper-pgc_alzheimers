manhattan_f <- function(
    sst, 
    annotate_loci = F, 
    loci = NULL, 
    loci_p_value = 5e-8, 
    filter_snps = T, 
    shadow = F) {
  
  ### input:
  ##     sst:               sumstats in GWAS catalog format
  ##     annotate_loci:     annotate loci with gene names?
  ##     loci:              dataframe with gene names for index SNPs
  ##     loci_p_value:      label SNPs with p_values below this
  ##     filter_snps:       filter majority of SNPs without signal to speed up
  ##                        plotting?
  ##     max_log_p_value:   cap y axis
  ##     shadow:            add white background to labels
  
  ## remove 99% of SNPs with p_value < 0.01 to speed up plotting
  if (filter_snps) {
    sst <- sst %>%
      mutate(
        p_value = as.numeric(p_value),
        keep = ifelse(p_value > 0.001, runif(n()) < 0.15, TRUE)
      ) %>%
      filter(keep)
  }
  
  max_log_p_value <- max(-log10(sst$p_value[!(sst$chromosome == 19 & sst$base_pair_location > 4e7 & sst$base_pair_location < 5e7)])) + 15
  
  ## add gene labels to index variants + filter SNPs + add APOE label
  if (annotate_loci & any(sst$p_value < 5e-8)) {
    sst <- sst %>%
      left_join(loci, by = "variant_id") %>%
      filter(-log10(p_value) < max_log_p_value) %>%
      group_by(region = (chromosome == 19 & base_pair_location > 4e7 & base_pair_location < 5e7)) %>%
      mutate(gene = if_else(region & p_value == min(p_value), "APOE", gene)) %>%
      ungroup() %>%
      dplyr::select(-region)
  }

  don <- sst %>% 
    # compute chromosome size
    group_by(chromosome) %>% 
    summarise(chr_len=max(base_pair_location)) %>% 
    # calculate cumulative position of each chromosome
    mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
    dplyr::select(-chr_len) %>%
    # add this info to the initial dataset
    left_join(sst, ., by=c("chromosome"="chromosome")) %>%
    # add a cumulative position of each SNP
    arrange(chromosome, base_pair_location) %>%
    mutate( BPcum=base_pair_location+tot)
  
  axisdf = don %>%
    group_by(chromosome) %>%
    summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  p <- don %>%
    mutate(p_value = as.numeric(p_value)) %>%
    ggplot(aes(x=BPcum, y=-log10(p_value))) +
    geom_point( aes(color=as.factor(chromosome)), alpha=0.8, size=0.8) +
    scale_color_manual(values = rep(c("#3690c0", "#fc8d59"), 22)) +
    geom_hline(yintercept = -log10(5e-8), colour = "black", linetype = "dashed") +
    # custom X axis:
    scale_x_continuous(label = axisdf$chromosome, breaks= axisdf$center ) +
    scale_y_continuous(limit = c(0, max_log_p_value + 5), expand = c(0.01, 0) ) + # remove space between plot area and x axis
    theme_classic() +
    theme( 
      legend.position="none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(size=8),
      plot.margin = mar gin(30, 10, 10, 10) ) +
    coord_cartesian(clip = 'off') +
    labs(
      x = "Chromosome",
      y = expression(-log[10](italic(p)))
    )
  
  ## add gene labels to index variants
  if (annotate_loci & any(sst$p_value < 5e-8)) {
    
    if (shadow) {
      p <- p +
        shadowtext::geom_shadowtext(
          aes(label = ifelse(p_value < loci_p_value & !is.na(gene), gene, "")),
          angle = 90,
          hjust = -0.2,
          size = 1.7,
          colour = "#bd0026",
          bg.colour = "white",
          bg.r = 0.2
        )
    } else {
      
     p <- p +
       geom_text(
         aes(label = ifelse(p_value < loci_p_value & !is.na(gene), gene, "")),
         angle = 90,
         hjust = -0.3,
         size = 1.7,
         colour = "#bd0026"
       )

    }
  }
  
  return(p)
  
}


