snp_modifyBuild_local <- function (
    info_snp,
    liftOver,
    from = "hg18",
    to = "hg19",
    check_reverse = TRUE)
  {
  ### creating a local version if Prive's snp_modifyBuild() function
  ### because somehow it keeps failing downloading the chain files

    if (!all(c("chr", "pos") %in% names(info_snp)))
      stop2("Expecting variables 'chr' and 'pos' in input 'info_snp'.")
    liftOver <-normalizePath(liftOver)
    info_BED <- with(info_snp,
                     data.frame(
                       chrom = paste0("chr",
                                      sub("^0", "", chr)),
                       start = pos - 1L,
                       end = pos,
                       id = seq_along(pos)
                     ))
    BED <- tempfile(fileext = ".BED")
    bigreadr::fwrite2(
      stats::na.omit(info_BED),
      BED,
      col.names = FALSE,
      sep = " ",
      scipen = 50
    )
    chain <- here("data", "reference_data", "chains", paste0(from, "To", tools::toTitleCase(to), ".over.chain.gz"))
    lifted <- tempfile(fileext = ".BED")
    system2(liftOver, c(BED, chain, lifted, tempfile(fileext = ".txt")))
    new_pos <- bigreadr::fread2(lifted, nThread = 1)
    is_bad <-
      vctrs::vec_duplicate_detect(new_pos$V4) | (new_pos$V1 !=
                                                   info_BED$chrom[new_pos$V4])
    new_pos <- new_pos[which(!is_bad),]
    pos0 <- info_snp$pos
    info_snp$pos <- NA_integer_
    info_snp$pos[new_pos$V4] <- new_pos$V3
    if (check_reverse) {
      pos2 <- suppressMessages(Recall(
        info_snp,
        liftOver,
        from = to,
        to = from,
        check_reverse = FALSE
      )$pos)
      info_snp$pos[pos2 != pos0] <- NA_integer_
    }
    bigassertr::message2("%d variants have not been mapped.", sum(is.na(info_snp$pos)))
    info_snp
  }