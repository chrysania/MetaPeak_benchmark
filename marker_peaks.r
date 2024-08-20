### FUNCTION: get_marker_peaks
###           get top n foldchange marker peaks for a cell type
###           manually run step 1-2 first to determine chromosome
###           have ccre_gr and cpeaks_gr granges ready
ccre <- read.table("data/GRCh38-cCREs.bed", sep = "\t", header = FALSE)
colnames(ccre) <- c("chr", "start", "stop", "ident", "ident2", "class")
ccre_gr <- GRanges(
  seqnames = ccre$chr,
  ranges = IRanges(start = ccre$start, end = ccre$stop)
)

cpeaks <- read.table("data/combined_cre.bed", sep = "\t", header = FALSE)
colnames(cpeaks) <- c("chr", "start", "stop", "ID", "class")
cpeaks_gr <- GRanges(
  seqnames = cpeaks$chr,
  ranges = IRanges(start = cpeaks$start, end = cpeaks$stop)
)
### params: seurat_obj: seurat object with cell annotation
###         cell_type: cell ident.1 for calculating FoldChange
###         chromosome: chromosome name string
###         n_candidates: number of foldchange peaks to be checked for <1MB distance pairing

get_marker_peaks <- function(seurat_obj, cell_type, chromosome, n_candidates) {
    # 1. calculate fold change for cell_type
    fc <- FoldChange(seurat_obj, ident.1 = cell_type)
    fc <- fc[order(fc$avg_log2FC, decreasing = TRUE), ]

    # 2. get top n fc for candidates
    fc_chr <- fc[grep(paste0("^", chromosome, "-"), rownames(fc)), ]
    chr_marker_candidates <- rownames(head(fc_chr, n_candidates))

    # 3. create fc dataframe
    split_info <- strsplit(chr_marker_candidates, "-")
    split_info <- do.call(rbind, split_info)

    chr_df <- data.frame(
      chr = split_info[, 1],
      start = as.integer(split_info[, 2]),
      stop = as.integer(split_info[, 3]),
      rowname = chr_marker_candidates,
      midpoint = (as.integer(split_info[, 2]) + as.integer(split_info[, 3])) / 2
    )
    chr_gr <- GRanges(seqnames = chr_df$chr, ranges = IRanges(start = chr_df$start, end = chr_df$stop)) # get granges obj

    # 4. generate pairs of peaks within 1MB
    row_indices <- combn(nrow(chr_df), 2)
    valid_pairs <- list()
    for (i in 1:ncol(row_indices)) {
        idx1 <- row_indices[1, i]
        idx2 <- row_indices[2, i]
        if (abs(chr_df$midpoint[idx1] - chr_df$midpoint[idx2]) < 1000000) {
            valid_pairs[[length(valid_pairs) + 1]] <- c(idx1, idx2)
        }
    }

    if (length(valid_pairs) > 0) {
        valid_pairs <- do.call(rbind, valid_pairs)
        pairs_df <- data.frame(idx1 = valid_pairs[, 1],
                               idx2 = valid_pairs[, 2],
                               mid1 = chr_df$midpoint[valid_pairs[, 1]],
                               mid2 = chr_df$midpoint[valid_pairs[, 2]],
                               peak1 = chr_df$rowname[valid_pairs[, 1]],
                               peak2 = chr_df$rowname[valid_pairs[, 2]])
    } else {
        message("No pairs found with distance less than 1MB")
        return(NA)
    }

    # 5. check overlaps with ccre, cpeaks+ccre
    overlaps_ccre <- findOverlaps(chr_gr, ccre_gr, type="any")
    overlaps_ccre_df <- as.data.frame(overlaps_ccre)

    overlaps_cpeaks <- findOverlaps(chr_gr, cpeaks_gr, type="any")
    overlaps_cpeaks_df <- as.data.frame(overlaps_cpeaks)

    # 6. match seurat obj peaks with ccre, cpeaks+ccre
    pairs_df$ccre1 <- ""
    pairs_df$ccre2 <- ""
    pairs_df$cpeaks1 <- ""
    pairs_df$cpeaks2 <- ""

    for (i in 1:nrow(pairs_df)) {  # ccre
        idx1 <- pairs_df$idx1[i]
        idx2 <- pairs_df$idx2[i]
  
        peakset1 <- paste(overlaps_ccre_df$subjectHits[overlaps_ccre_df$queryHits == idx1], collapse = ",")
        peakset2 <- paste(overlaps_ccre_df$subjectHits[overlaps_ccre_df$queryHits == idx2], collapse = ",")
  
        if (nchar(peakset1) > 0) {
            pairs_df$ccre1[i] <- peakset1
        }
        if (nchar(peakset2) > 0) {
            pairs_df$ccre2[i] <- peakset2
        }
    }

    for (i in 1:nrow(pairs_df)) {  # cpeaks+ccre
        idx1 <- pairs_df$idx1[i]
        idx2 <- pairs_df$idx2[i]
  
        peakset1 <- paste(overlaps_cpeaks_df$subjectHits[overlaps_cpeaks_df$queryHits == idx1], collapse = ",")
        peakset2 <- paste(overlaps_cpeaks_df$subjectHits[overlaps_cpeaks_df$queryHits == idx2], collapse = ",")
  
        if (nchar(peakset1) > 0) {
            pairs_df$cpeaks1[i] <- peakset1
        }
        if (nchar(peakset2) > 0) {
            pairs_df$cpeaks2[i] <- peakset2
        }
    }

    # 7. Filter out peaks with no overlaps with ccre and/or ccre+cpeaks
    pairs_df <- pairs_df[!(pairs_df$ccre1 == "" | pairs_df$ccre2 == "" | pairs_df$cpeaks1 == "" | pairs_df$cpeaks2 == ""), ]

    return(pairs_df)
}