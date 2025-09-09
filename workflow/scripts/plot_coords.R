library(ggplot2)
library(viridis)
library(here)

args <- commandArgs(trailingOnly = TRUE)
(in_dir <- args[1])
(out_dir <- args[2])
(reference <- args[3])
(query <- args[4])

#Uncomment below for testing
#in_dir <- "filter20kb_coords"
#out_dir <- "plots"
#reference <- "Chinese Spring"
#query <- "Hilliard"
#dir.create(out_dir)

chrs <- c("1A", "1B", "1D", 
          "2A", "2B", "2D", 
          "3A", "3B", "3D", 
          "4A", "4B", "4D", 
          "5A", "5B", "5D", 
          "6A", "6B", "6D", 
          "7A", "7B", "7D")

filenames <- list.files(in_dir)

hdr <- readLines(here(in_dir, filenames[1]), n = 20)
skip <- which(grepl("^=+", hdr))[1] + 1

list_dfs <- lapply(here(in_dir, filenames), 
                   read.table, 
                   sep = "|",
                   skip = skip,
                   header = FALSE,
                   fill = TRUE,
                   strip.white = TRUE,
                   comment.char = "",
                   quote = "")

split_cols <- function(x, n, names) {
  tmp <- strsplit(trimws(x), "\\s+")
  # pad short rows
  tmp <- lapply(tmp, function(v) { length(v) <- n; v })
  out <- do.call(rbind, tmp)
  out <- as.data.frame(out, stringsAsFactors = FALSE)
  names(out) <- names
  out[] <- lapply(out, type.convert, as.is = TRUE)
  out
}

coords_list <- list()

for (i in 1:length(list_dfs)) {
  coords_raw <- list_dfs[[i]]
  g1 <- split_cols(coords_raw[[1]], 2, c("S1","E1"))
  g2 <- split_cols(coords_raw[[2]], 2, c("S2","E2"))
  g3 <- split_cols(coords_raw[[3]], 2, c("LEN1","LEN2"))
  # Some coords files show %IDY as a single value; if yours has two, change '1'->'2' and names accordingly
  g4 <- split_cols(coords_raw[[4]], 1, c("PctIDY"))
  g5 <- split_cols(coords_raw[[5]], 2, c("LenR","LenQ"))
  g6 <- split_cols(coords_raw[[6]], 2, c("CovR","CovQ"))
  
  # Tags column often contains two names separated by a tab
  tags <- do.call(rbind, strsplit(trimws(coords_raw[[7]]), "\t"))
  tags <- as.data.frame(tags, stringsAsFactors = FALSE)
  names(tags) <- c("TagR","TagQ")
  
  coords_list[[i]] <- cbind(g1,g2,g3,g4,g5,g6,tags)
}

names(coords_list) <- filenames

for (i in 1:length(coords_list)) {
  coords <- coords_list[[i]]
  #Remove alignments from different chr
  coords_chr <- coords[coords$TagR == coords$TagQ,] 
  
  p1 <- ggplot(coords_chr, aes(x = S1 / 1000000, y = PctIDY)) + 
    geom_point(shape = 1) +
    labs(title = chrs[i], 
         subtitle = paste0(query, " aligned to ", reference), 
         x = "Mbp", 
         y = "% ID")
  ggsave(paste0(out_dir, "/seqid_", names(coords_list[i]), ".jpeg"), p1)
  
  coords_floor <- coords_chr[coords_chr$PctIDY >= 97,]
  
  p2 <- ggplot(coords_floor, aes(x = S1 / 1000000, y = S2 / 1000000, color = PctIDY)) + 
    geom_point(alpha=0.5) +
    scale_color_viridis() +
    labs(x = paste0(reference, " (Mbp)"),
         y = paste0(query, " (Mbp)"), 
         color = "% ID", 
         title = chrs[i], 
         subtitle = paste0(query, " aligned to ", reference))
  ggsave(paste0(out_dir, "/pairwise_", names(coords_list[i]), ".jpeg"), p2)
}
