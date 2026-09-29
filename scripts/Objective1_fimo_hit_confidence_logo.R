## Visualize FIMO hits for one motif as:
##  (A) a threshold table -- how many hits survive at each p-value cutoff
##  (B) the HOCOMOCO reference motif, rendered as a real bits-scaled MEME-style
##      logo (via motifStack, same engine as the existing PPRE reference panel)
##  (C) a "sequence ladder" -- the most frequent matched sequences, drawn as
##      colored-letter rows and placed on the y-axis by their mean FIMO score,
##      so more-confident variants sit higher than looser ones.
##
## Answers the WhatsApp thread: "is there a chance to see the hits based on
## how far they are on the vertical axis, while considering the sequence?",
## "how many hits would fit where on such matrices?", and "can you give the
## sequence how it's reported online for this motif_id to compare".

library(dplyr)
library(readr)
library(ggplot2)
library(scales)
library(gridExtra)
library(motifStack)
library(grid)

fimo_dir <- "data/PPAR_RXR_motifs_Hocomoco_Jaspar/FIMO"
hocomoco_meme_dir <- "data/PPAR_RXR_motifs_Hocomoco_Jaspar/HOCOMOCO"
motif_id <- "PPARG.H13CORE.0.P.B"   # change to any motif_id present in the FIMO tsvs
top_n <- 12                          # how many distinct sequence variants to show
p_thresholds <- c(1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10)

nuc_colors <- c(A = "#008000", C = "#0000CC", G = "#FFA500", T = "#CC0000")

## ---- parse a MEME-format letter-probability matrix into a 4 x width matrix ----
parse_meme_pfm <- function(path) {
  lines <- readLines(path)
  header_idx <- grep("^letter-probability matrix", lines)
  stopifnot(length(header_idx) == 1)
  w <- as.integer(sub(".*w=\\s*(\\d+).*", "\\1", lines[header_idx]))
  mat_lines <- lines[(header_idx + 1):(header_idx + w)]
  mat <- t(do.call(rbind, lapply(mat_lines, function(l) as.numeric(strsplit(trimws(l), "\\s+")[[1]]))))
  rownames(mat) <- c("A", "C", "G", "T")   # matches ALPHABET= ACGT declaration
  mat
}

## ---- load & filter hits for this motif across all FIMO tsv shards ----
fimo_files <- list.files(fimo_dir, pattern = "\\.tsv$", full.names = TRUE)

fimo_col_types <- cols(
  motif_id = col_character(), motif_alt_id = col_character(),
  sequence_name = col_character(), start = col_double(), stop = col_double(),
  strand = col_character(), score = col_double(), `p-value` = col_double(),
  `q-value` = col_double(), matched_sequence = col_character()
)

hits <- lapply(fimo_files, function(f) {
  read_tsv(f, comment = "#", col_types = fimo_col_types) %>%
    filter(motif_id == !!motif_id)
}) %>% bind_rows()

stopifnot(nrow(hits) > 0)
hits <- hits %>% mutate(matched_sequence = toupper(matched_sequence))
motif_width <- nchar(hits$matched_sequence[1])

## ---- (A) threshold table: how many hits survive at each p-value cutoff ----
threshold_tbl <- tibble(p_value = p_thresholds) %>%
  rowwise() %>%
  mutate(n_hits = sum(hits$`p-value` <= p_value)) %>%
  ungroup()

p_threshold <- ggplot(threshold_tbl, aes(x = factor(p_value, levels = p_thresholds), y = n_hits)) +
  geom_col(fill = "grey40") +
  geom_text(aes(label = comma(n_hits)), vjust = -0.4, size = 3.2) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15))) +
  labs(x = "p-value threshold", y = "number of FIMO hits",
       title = paste0(motif_id, " -- hits surviving each p-value cutoff"),
       subtitle = paste0("total hits at FIMO run threshold: ", comma(nrow(hits)))) +
  theme_minimal(base_size = 11)

## ---- (B) HOCOMOCO reference motif as a real bits-scaled logo ----
hocomoco_meme_file <- file.path(hocomoco_meme_dir, paste0(motif_id, "_meme_format.meme"))
ref_logo_grob <- NULL
if (file.exists(hocomoco_meme_file)) {
  ref_mat <- parse_meme_pfm(hocomoco_meme_file)
  ref_pfm <- new("pfm", mat = ref_mat, name = paste0(motif_id, " (HOCOMOCO reference)"))
  ref_logo_grob <- grid::grid.grabExpr({
    plotMotifLogo(pfm = ref_pfm, motifName = paste0(motif_id, " -- HOCOMOCO reference (online record)"),
                  p = ref_pfm@background, colset = ref_pfm@color,
                  yaxis = TRUE, ylab = "Information content", xlab = NA,
                  margins = c(1.5, 4.1, 3.0, 1.0), newpage = FALSE)
  })
} else {
  message("No HOCOMOCO meme file found at ", hocomoco_meme_file, " -- skipping reference panel.")
}

## ---- (C) sequence ladder: top variants placed by mean score ----
variant_tbl <- hits %>%
  group_by(matched_sequence) %>%
  summarise(n = n(), mean_score = mean(score), mean_neglog10p = mean(-log10(`p-value`)), .groups = "drop") %>%
  arrange(desc(n)) %>%
  slice_head(n = top_n) %>%
  arrange(desc(mean_score)) %>%
  mutate(row_y = row_number())   # evenly spaced rows, ordered by confidence (mean score)

letters_df <- variant_tbl %>%
  mutate(row = row_number()) %>%
  tidyr::separate_rows(matched_sequence, sep = "(?<=.)(?=.)") %>%
  group_by(row) %>%
  mutate(position = row_number()) %>%
  ungroup() %>%
  rename(base = matched_sequence)

p_ladder <- ggplot(letters_df, aes(x = position, y = -row_y, group = row)) +
  geom_text(aes(label = base, color = base), fontface = "bold", size = 5,
            family = "mono") +
  geom_text(data = variant_tbl, aes(x = motif_width + 1.3, y = -row_y,
                                     label = sprintf("score=%.1f, n=%s", mean_score, comma(n))),
            inherit.aes = FALSE, hjust = 0, size = 3, color = "grey30") +
  scale_color_manual(values = nuc_colors, guide = "none") +
  scale_x_continuous(limits = c(0.3, motif_width + 5), breaks = seq_len(motif_width)) +
  labs(x = "position", y = NULL,
       title = paste0(motif_id, " -- top ", top_n, " matched sequences by frequency"),
       subtitle = "rows ranked top-to-bottom by mean FIMO score (most confident sequence on top)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        axis.text.y = element_blank())

## ---- (D) rank-abundance plot: full distribution across all distinct sequences ----
variant_tbl_full <- hits %>%
  group_by(matched_sequence) %>%
  summarise(n = n(), mean_score = mean(score), .groups = "drop") %>%
  arrange(desc(n)) %>%
  mutate(rank = row_number())

p_rank_abundance <- ggplot(variant_tbl_full, aes(x = rank, y = n)) +
  geom_line(color = "grey30", linewidth = 0.6) +
  scale_x_log10(labels = comma) +
  scale_y_log10(labels = comma) +
  annotation_logticks(sides = "bl") +
  labs(x = "sequence rank (log scale)", y = "hit count for that exact sequence (log scale)",
       title = paste0(motif_id, " -- rank-abundance across all ", comma(nrow(variant_tbl_full)), " distinct sequences"),
       subtitle = paste0(comma(nrow(hits)), " total hits; steep drop-off = a few sequences dominate, long flat tail = many one-off variants")) +
  theme_minimal(base_size = 11)

## ---- (E) score histogram, colored by how many distinct sequences fall in each bin ----
score_binwidth <- 0.5
score_bins <- hits %>%
  mutate(score_bin = floor(score / score_binwidth) * score_binwidth) %>%
  group_by(score_bin) %>%
  summarise(n_hits = n(), n_distinct_seqs = n_distinct(matched_sequence), .groups = "drop")

p_score_hist <- ggplot(score_bins, aes(x = score_bin, y = n_hits, fill = n_distinct_seqs)) +
  geom_col(width = score_binwidth * 0.9) +
  scale_fill_viridis_c(name = "distinct\nsequences", labels = comma) +
  scale_y_continuous(labels = comma) +
  labs(x = "FIMO score (bin width 0.5)", y = "number of hits",
       title = paste0(motif_id, " -- score distribution, colored by sequence diversity per bin"),
       subtitle = "dark/low bins = few distinct sequences repeating a lot; bright/high bins = many different sequences, each rare") +
  theme_minimal(base_size = 11)

rank_hist_path <- file.path(dirname(fimo_dir), "Objective1_results", paste0(motif_id, "_rank_abundance_and_score_histogram.png"))
rank_hist_combined <- arrangeGrob(grobs = list(ggplotGrob(p_rank_abundance), ggplotGrob(p_score_hist)), ncol = 1, heights = c(1, 1))
ggsave(rank_hist_path, rank_hist_combined, width = 9, height = 8, dpi = 150, bg = "white")
message("Saved: ", rank_hist_path)

## ---- combine and save ----
out_path <- file.path(dirname(fimo_dir), "Objective1_results", paste0(motif_id, "_hit_confidence_ladder.png"))
dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)

panels <- list(ggplotGrob(p_threshold))
heights <- c(1)
if (!is.null(ref_logo_grob)) {
  panels <- c(panels, list(ref_logo_grob))
  heights <- c(heights, 0.8)
}
panels <- c(panels, list(ggplotGrob(p_ladder)))
heights <- c(heights, 1.4)

combined <- arrangeGrob(grobs = panels, ncol = 1, heights = heights)
ggsave(out_path, combined, width = 9, height = sum(heights) * 3.4, dpi = 150, bg = "white")
message("Saved: ", out_path)
