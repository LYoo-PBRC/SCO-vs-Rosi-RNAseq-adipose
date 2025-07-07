#######################################################
# SECTION 1: Setup - Load Required Libraries and Create Output Folders
#######################################################
library(edgeR)
library(limma)
library(ggplot2)
library(dplyr)
library(openxlsx)
library(reshape2)
library(readr)
library(ggrepel)
library(readxl)
library(enrichR)

# Custom QC/DE plot colors
qc_colors <- c("DMSO"="#00c0c0", "SCO"="#00E200", "CAP7"="#7f8c8d", "ROSI"="#ff6666")

output_folder <- "output"
subfolders <- c("QC_plots", "limma_DE_results", "limma_Venn", "limma_Reversal", "Pathway_enrichment", "fgsea")
for (f in subfolders) dir.create(file.path(output_folder,f), recursive=TRUE, showWarnings=FALSE)

#######################################################
# SECTION 2: Load Count Data and Setup Gene Names
#######################################################
count_data  <- read.xlsx("annotated_counts.xlsx")
sample_cols <- paste0("Sample_", sprintf("%02d", 1:32))

gene_col <- intersect(c("gene_name", "GeneName", "Gene", "gene"), colnames(count_data))[1]
if (is.na(gene_col)) gene_col <- colnames(count_data)[1]
possible_ensembl_cols <- c("ensembl_gene_id", "Ensembl_ID", "ensembl_id", "gene_id", "GeneID", "Gene_Id")
ensembl_col <- intersect(possible_ensembl_cols, colnames(count_data))[1]
if (is.na(ensembl_col)) stop("No Ensembl-ID column found in annotated_counts.xlsx")
gene_names <- as.character(count_data[[gene_col]])
gene_names[gene_names == ""] <- NA
gene_names[is.na(gene_names)] <- count_data[[ensembl_col]][is.na(gene_names)]
dup_idx <- duplicated(gene_names)
iter <- 1
while (any(dup_idx)) {
  gene_names[dup_idx] <- paste0(count_data[[ensembl_col]][dup_idx], "_", gene_names[dup_idx])
  dup_idx <- duplicated(gene_names)
  iter <- iter + 1
  if (iter > 5) stop("Still have duplicate names after 5 iterations")
}
rownames(count_data) <- gene_names
counts <- as.matrix(count_data[ , sample_cols ])
mode(counts) <- "integer"
annotation_cols <- setdiff(colnames(count_data), sample_cols)
annotations     <- count_data[ , annotation_cols, drop = FALSE]
annotations$GeneName <- rownames(count_data)

#######################################################
# SECTION 3: Setup Metadata for Samples
#######################################################
condition <- rep(c("Basal","TNF"), each = 16)
treatment <- rep(rep(c("DMSO","SCO","CAP7","ROSI"), each = 4), times = 2)
coldata <- data.frame(
  condition = factor(condition),
  treatment = factor(treatment, levels = c("DMSO","SCO","CAP7","ROSI")),
  row.names = colnames(counts)
)
coldata$group <- factor(
  paste(coldata$condition, coldata$treatment, sep = "_"),
  levels = c("Basal_DMSO","Basal_SCO","Basal_CAP7","Basal_ROSI",
             "TNF_DMSO","TNF_SCO","TNF_CAP7","TNF_ROSI"),
  labels = c("B.DMSO","B.SCO","B.CAP7","B.ROSI",
             "T.DMSO","T.SCO","T.CAP7","T.ROSI")
)

#######################################################
# SECTION 4: edgeR Normalization then Filtering
#######################################################
y <- DGEList(counts = counts, group = coldata$group)
y <- calcNormFactors(y)
keep <- rowSums(cpm(y) > 1) >= 16
y <- y[keep , , keep.lib.sizes = FALSE]
cat("Genes after filter :", nrow(y$counts), "\n")

#######################################################
# SECTION 5: Quality Control (QC) Plots 
#######################################################
logCPM <- cpm(y, log = TRUE, prior.count = 3)
df_den <- reshape2::melt(logCPM, varnames = c("Gene", "Sample"), value.name = "logCPM")
df_den$Treatment <- coldata$treatment[match(df_den$Sample, rownames(coldata))]
p_den <- ggplot(df_den, aes(x = logCPM, color = Treatment)) +
  geom_density(alpha = 0.7, size = 1.2) + theme_minimal() +
  labs(title = "Density Plot (logCPM)") +
  scale_color_manual(values = qc_colors) +
  theme(legend.title = element_blank())
ggsave(file.path(output_folder,"QC_plots","Density_logCPM.png"), p_den, width = 8, height = 6, dpi = 600)

png(file.path(output_folder,"QC_plots","MDS.png"), width = 2000, height = 1600, res = 300)
plotMDS(y, col = qc_colors[as.character(coldata$treatment)], pch = as.numeric(coldata$condition) + 15)
legend("topright", levels(coldata$treatment), col = qc_colors, pch = 16, bty = "n")
dev.off()

pca <- prcomp(t(logCPM), scale. = TRUE)
pvar <- round((pca$sdev^2)/sum(pca$sdev^2)*100,1)
pca_df <- data.frame(
  PC1 = pca$x[,1], PC2 = pca$x[,2],
  Treatment = coldata$treatment, Condition = coldata$condition
)
p_pca <- ggplot(pca_df, aes(PC1, PC2, color = Treatment, shape = Condition)) +
  geom_point(size = 4, alpha = 0.8) + theme_minimal() +
  labs(title = "PCA Plot (logCPM)", x = paste0("PC1 (", pvar[1], "%)"), y = paste0("PC2 (", pvar[2], "%)")) +
  scale_color_manual(values = qc_colors)
ggsave(file.path(output_folder,"QC_plots","PCA.png"), p_pca, width = 8, height = 6, dpi = 600)

#######################################################
# SECTION 6: Differential Expression Analysis (limma-voom)
#######################################################
design <- model.matrix(~0 + group, data = coldata)
colnames(design) <- levels(coldata$group)
v <- voom(y, design)
fit <- lmFit(v, design)
CM <- makeContrasts(
  B.SCO  = B.SCO - B.DMSO,
  B.ROSI = B.ROSI - B.DMSO,
  B.CAP7 = B.CAP7 - B.DMSO,
  T.DMSO = T.DMSO - B.DMSO,
  T.SCO  = T.SCO  - T.DMSO,
  T.ROSI = T.ROSI - T.DMSO,
  T.CAP7 = T.CAP7 - T.DMSO,
  levels = design
)
fit2 <- eBayes(contrasts.fit(fit, CM))

#######################################################
# SECTION 7: Write DE Tables 
#######################################################
sig_cut  <- 0.05
limma_res <- list()
for (cn in colnames(CM)) {
  res_full <- topTable(fit2, coef = cn, n = Inf, adjust.method = "BH")
  if (nrow(res_full) == 0) next
  res_full$GeneID    <- rownames(res_full)
  res_full$gene_name <- annotations$GeneName[match(rownames(res_full), rownames(annotations))]
  res_full <- res_full[, c("GeneID","gene_name", setdiff(colnames(res_full), c("GeneID","gene_name")))]
  res_sig <- res_full[ res_full$adj.P.Val < sig_cut , ]
  # Export ONLY significant results:
  write.csv(
    res_sig,
    file.path(output_folder, "limma_DE_results",
              sprintf("limma_DE_%s_sig_%.2g.csv", cn, sig_cut)),
    row.names = FALSE
  )
  limma_res[[cn]] <- res_sig
  message(sprintf(
    "Contrast %s ??? %d significant genes (adj.P.Val < %.2g)",
    cn, nrow(res_sig), sig_cut))
}

#######################################################
# SECTION 8: Reversal Analysis (CSV outputs)
#######################################################
tnf_threshold <- 0.05
treatment_threshold <- 0.05

limma_TNF  <- topTable(fit2, coef="T.DMSO", n=Inf)
limma_ROSI <- topTable(fit2, coef="T.ROSI", n=Inf)
limma_SCO  <- topTable(fit2, coef="T.SCO",  n=Inf)

sig_TNF  <- limma_TNF[limma_TNF$adj.P.Val < tnf_threshold, ]
sig_ROSI <- limma_ROSI[limma_ROSI$adj.P.Val < treatment_threshold, ]
sig_SCO  <- limma_SCO[limma_SCO$adj.P.Val < treatment_threshold, ]

shared_ROSI <- intersect(rownames(sig_TNF), rownames(sig_ROSI))
merged_ROSI <- merge(sig_TNF[shared_ROSI,], sig_ROSI[shared_ROSI,], by="row.names", suffix=c("_TNF","_ROSI"))
reversed_ROSI <- merged_ROSI[
  (merged_ROSI$logFC_TNF > 0 & merged_ROSI$logFC_ROSI < 0) |
    (merged_ROSI$logFC_TNF < 0 & merged_ROSI$logFC_ROSI > 0),
  "Row.names"]

shared_SCO <- intersect(rownames(sig_TNF), rownames(sig_SCO))
merged_SCO <- merge(sig_TNF[shared_SCO,], sig_SCO[shared_SCO,], by="row.names", suffix=c("_TNF","_SCO"))
reversed_SCO <- merged_SCO[
  (merged_SCO$logFC_TNF > 0 & merged_SCO$logFC_SCO < 0) |
    (merged_SCO$logFC_TNF < 0 & merged_SCO$logFC_SCO > 0),
  "Row.names"]

sco_only   <- setdiff(reversed_SCO, reversed_ROSI)
rosi_only  <- setdiff(reversed_ROSI, reversed_SCO)
shared_rev <- intersect(reversed_SCO, reversed_ROSI)

write.csv(merged_SCO,  file.path(output_folder, "limma_Reversal", "SCO_reversal_genes.csv"), row.names=FALSE)
write.csv(merged_ROSI, file.path(output_folder, "limma_Reversal", "ROSI_reversal_genes.csv"), row.names=FALSE)

# --- Directionality stacked barplot for reversal genes ---
opposite_direction_SCO <- sum(sign(merged_SCO$logFC_TNF) != sign(merged_SCO$logFC_SCO))
same_direction_SCO <- sum(sign(merged_SCO$logFC_TNF) == sign(merged_SCO$logFC_SCO))
opposite_direction_ROSI <- sum(sign(merged_ROSI$logFC_TNF) != sign(merged_ROSI$logFC_ROSI))
same_direction_ROSI <- sum(sign(merged_ROSI$logFC_TNF) == sign(merged_ROSI$logFC_ROSI))

legend_levels <- c(
  "ROSI_Opposite Direction",
  "ROSI_Same Direction",
  "SCO_Opposite Direction",
  "SCO_Same Direction"
)
bar_data <- data.frame(
  Treatment = rep(c("ROSI", "SCO"), each = 2),
  Direction = rep(c("Opposite Direction", "Same Direction"), 2),
  Count = c(opposite_direction_ROSI, same_direction_ROSI,
            opposite_direction_SCO, same_direction_SCO)
)
bar_data$Group <- interaction(bar_data$Treatment, bar_data$Direction, sep = "_")
bar_data$Group <- factor(bar_data$Group, levels = legend_levels)

bar_plot <- ggplot(bar_data, aes(x = Treatment, y = Count, fill = Group)) +
  geom_bar(
    stat = "identity",
    position = "stack",
    color = "black",
    width = 0.4
  ) +
  scale_fill_manual(
    values = c(
      "ROSI_Opposite Direction" = "#ff6666",
      "ROSI_Same Direction"     = "#FFC1C1",
      "SCO_Opposite Direction"  = "#00E200",
      "SCO_Same Direction"      = "#BDFCC9"
    ),
    labels = c("", "", "", ""),
    name = NULL
  ) +
  scale_x_discrete(expand = expansion(add = 0.36), labels = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme_minimal(base_size = 15) +
  theme(
    axis.text.x = element_blank(),
    plot.title = element_blank(),
    legend.title = element_blank(),
    legend.text = element_blank()
  ) +
  guides(fill = guide_legend(order = 1))

ggsave(
  file.path(output_folder, "limma_Reversal", "stacked_bar_overlap_directionality.png"),
  bar_plot,
  width = 3.5,
  height = 6,
  dpi = 400
)


#######################################################
# SECTION 9: Venn CSVs Only (2-way and 3-way)
#######################################################
get_sig_genes <- function(df, pcol, fcol) {
  sig <- df[df[[pcol]] < 0.05, ]
  if ("gene_name" %in% colnames(sig)) sig$gene_name else rownames(sig)
}
lSets <- list(
  B.SCO  = get_sig_genes(limma_res[["B.SCO"]],  "adj.P.Val", "logFC"),
  B.ROSI = get_sig_genes(limma_res[["B.ROSI"]], "adj.P.Val", "logFC"),
  B.CAP7 = get_sig_genes(limma_res[["B.CAP7"]], "adj.P.Val", "logFC"),
  T.DMSO = get_sig_genes(limma_res[["T.DMSO"]], "adj.P.Val", "logFC"),
  T.SCO  = get_sig_genes(limma_res[["T.SCO"]],  "adj.P.Val", "logFC"),
  T.ROSI = get_sig_genes(limma_res[["T.ROSI"]], "adj.P.Val", "logFC")
)

# 2-way venns
pairs <- list(
  c("B.SCO","B.ROSI"), c("B.SCO","B.CAP7"), c("B.CAP7","B.ROSI"),
  c("T.DMSO","T.SCO"), c("T.DMSO","T.ROSI"),
  c("T.SCO","B.SCO"),  c("T.ROSI","B.ROSI")
)
plot_venn_export_csv <- function(sets, csv_file) {
  n <- length(sets)
  names(sets) <- make.names(names(sets))
  if (n == 2) {
    a <- sets[[1]]; b <- sets[[2]]
    shared <- intersect(a, b)
    only_a <- setdiff(a, b)
    only_b <- setdiff(b, a)
    max_len <- max(length(shared), length(only_a), length(only_b))
    df <- data.frame(
      Shared = c(shared, rep(NA, max_len - length(shared))),
      Only_A = c(only_a, rep(NA, max_len - length(only_a))),
      Only_B = c(only_b, rep(NA, max_len - length(only_b)))
    )
    colnames(df) <- c("Shared", names(sets)[1], names(sets)[2])
    write.csv(df, csv_file, row.names = FALSE)
  }
}
for (pr in pairs) {
  s1 <- pr[1]; s2 <- pr[2]
  sets <- list(lSets[[s1]], lSets[[s2]])
  names(sets) <- c(s1, s2)
  csvf <- file.path(output_folder, "limma_Venn", paste0("Venn_", s1, "__", s2, ".csv"))
  plot_venn_export_csv(sets, csvf)
}

# 3-way venn for T.DMSO, T.ROSI, T.SCO (mutually exclusive slices, wide csv)
three_sets <- list(
  T.DMSO = lSets[["T.DMSO"]],
  T.ROSI = lSets[["T.ROSI"]],
  T.SCO  = lSets[["T.SCO"]]
)
A <- three_sets[[1]]; B <- three_sets[[2]]; C <- three_sets[[3]]
venn_regions <- list(
  only_T.DMSO       = setdiff(A, union(B, C)),
  only_T.ROSI       = setdiff(B, union(A, C)),
  only_T.SCO        = setdiff(C, union(A, B)),
  T.DMSO_and_T.ROSI = setdiff(intersect(A, B), C),
  T.DMSO_and_T.SCO  = setdiff(intersect(A, C), B),
  T.ROSI_and_T.SCO  = setdiff(intersect(B, C), A),
  all_three         = Reduce(intersect, list(A, B, C))
)
region_names <- names(venn_regions)
max_len <- max(sapply(venn_regions, length))
venn_wide <- as.data.frame(
  sapply(venn_regions, function(x) { x <- as.character(x); length(x) <- max_len; x })
)
colnames(venn_wide) <- region_names
venn_csv_file <- file.path(output_folder, "limma_Venn", "Venn_T_DMSO_ROSI_SCO.csv")
write.csv(venn_wide, venn_csv_file, row.names = FALSE)

#######################################################
# SECTION 10: FGSEA Barplots (WikiPathways only)
#######################################################
fgsea_outdir <- file.path(output_folder, "fgsea")
if (!dir.exists(fgsea_outdir)) dir.create(fgsea_outdir, recursive = TRUE)

fgsea_file <- "fgsea_sigs.xlsx"
sheetnames <- readxl::excel_sheets(fgsea_file)
print(sheetnames)

df_sco_dmso <- readxl::read_xlsx(fgsea_file, sheet = "SCOvsDMSO_p0.01") %>%
  filter(Category == "WikiPathways", padj < 0.01)
df_rosi_dmso <- readxl::read_xlsx(fgsea_file, sheet = "ROSIvsDMSO_p0.01") %>%
  filter(Category == "WikiPathways", padj < 0.01)
df_rosi_sco <- readxl::read_xlsx(fgsea_file, sheet = "ROSIvsSCO_p0.01") %>%
  filter(Category == "WikiPathways", padj < 0.01)

write.csv(df_sco_dmso,  file.path(fgsea_outdir, "WikiPathways_SCO_vs_DMSO_fgsea.csv"),  row.names = FALSE)
write.csv(df_rosi_dmso, file.path(fgsea_outdir, "WikiPathways_ROSI_vs_DMSO_fgsea.csv"), row.names = FALSE)
write.csv(df_rosi_sco,  file.path(fgsea_outdir, "WikiPathways_ROSI_vs_SCO_fgsea.csv"),  row.names = FALSE)

plot_fgsea_barplot_nolabel <- function(df, title, outfile, low_col, high_col) {
  df$pathway_clean <- gsub("^WP_", "", df$pathway)
  df$log10padj <- -log10(df$padj)
  df <- df %>% arrange(padj, desc(abs(NES)))
  df$pathway_clean <- factor(df$pathway_clean, levels = rev(df$pathway_clean))
  nes_min <- min(df$NES, na.rm = TRUE)
  nes_max <- max(df$NES, na.rm = TRUE)
  padding <- 0.25
  upper_limit <- max(nes_max + padding, 3)
  lower_limit <- min(nes_min - padding, -2)
  nes_limits <- c(lower_limit, upper_limit)
  p <- ggplot(df) +
    geom_col(aes(x = pathway_clean, y = NES, fill = log10padj),
             color = "black", width = 0.9) +
    coord_flip() +
    scale_x_discrete(expand = c(0,0), labels = NULL) +
    scale_y_continuous(limits = nes_limits, expand = c(0,0)) +
    scale_fill_gradient(
      low = low_col, high = high_col, name = NULL,
      guide = guide_colorbar(title.position = "left", title.hjust = 1, title.vjust = .85)
    ) +
    theme_minimal(base_size = 15) +
    theme(
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x  = element_text(size = 14),
      axis.title.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.title = element_blank(),
      legend.position = "bottom",
      legend.key.height = unit(1.4, "lines"),
      legend.key.width  = unit(1.8, "lines"),
      legend.text  = element_text(size=13),
      legend.title = element_text(size=15, face="bold"),
      legend.spacing = unit(2, "lines"),
      plot.margin = margin(10, 12, 45, 12)
    ) +
    labs(x = NULL, y = NULL, title = NULL)
  ggsave(outfile, p, width = 5, height = max(4, nrow(df)*0.4), dpi = 400)
  print(p)
}
plot_fgsea_barplot_nolabel(df_sco_dmso,  "", file.path(fgsea_outdir, "WikiPathways_SCO_vs_DMSO_fgsea_barplot_gapless.png"), "#eaffea", "#00E200")
plot_fgsea_barplot_nolabel(df_rosi_dmso, "", file.path(fgsea_outdir, "WikiPathways_ROSI_vs_DMSO_fgsea_barplot_gapless.png"), "#ffeaea", "#ff6666")
plot_fgsea_barplot_nolabel(df_rosi_sco,  "", file.path(fgsea_outdir, "WikiPathways_ROSI_vs_SCO_fgsea_barplot_gapless.png"), "#e7ffff", "#00c0c0")

#######################################################
# SECTION 11: ENRICHR Barplots for Reversal & Venn Slices (WikiPathways only)
#######################################################
strip_prefix <- function(x) sub("^ENSMUSG[0-9]+_", "", x)
sco_genes   <- strip_prefix(lSets[["B.SCO"]])
rosi_genes  <- strip_prefix(lSets[["B.ROSI"]])
sco_unique  <- setdiff(sco_genes,  rosi_genes)
rosi_unique <- setdiff(rosi_genes, sco_genes)
shared_genes<- intersect(sco_genes, rosi_genes)
gene_sets_basal <- list(
  SCO_unique   = sco_unique,
  ROSI_unique  = rosi_unique,
  Shared_basal = shared_genes
)

enr_SCO    <- enrichr(sco_only,   "WikiPathways_2024_Mouse")$WikiPathways_2024_Mouse
enr_ROSI   <- enrichr(rosi_only,  "WikiPathways_2024_Mouse")$WikiPathways_2024_Mouse
enr_Shared <- enrichr(shared_rev, "WikiPathways_2024_Mouse")$WikiPathways_2024_Mouse
outdir_reversal <- file.path(output_folder, "Pathway_enrichment", "Reversal")
if (!dir.exists(outdir_reversal)) dir.create(outdir_reversal, recursive = TRUE)

plot_enrichr_barplot_nolabel <- function(
    df, out_prefix, db_label, outdir,
    n_plot = 15,
    plot_width = 7,
    bar_height = 1.5,
    min_height = 2.2,
    bar_width = 0.9,
    legend_key_height = 1.4,
    legend_key_width  = 1.8,
    legend_title_size = 15,
    legend_text_size  = 13
) {
  padj_col <- "Adjusted.P.value"
  p_col    <- "P.value"
  odds_col <- "Odds.Ratio"
  df[[padj_col]] <- as.numeric(df[[padj_col]])
  df_sig <- df[df[[padj_col]] < 0.01 & !is.na(df[[padj_col]]), ]
  if (nrow(df_sig) == 0) return(invisible(NULL))
  if (out_prefix %in% c("Shared", "Shared_basal", "Shared_rev")) {
    top_df <- df_sig
  } else {
    top_df <- head(df_sig, n_plot)
  }
  top_df$Term <- factor(top_df$Term, levels = rev(top_df$Term))
  if (grepl("ROSI", out_prefix, ignore.case = TRUE)) {
    low_col <- "#ffeaea"; high_col <- "#ff6666"
  } else if (grepl("SCO", out_prefix, ignore.case = TRUE)) {
    low_col <- "#eaffea"; high_col <- "#00E200"
  } else if (grepl("Shared", out_prefix, ignore.case = TRUE)) {
    low_col <- "#e5e5e5"; high_col <- "#333333"
  } else {
    low_col <- "#AAAAAA"; high_col <- "#333333"
  }
  odds_min <- min(top_df[[odds_col]], na.rm = TRUE)
  odds_max <- max(top_df[[odds_col]], na.rm = TRUE)
  padding  <- 0.25
  upper_limit <- max(odds_max + padding, 6)
  lower_limit <- min(odds_min - padding, 0)
  odds_limits <- c(lower_limit, upper_limit)
  outfile <- file.path(outdir, paste0("Enrichr_", db_label, "_Barplot_", out_prefix, "_nolabel.png"))
  plot_height <- max(min_height, length(unique(top_df$Term)) * bar_height)
  p <- ggplot(top_df) +
    geom_col(
      aes(x = Term, y = !!as.name(odds_col), fill = -log10(!!as.name(p_col))),
      color = "black", width = bar_width
    ) +
    coord_flip() +
    scale_x_discrete(expand = c(0,0), labels = NULL) +
    scale_y_continuous(limits = odds_limits, expand = c(0,0)) +
    scale_fill_gradient(
      low = low_col, high = high_col, name = NULL,
      guide = guide_colorbar(
        title.position = "left", title.hjust = 1, title.vjust = 0.85
      )
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x  = element_text(size = 54, face = "bold"),
      axis.title.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.title = element_blank(),
      legend.position = "bottom",
      legend.key.height = unit(legend_key_height, "lines"),
      legend.key.width  = unit(legend_key_width, "lines"),
      legend.text  = element_text(size=legend_text_size),
      legend.title = element_text(size=legend_title_size, face="bold"),
      legend.spacing = unit(2, "lines"),
      plot.margin = margin(10, 12, 45, 12)
    )
  ggsave(outfile, p, width = plot_width, height = plot_height, dpi = 400)
  print(p)
}
plot_enrichr_barplot_nolabel(enr_SCO,    "SCO_only",   "WikiPathways2024", outdir_reversal)
plot_enrichr_barplot_nolabel(enr_ROSI,   "ROSI_only",  "WikiPathways2024", outdir_reversal)
plot_enrichr_barplot_nolabel(enr_Shared, "Shared_rev", "WikiPathways2024", outdir_reversal)

#######################################################
# SECTION 12: Sanity Check / Summary Report
#######################################################
limma_pcut    <- 0.05
rev_tnf_cut   <- tnf_threshold
rev_treat_cut <- treatment_threshold
contrast_desc <- c(
  B.SCO  = "DE genes in Basal SCO vs Basal DMSO",
  B.ROSI = "DE genes in Basal ROSI vs Basal DMSO",
  B.CAP7 = "DE genes in Basal CAP7 vs Basal DMSO",
  T.DMSO = "DE genes in TNF vs Basal (DMSO vehicle)",
  T.SCO  = "DE genes in TNF+SCO vs TNF+DMSO",
  T.ROSI = "DE genes in TNF+ROSI vs TNF+DMSO",
  T.CAP7 = "DE genes in TNF+CAP7 vs TNF+DMSO"
)
pair_desc <- function(a, b) sprintf("overlap of %s and %s", contrast_desc[a], contrast_desc[b])
lines <- c(
  "============================================",
  "           limma Analysis Summary            ",
  "============================================",
  sprintf("Genes in input matrix:        %5d", nrow(count_data)),
  sprintf("Genes after edgeR filtering:  %5d", nrow(y$counts)),
  sprintf("=== limma DE genes per contrast (adj.P.Val < %.2g) ===", limma_pcut)
)
for (cn in names(limma_res)) {
  nDE <- if (!is.null(limma_res[[cn]])) nrow(limma_res[[cn]]) else 0
  lines <- c(lines, sprintf("  %-8s : %5d (%s)", cn, nDE, contrast_desc[cn]))
}
lines <- c(lines, "", "=== Pairwise Venn slices (adj.P.Val < 0.05) ===")
for (pr in pairs) {
  a <- pr[1]; b <- pr[2]
  A <- lSets[[a]]; B <- lSets[[b]]
  lines <- c(
    lines,
    sprintf("  %s vs %s: (%s)", a, b, pair_desc(a, b)),
    sprintf("    total %-6s : %5d", a, length(A)),
    sprintf("    total %-6s : %5d", b, length(B)),
    sprintf("    %s only      : %5d", a, length(setdiff(A, B))),
    sprintf("    %s only      : %5d", b, length(setdiff(B, A))),
    sprintf("    shared       : %5d", length(intersect(A, B))),
    ""
  )
}
A <- lSets[["T.DMSO"]]; B <- lSets[["T.ROSI"]]; C <- lSets[["T.SCO"]]
onlyA <- setdiff(A, union(B, C))
onlyB <- setdiff(B, union(A, C))
onlyC <- setdiff(C, union(A, B))
AB <- setdiff(intersect(A, B), C)
AC <- setdiff(intersect(A, C), B)
BC <- setdiff(intersect(B, C), A)
ABC <- Reduce(intersect, list(A, B, C))
lines <- c(lines,
           "=== 3-way Venn regions (adj.P.Val < 0.05) ===",
           sprintf("  only   T.DMSO : %5d", length(onlyA)),
           sprintf("  only   T.ROSI : %5d", length(onlyB)),
           sprintf("  only   T.SCO  : %5d", length(onlyC)),
           sprintf("  T.DMSO???T.ROSI : %5d", length(AB)),
           sprintf("  T.DMSO???T.SCO  : %5d", length(AC)),
           sprintf("  T.ROSI???T.SCO  : %5d", length(BC)),
           sprintf("  shared all 3  : %5d", length(ABC)),
           ""
)
total_rev <- length(unique(c(reversed_SCO, reversed_ROSI)))
lines <- c(lines,
           sprintf("=== Reversal genes (TNF p<%.2g; treat p<%.2g) ===", rev_tnf_cut, rev_treat_cut),
           sprintf("  Total reversed  : %5d", total_rev),
           sprintf("    SCO-only      : %5d", length(sco_only)),
           sprintf("    ROSI-only     : %5d", length(rosi_only)),
           sprintf("    Shared        : %5d", length(shared_rev)),
           ""
)
get_overlap_directionality <- function(sig_tnf, sig_treat) {
  overlap <- intersect(rownames(sig_tnf), rownames(sig_treat))
  logFC_tnf <- sig_tnf[overlap, "logFC"]
  logFC_trt <- sig_treat[overlap, "logFC"]
  opposite <- sum(sign(logFC_tnf) != sign(logFC_trt))
  same <- length(overlap) - opposite
  list(
    overlap = length(overlap),
    opposite = opposite,
    same = same,
    pct_opposite = 100 * opposite / length(overlap),
    pct_same = 100 * same / length(overlap)
  )
}
sco_stats <- get_overlap_directionality(sig_TNF, sig_SCO)
rosi_stats <- get_overlap_directionality(sig_TNF, sig_ROSI)
lines <- c(lines,
           "=== Overlap and Directionality Between TNF and SCO (adj.P.Val < 0.05) ===",
           sprintf("Total genes regulated by TNF: %d", nrow(sig_TNF)),
           sprintf("TNF-regulated genes also regulated by SCO: %d (%.1f%%)", sco_stats$overlap, 100 * sco_stats$overlap / nrow(sig_TNF)),
           sprintf("  Opposite direction vs TNF: %d (%.1f%%)", sco_stats$opposite, sco_stats$pct_opposite),
           sprintf("  Same direction as TNF:     %d (%.1f%%)", sco_stats$same, sco_stats$pct_same),
           "",
           "=== Overlap and Directionality Between TNF and ROSI (adj.P.Val < 0.05) ===",
           sprintf("Total genes regulated by TNF: %d", nrow(sig_TNF)),
           sprintf("TNF-regulated genes also regulated by ROSI: %d (%.1f%%)", rosi_stats$overlap, 100 * rosi_stats$overlap / nrow(sig_TNF)),
           sprintf("  Opposite direction vs TNF: %d (%.1f%%)", rosi_stats$opposite, rosi_stats$pct_opposite),
           sprintf("  Same direction as TNF:     %d (%.1f%%)", rosi_stats$same, rosi_stats$pct_same),
           ""
)
cat(paste(lines, collapse = "\n"), "\n")
# ---- END PIPELINE ----
