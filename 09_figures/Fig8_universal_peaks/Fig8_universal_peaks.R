#!/usr/bin/env Rscript
# Fig 8 - Universal peak set composition 
# Input: upset_input.csv from merge_universal_bedtools.sh

library(ggplot2)
library(dplyr)
library(ComplexUpset)   # install if needed: install.packages("ComplexUpset")
library(patchwork)

## Paths
upset.csv <- "/QRISdata/Q8448/Mouse_disease_data/universal_peaks_v2/upset_input.csv"
out.dir   <- "/QRISdata/Q8448/Mouse_disease_data/QC_figures"
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

## Load data
df <- read.csv(upset.csv, stringsAsFactors = FALSE)
cat("Total universal peaks:", nrow(df), "\n")

tissues <- c("Kidney", "Lung", "Aorta", "Tcells")

# Convert 0/1 to logical (required by ComplexUpset)
for (t in tissues) df[[t]] <- as.logical(df[[t]])

## Per-tissue peak counts (for Fig 8 bar chart)
tissue.counts <- data.frame(
  Tissue = c("Kidney", "Lung", "Aorta", "T cells", "Universal"),
  Peaks  = c(
    sum(df$Kidney),
    sum(df$Lung),
    sum(df$Aorta),
    sum(df$Tcells),
    nrow(df)
  )
)
tissue.counts$Tissue <- factor(
  tissue.counts$Tissue,
  levels = c("Kidney", "Lung", "Aorta", "T cells", "Universal")
)
tissue.counts$Type <- ifelse(tissue.counts$Tissue == "Universal",
                             "Universal", "Tissue")

cat("\n=== Per-tissue peak counts ===\n")
print(tissue.counts)

# Fig 8 - Bar chart
p5 <- ggplot(tissue.counts,
             aes(x = Tissue, y = Peaks / 1000, fill = Type)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = formatC(Peaks, format = "d", big.mark = ",")),
            vjust = -0.4, size = 3.5) +
  geom_vline(xintercept = 4.5, linetype = "dashed",
             colour = "grey50", linewidth = 0.6) +
  scale_fill_manual(values = c("Tissue" = "#51247A", "Universal" = "#962A8B"),
                    guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "A  Universal peak set composition",
    x     = NULL,
    y     = "Number of peaks (x1,000)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title  = element_text(face = "bold", size = 13),
    axis.text.x = element_text(size = 11)
  )

ggsave(file.path(out.dir, "Fig8_universal_peaks.pdf"), p5,
       width = 7, height = 5)
ggsave(file.path(out.dir, "Fig8_universal_peaks.png"), p5,
       width = 7, height = 5, dpi = 300)
message("Saved: Fig8_universal_peaks")

## Fig 8 (UpSet) - UpSet plot (ComplexUpset)
# Tissue membership summary for console
cat("\n=== Tissue membership breakdown ===\n")
cat("Tissue-specific peaks:\n")
for (t in tissues) {
  n <- sum(df[[t]] & df$n_tissues == 1)
  cat(sprintf("  %-10s: %d\n", t, n))
}
cat("Shared 2 tissues:", sum(df$n_tissues == 2), "\n")
cat("Shared 3 tissues:", sum(df$n_tissues == 3), "\n")
cat("Shared 4 tissues:", sum(df$n_tissues == 4), "\n")

# Display labels
tissue.labels <- c(
  Kidney = "Kidney", Lung = "Lung",
  Aorta = "Aorta", Tcells = "T cells"
)

p6 <- upset(
  df,
  intersect    = tissues,
  labeller     = labeller(intersection = tissue.labels),
  width_ratio  = 0.2,
  min_size     = 100,           # hide intersections with <100 peaks
  sort_intersections_by = "cardinality",
  base_annotations = list(
    "Intersection size" = intersection_size(
      counts      = TRUE,
      text        = list(size = 3),
      bar_number_threshold = 0,
      mapping     = aes(fill = "bar")
    ) +
      scale_fill_manual(values = c(bar = "#51247A"), guide = "none") +
      labs(y = "Intersection size") +
      theme(axis.text.y = element_text(size = 9))
  ),
  set_sizes = (
    upset_set_size(
      geom = geom_bar(fill = "#51247A", width = 0.6)
    ) +
      labs(x = "Set size") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
  ),
  themes = upset_modify_themes(list(
    "intersections_matrix" = theme(text = element_text(size = 11)),
    "overall_sizes"        = theme(text = element_text(size = 10))
  ))
) +
  labs(title = "B  Tissue peak overlap in universal peak set") +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0))

ggsave(file.path(out.dir, "Fig8_upset.pdf"), p6,
       width = 12, height = 6)
ggsave(file.path(out.dir, "Fig8_upset.png"), p6,
       width = 12, height = 6, dpi = 300)
message("Saved: Fig8_upset")

## Print final numbers for thesis text
cat("\n=== Numbers for thesis text ===\n")
cat(sprintf("Total universal peaks:   %s\n",
            formatC(nrow(df), format="d", big.mark=",")))
cat(sprintf("Kidney-only:             %s\n",
            formatC(sum(df$Kidney & df$n_tissues==1), format="d", big.mark=",")))
cat(sprintf("Lung-only:               %s\n",
            formatC(sum(df$Lung   & df$n_tissues==1), format="d", big.mark=",")))
cat(sprintf("Aorta-only:              %s\n",
            formatC(sum(df$Aorta  & df$n_tissues==1), format="d", big.mark=",")))
cat(sprintf("Tcells-only:             %s\n",
            formatC(sum(df$Tcells & df$n_tissues==1), format="d", big.mark=",")))
cat(sprintf("Shared all 4 tissues:    %s\n",
            formatC(sum(df$n_tissues==4), format="d", big.mark=",")))

message("\nDone. Outputs -> ", out.dir)
