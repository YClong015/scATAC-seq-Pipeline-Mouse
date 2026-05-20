suppressPackageStartupMessages({
  library(magick)
})

pdf_path <- "/Users/zhengyanchen/Desktop/HOMER_Plots/DAR_closing_vs_opening/Cell_type_identity.pdf"
out_dir  <- "/Users/zhengyanchen/my_claude/scripts/CMpaper_heatmap"

# Selected pages (1-indexed):
# Page 1:  Tcell_Tcell   — TCF3/LEF1/TCF7, canonical T cell identity
# Page 2:  Kidney_DCT    — HNF4a/HNF1b/PPARa, distal tubule
# Page 4:  Kidney_PT     — HNF1b/PPARa/HNF4a, proximal tubule (most disease-relevant)
# Page 6:  Lung_AT2      — Smad2/HIF-1a, alveolar epithelial
# Page 10: Lung_Fib      — KLF/ETS, fibrosis-relevant fibroblast
# Page 17: Aorta_SMC     — WT1/AP-2alpha/Klf4, primary vascular disease cell
#
# Excluded: Eosinophils/NK/B (weak signal), Aorta_Mac (Inconclusive),
#           PC/TAL/Pericytes (redundant with selected)

selected <- list(
  list(page = 4,  label = "A  Kidney — PT"),
  list(page = 2,  label = "B  Kidney — DCT"),
  list(page = 6,  label = "C  Lung — AT2"),
  list(page = 10, label = "D  Lung — Fib"),
  list(page = 17, label = "E  Aorta — SMC"),
  list(page = 1,  label = "F  T cell")
)

message("Reading pages from PDF...")
imgs <- lapply(selected, function(s) {
  img <- image_read_pdf(pdf_path, pages = s$page, density = 200)
  # Trim whitespace, then add a small border for spacing
  img <- image_trim(img)
  img <- image_border(img, "white", "10x10")
  img
})

# Layout: 2 columns × 3 rows
# Row 1: Kidney_PT + Kidney_DCT
# Row 2: Lung_AT2  + Lung_Fib
# Row 3: Aorta_SMC + Tcell

# Make each image the same width before appending
max_w <- max(sapply(imgs, function(i) image_info(i)$width))
imgs_resized <- lapply(imgs, function(i) {
  info <- image_info(i)
  if (info$width != max_w) {
    image_resize(i, paste0(max_w, "x"))
  } else {
    i
  }
})

row1 <- image_append(c(imgs_resized[[1]], imgs_resized[[2]]), stack = FALSE)
row2 <- image_append(c(imgs_resized[[3]], imgs_resized[[4]]), stack = FALSE)
row3 <- image_append(c(imgs_resized[[5]], imgs_resized[[6]]), stack = FALSE)

combined <- image_append(c(row1, row2, row3), stack = TRUE)

out_png <- file.path(out_dir, "Fig6_cell_identity_selected.png")
out_pdf <- file.path(out_dir, "Fig6_cell_identity_selected.pdf")

image_write(combined, path = out_png, format = "png", density = "300")
image_write(combined, path = out_pdf, format = "pdf")

message("Saved: ", out_png)
message("Saved: ", out_pdf)
message("Dimensions: ", image_info(combined)$width, " x ", image_info(combined)$height, " px")
