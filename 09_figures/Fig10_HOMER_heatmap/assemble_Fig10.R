#!/usr/bin/env Rscript
# Combine the NS + Stable integrated HOMER heatmap PNGs side by side into Fig 10.
# Reads the two PNGs from Fig10_integrate_heatmap.R (override dir with HEATMAP_DIR).

library(magick)

heatmap.dir <- Sys.getenv(
  "HEATMAP_DIR",
  "/QRISdata/Q8448/Mouse_disease_data/DAR/Integrated_HOMER_Heatmaps"
)
ns.png     <- file.path(heatmap.dir, "Integrated_NS_heatmap.png")
stable.png <- file.path(heatmap.dir, "Integrated_Stable_heatmap.png")
out.dir    <- heatmap.dir

# Trim whitespace, add uniform border
img.ns     <- image_read(ns.png)     |> image_trim() |> image_border("white", "40x40")
img.stable <- image_read(stable.png) |> image_trim() |> image_border("white", "40x40")

# Match heights before placing side by side
h <- max(image_info(img.ns)$height, image_info(img.stable)$height)
img.ns     <- image_resize(img.ns,     paste0("x", h))
img.stable <- image_resize(img.stable, paste0("x", h))

# Separator between panels
w.sep <- 80
separator <- image_blank(w.sep, h, color = "white")

# Place: NS on left, Stable on right
combined <- image_append(c(img.ns, separator, img.stable), stack = FALSE)

image_write(combined, path = file.path(out.dir, "Fig10_integrated_heatmap.png"), format = "png")
image_write(combined, path = file.path(out.dir, "Fig10_integrated_heatmap.pdf"), format = "pdf")

message("Saved: Fig10_integrated_heatmap.png / .pdf")
message("Dimensions: ", image_info(combined)$width, " x ", image_info(combined)$height, " px")
