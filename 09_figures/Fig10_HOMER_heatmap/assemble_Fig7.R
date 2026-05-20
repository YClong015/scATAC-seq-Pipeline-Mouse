suppressPackageStartupMessages(library(magick))

ns_png     <- "/Users/zhengyanchen/Desktop/project_outline_figures/Integrated_NS_heatmap.png"
stable_png <- "/Users/zhengyanchen/Desktop/project_outline_figures/Integrated_Stable_heatmap.png"
out_dir    <- "/Users/zhengyanchen/my_claude/scripts/CMpaper_heatmap"

# Trim whitespace, add uniform border
img_ns     <- image_read(ns_png)     |> image_trim() |> image_border("white", "40x40")
img_stable <- image_read(stable_png) |> image_trim() |> image_border("white", "40x40")

# Match heights before placing side by side
h <- max(image_info(img_ns)$height, image_info(img_stable)$height)
img_ns     <- image_resize(img_ns,     paste0("x", h))
img_stable <- image_resize(img_stable, paste0("x", h))

# Separator between panels
w_sep <- 80
separator <- image_blank(w_sep, h, color = "white")

# Place: NS on left, Stable on right
combined <- image_append(c(img_ns, separator, img_stable), stack = FALSE)

image_write(combined,
            path   = file.path(out_dir, "Fig7_integrated_heatmap.png"),
            format = "png")
image_write(combined,
            path   = file.path(out_dir, "Fig7_integrated_heatmap.pdf"),
            format = "pdf")

message("Saved: Fig7_integrated_heatmap.png / .pdf")
message("Dimensions: ",
        image_info(combined)$width, " x ", image_info(combined)$height, " px")
