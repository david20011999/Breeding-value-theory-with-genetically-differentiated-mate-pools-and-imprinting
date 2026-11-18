# Requires map and diff_50_abs from the cattle analysis in the R session.
# From the project root: source("Figures/RealData_Figures.R")
figure_dir <- normalizePath("Figures")

#Figure 8

map$diff <- diff_50_abs
map$V3 <- 1:nrow(map)
map$is_highlight <- NA

map$is_highlight[map$V4>747465 & map$V4<8824176 & map$V1 == 2] <- "yes"
map$is_highlight[map$V4>114416811  & map$V4<114439054 & map$V1 == 4] <- "yes"
map$is_highlight[map$V4>70680628 & map$V4<71749458 & map$V1 == 6] <- "yes"
map$is_highlight[map$V4>74723634 & map$V4<74925889 & map$V1 == 17] <- "yes"

aux <- as.data.frame(tapply(map$V4, list(map$V1, map$is_highlight=="yes"), quantile,0.50))
colnames(aux) <- "pos"
aux$chr <- row.names(aux)
aux <- aux[!is.na(aux$pos),]
aux$label <- c("MSTN", "KIT", NA,NA)

map$label <- NA
map$label[map$V4==4421299 & map$V1 == 2] <- "MSTN"
map$label[map$V4==114437170 & map$V1 == 4] <- "KIT"


source(file.path(figure_dir, "figure_style.R"))
p1 <- plot_cattle_differences(list(map=map))
if (interactive()) print(p1)
save_paper_plot(p1, "diff_50", paper_width_mm, 62)
