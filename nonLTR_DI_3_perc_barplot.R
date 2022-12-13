## Code for the barplots with lifespan and cancer incidence and density of insertion

filename = "nonLTR_DI_3_perc_barplot.out"

library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# read input
DATA = fread(file = filename, stringsAsFactors = T, fill = T, skip = 1, sep ="\t")

# give names to columns
names(DATA) = c("Species", "Density", "Category", "Lifespan")

# categories of lifespan and cancer incidence
Categories = c("Short lifespan prone to cancer","Short lifespan resistant to cancer","Long lifespan resistant to cancer")
colors = c("#FF0000","#ff9a00ff","#00be00ff")

# modify levels for plotting in the right order and colors
DATA$Species = factor(x = DATA$Species, levels = DATA$Species[order(DATA$Lifespan)])
DATA$Category = factor(x = DATA$Category, levels = Categories)

# order by increasing lifespan
DATA = DATA[order(DATA$Lifespan, decreasing = FALSE),]

# make plot
DI_plot = ggplot(data = DATA, aes(x = Species, y = Density, fill = Category )) + geom_bar(stat = "identity", width = 0.40, ) + scale_fill_manual(name = "Category", values = colors) + theme_bw() + theme(panel.grid.major.x = element_blank()) # remove the vertical grid lines

# save plot as PDF
ggsave(filename = 'nonLTR_DI_3_perc_barplot.pdf', plot = DI_plot, device = "pdf", width = 65, height = 27, units = "cm", scale = .5)
