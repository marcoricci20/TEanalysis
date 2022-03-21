# usage: Rscript --vanilla TE_landscape_plots.R RMSK.out genomesize prefix4plot

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 3) {
  stop("Three arguments must be supplied", call.=FALSE)
} else {
  
  filename = args[1]
  genomesize = as.numeric(args[2])
  plotname = args[3]

}

# libraries required
library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# read RepeatMasker output file
DATA = fread(file = filename, header = F, stringsAsFactors = F, skip = 3, fill = T)

# subset relevant columns
DATA = DATA[,c(5,6,7,9,10,11,2,15)]
names(DATA) = c("Scaffold", "Begin", "End", "Strand", "Element", "Family", "Divergence", "ID")

# add length of the hits
DATA$Length = DATA$End - DATA$Begin + 1

# delete elements with divergence greater than 100 (there could be artefacts sometimes)
DATA = DATA[DATA$Divergence < 100,]

# re-order the column
DATA = DATA[,c(1:6,9,7,8)]

# replace "C" with "-" in the Strand column
DATA$Strand = sub(pattern = "C", replacement = "-", x = DATA$Strand)

# Discard elements we are not interested into
DISCARD = c("Low_complexity", "Satellite", "Simple_repeat", "tRNA", "rRNA", "scRNA", "srpRNA", "snRNA", "Unspecified", "RNA")
DATA = DATA[!DATA$Family %in% DISCARD,]
boo = grepl(pattern = "DNA\\?|LTR\\?|RC\\?|SINE\\?|Retrogene|Retroposon|RNA|rRNA|scRNA|Segmenta|Simple_repeat|snRNA|tRNA|Satellite", x = DATA$Family)
DATA = DATA[!boo,]

# round the divergence values
DATA$RoundDiv = floor(DATA$Divergence)

# group elements in major classes
elements = sapply(strsplit(x = DATA$Family,split = "/"), "[", 1)
DATA$newFamily = elements

# create a factor with the name of the subfamily/element and the divergence associated to it
# so we can get the number of bps associated to that particular subfamily at that particular divergence
DATA$Factor = paste(DATA$newFamily, DATA$RoundDiv, sep = "$")

# general landscape - bps occupied
DATA_bps = aggregate(Length ~ Factor, DATA, sum)
DATA_bps$Element = sapply(strsplit(DATA_bps$Factor, "\\$"), "[[", 1)
DATA_bps$Divergence = sapply(strsplit(DATA_bps$Factor, "\\$"), "[[", 2)

# conversion in megabases
DATA_bps$Mb = DATA_bps$Length / 1000000
# conversion in percentage of genome aggiungi le megabasi del genoma di volta in volta
DATA_bps$Percentage = (DATA_bps$Mb / genomesize) * 100

elements = c("DNA","RC","SINE","LINE","LTR","Unknown")
colors = data.frame(Element = elements, Color = c("#0070C0","#00FFCC","#61E101", "#048A1E", "#FF0000", "#7030A0"))

colors = colors[colors$Element %in% unique(DATA_bps$Element),]

DATA_bps$Element = factor(x = DATA_bps$Element, levels = elements)

# plot
DATA_plot = ggplot(data = DATA_bps, aes(x = as.integer(Divergence), y = Percentage, fill = Element)) + geom_bar(stat = "identity") + theme_bw() + theme(panel.border = element_rect(colour = 'darkgrey', fill = NA)) + xlab("Divergence") + scale_fill_manual(name = "Categories", values = colors) + ylab("Percentage") + xlim(NA,50)+ ylim(NA,2.7)

ggsave(filename = paste0(plotname, ".pdf"), plot = DATA_plot, device = "pdf", width = 50, height = 27, units = "cm", scale = .5)
