args = commandArgs(trailingOnly=TRUE)

rmskout = args[1]
genomesize = as.numeric(args[2]) / 1000
species = sub(pattern = ".fasta.out", x = rmskout, replacement = "")

library(data.table)

DATA = fread(file = rmskout, header = F, stringsAsFactors = F, skip = 3, fill = T)

DATA = DATA[,c(5,6,7,9,10,11,2,15)]
names(DATA) = c("Scaffold", "Begin", "End", "Strand", "Element", "Family", "Divergence", "ID")

# select only hits with divergence strictly lesser than 3%
DATA = DATA[DATA$Divergence < 3,]

# select nonLTR
boo = grepl(pattern = "LINE|SINE", x = DATA$Family)
DATA = DATA[boo,]

# calculate density of insertion as number of insertions divided by genome size in Gb
# use unique IDs
total_number = length(unique(DATA$ID))
# genome size is originally in Mb, here convert it in Gb
density_insertion = total_number / genomesize

summary = data.frame(Species = species, Density = density_insertion)

write.table(x = summary, file = "nonLTR_DI_3_perc.out", quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE, append = TRUE)