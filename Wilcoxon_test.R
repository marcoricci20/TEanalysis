args = commandArgs(trailingOnly=TRUE)

filename = args[1]

data = read.table(file = filename, header = TRUE, stringsAsFactors = FALSE, sep = "\t")

# Run the test
res = wilcox.test(x = data$DI_LONG, y = data$DI_SHORT, paired = TRUE, alternative = "less")

sink(file = "Wilcoxon_test.out")
cat(res)
sink()