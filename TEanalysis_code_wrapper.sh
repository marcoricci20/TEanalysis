## Ricci et al 2021 - Comparative analysis of bats and rodents' genomes suggests a relation between
## non-LTR retrotransposons, cancer incidence, and ageing


# -------------------------------------------------

## Code for making the landscape plots: Figure 3 and Figure 4

# RepeatMasker out must be in the format: hetGla.fasta.out, ratNor.fasta.out etc etc

# create file for genome sizes (necessary to calculate the percentage of repeats at each divergence level)
printf 'cavPor\t2723.20284\nmusMus\t2728.222451\nratNor\t2647.915728\nhetGla\t3041.864219\nmyoLuc\t2034.5753\nmyoMyo\t2002.797769\nrhiFer\t2075.104659\npteVam\t1996.07641\nrouAeg\t1893.602072\nmolMol\t2319.008189\n' > genome_size.txt

# load R libraries
module load R_packages/3.6.0

# the script needs to find hetGla.fasta.out, ratNor.fasta.out etc etc files in the folder where it is run
# it also needs genome_size.txt in the same folder
Rscript --vanilla TE_lanscape_plots.R 

# -------------------------------------------------

## Code for making the barplot of Figure 5

# this script calculates the number of insertions for non-LTR elements with a divergence < 3%
# it saves the density of insertions in the file nonLTR_DI_3_perc.out
for OUT in $( ls *.fasta.out )
do
 SPECIES=${OUT%.*.*}
 GS=`grep $SPECIES genome_size.txt | cut -f2`
 Rscript --vanilla nonLTR_DI_3_perc.R $OUT $GS
done

# add lifespans and category for the barplot
# order by lifespan
# color by category
printf 'Short lifespan prone to cancer\t12\nLong lifespan resistant to cancer\t31\nShort lifespan resistant to cancer\t5.6\nShort lifespan prone to cancer\t4\nLong lifespan resistant to cancer\t34\nLong lifespan resistant to cancer\t37.1\nLong lifespan resistant to cancer\t20.9\nShort lifespan prone to cancer\t3.8\nLong lifespan resistant to cancer\t30.5\nLong lifespan resistant to cancer\t22.9\n' > lifespan_categories.txt
paste -d '\t' nonLTR_DI_3_perc.out lifespan_categories.txt > nonLTR_DI_3_perc_barplot.out

# this script creates the barplot of Figure 5
# it requires nonLTR_DI_3_perc_barplot.out
# it outputs nonLTR_DI_3_perc_barplot.pdf
Rscript --vanilla nonLTR_DI_3_perc_barplot.R 

# -------------------------------------------------

## Code to run Wilcoxon paired ranked test on lifespan and density of insertion

# The script needs dataset_wilcoxon.txt file (found in the supplementary materials and Github)
Rscript --vanilla Wilcoxon_test.R dataset_wilcoxon.txt

# -------------------------------------------------

