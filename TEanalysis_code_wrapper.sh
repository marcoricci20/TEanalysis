## Ricci et al 2021 - Comparative analysis of bats and rodents' genomes suggests a relation between
## non-LTR retrotransposons, cancer incidence, and ageing

# -------------------------------------------------

## Code used for running RepeatModeler2 on Heterocephalus glaber genome assembly

module load bioinfo-tools
module load RepeatModeler/2.0.3

BuildDatabase -name hetGlaDB hetGla.fasta
RepeatModeler -database hetGlaDB -pa 20 -LTRStruct
# hetGla_rm1.0.lib is the library with the new consensus sequences for the species

# -------------------------------------------------

## Code for building the hetGla library

module load bioinfo-tools
module load RepeatMasker/4.1.0

queryRepeatDatabase.pl -species rodents > rodents.lib
grep 'L2|MIR' | cut -c2- rodents.lib > list_ancient_rodent_TEs.txt
perl extractFromFasta.pl rodents.lib list list_ancient_rodent_TEs.txt > list_ancient_rodent_TEs.fasta
cat hetGla_rm1.0.lib list_ancient_rodent_TEs.fasta > hetGla_rm1.1.lib

# -------------------------------------------------

## The repeat libraries for the other rodents were already established and directly specified in the RepeatMasker command
## The repeat library for bats was retrieved from Jebb et al 2020. Here the bat library is called bats.lib (includes Repbase)

# -------------------------------------------------

## Code used for running RepeatMasker on all the genome assemblies

module load bioinfo-tools
module load RepeatMasker/4.1.0

# cavPor
RepeatMasker -pa 20 -a -xsmall -gccalc -excln -species cavia porcellus cavPor.fasta

# hetGla
RepeatMasker -pa 20 -a -xsmall -gccalc -excln -lib hetGla_rm1.1.lib hetGla.fasta

# ratNor
RepeatMasker -pa 20 -a -xsmall -gccalc -excln -species rattus norvegicus ratNor.fasta

# musMus
RepeatMasker -pa 20 -a -xsmall -gccalc -excln -species mus musculus musMus.fasta

# molMol
RepeatMasker -pa 20 -a -xsmall -gccalc -excln -lib bats.lib molMol.fasta

# myoLuc
RepeatMasker -pa 20 -a -xsmall -gccalc -excln -lib bats.lib myoLuc.fasta

# myoMyo
RepeatMasker -pa 20 -a -xsmall -gccalc -excln -lib bats.lib myoMyo.fasta

# pteVam
RepeatMasker -pa 20 -a -xsmall -gccalc -excln -lib bats.lib pteVam.fasta

# rhiFer
RepeatMasker -pa 20 -a -xsmall -gccalc -excln -lib bats.lib rhiFer.fasta

# rouAeg
RepeatMasker -pa 20 -a -xsmall -gccalc -excln -lib bats.lib rouAeg.fasta

# the output files have names like: hetGla.fasta.out, molMol.fasta.out, etc...

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

