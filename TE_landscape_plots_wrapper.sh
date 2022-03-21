# RepeatMasker out must be in the format: hetGla.fasta.out, ratNor.fasta.out etc etc

printf 'cavPor\t2723.20284\nmusMus\t2728.222451\nratNor\t2647.915728\nhetGla\t3041.864219\nmyoLuc\t2034.5753\nmyoMyo\t2002.797769\nrhiFer\t2075.104659\npteVam\t1996.07641\nrouAeg\t1893.602072\nmolMol\t2319.008189\n' > genome_size.txt

for OUT in $( ls *.fasta.out )
do
 SPECIES=${OUT%.*.*}
 GS=`grep $SPECIES genome_size.txt | cut -f2`
 PREFIX=${SPECIES}_landscape
 Rscript --vanilla TE_lanscape_plots.R $OUT $GS $PREFIX
done 
