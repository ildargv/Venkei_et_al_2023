#1.trim the adapter
module load fastx_toolkit/0.0.14
for i in $(ls *.fastq); do bsub -q short -W 240 -n 4  -e $i.err "fastx_clipper -a TGGAATTCTCGGGTGCCAAGG -l 45 -c -v -i $i > $i.trimmed"; done

#2.filter reads at p20q100 (Phred score 20+ for all nucleotides)
module load fastx_toolkit/0.0.14
for i in $(find . -maxdepth 1 -name "*.trimmed"); do /pi/phillip.zamore-umw/home/ildar.gainetdinov-umw/common/pipelines/ms_short_8G "fastq_quality_filter -q 20 -p 100 -i $i -o $i.q20p100" ; done 

#3.UMI:remove duplicates and anything <18nt
for i in $(find . -maxdepth 1 -name "*.trimmed.q20p100"); do /pi/phillip.zamore-umw/home/ildar.gainetdinov-umw/common/pipelines/ms_short_64G "awk '(NR%4==1){name=\$1}(NR%4==2){total++;if((substr(\$1,length(\$1)-11,3)==\"GTC\")&&(substr(\$1,length(\$1)-5,3)==\"TAG\")&&(((substr(\$1,4,3)==\"CGA\")&&(substr(\$1,10,3)==\"TAC\"))||((substr(\$1,4,3)==\"ATC\")&&(substr(\$1,10,3)==\"AGT\")))){umis++;if(a[\$1]!=1){nondup++;a[\$1]=1;if(length(\$1)>47){longer18++; print name; print substr (\$1,16,length (\$1)-30);getline; print; getline; print substr (\$1,16,length (\$1)-30)}}}}END{print FILENAME\"\\t\"total\"\\t\"umis\"\\t\"nondup\"\\t\"longer18 > FILENAME\".dup\"}' $i > $i.deUMI.dedup.fq" ; done 

#4.remove rRNA from genbank M21017.1 with 1 mismatch allowed
module load bowtie/1.3.1
for i in $(find . -name "*.fq"); do /pi/phillip.zamore-umw/home/ildar.gainetdinov-umw/common/pipelines/ms_short_8G "bowtie --un $i.x_filter_1mm -k 1 -v 1 /pi/phillip.zamore-umw/home/ildar.gainetdinov-umw/common/dm6_set/dm6.filter $i > /dev/null" ; done 

#5.remove spikeins
module load bowtie/1.3.1
for i in $(find . -name "*.x_filter_1mm"); do /pi/phillip.zamore-umw/home/ildar.gainetdinov-umw/common/pipelines/ms_short_8G "bowtie --norc --un $i.x_spikein.fq -v 0 /pi/phillip.zamore-umw/home/ildar.gainetdinov-umw/common/dm6_set/final.set $i | awk '(substr(\$3,index(\$3,\"-\")+1)==length(\$5)){spike[\$3]++}END{for (i in spike){print i\"\\t\"spike[i]}}' > $i.spikein" ; done 

module load gawk/4.1.4
gawk '{a[$1][FILENAME]=$2}END{printf "name";for (i in a["114-26"]){printf "\t"i};printf("\n");for (r in a){printf r; for (f in a["114-26"]){if (a[r][f]==""){printf "\t0"}else{printf "\t"a[r][f]}};printf ("\n")}}' *.spikein > spikein.xls

#6.send to Tailor
for i in `find . -name "*.x_filter_1mm"`; do /pi/phillip.zamore-umw/home/ildar.gainetdinov-umw/common/pipelines/ms "awk 'BEGIN{FS=OFS=\" \"}(FNR%4==1){print $1;getline; print;getline;print \"+\";getline;print}' $i > $i.fq " ; done

for i in `find . -name "*.x_filter_1mm.fq"`; do /pi/phillip.zamore-umw/home/ildar.gainetdinov-umw/common/pipelines/ms_tailor "/pi/phillip.zamore-umw/home/ildar.gainetdinov-umw/common/pipelines/Tailor/run_tailing_pipeline.sh -q 20 -T 10 -t /pi/phillip.zamore-umw/home/ildar.gainetdinov-umw/common/pipelines/Tailor/annotation/dm6.genomic_features -i $i -c 8 -g /pi/phillip.zamore-umw/home/ildar.gainetdinov-umw/common/pipelines/Tailor/indexes/dm6.pipipes.fa -o $i.Tailor" ; done

#7.convert bed2 files from Tailor output to rpm files
for i in $(find . -name "*.p20.bed2"); do  /pi/phillip.zamore-umw/home/ildar.gainetdinov-umw/common/pipelines/ms_short_8G "awk '((\$9==0)&&(dealt[\$7]!=1)){nontailed[substr (\$7, 1, length (\$7)-\$9)]+=\$4; all[substr (\$7, 1, length (\$7)-\$9)]=1;coordinates[substr (\$7, 1, length (\$7)-\$9)]=\$1\":\"\$2\"-\"\$3\"\\t\"\$6; alluniq+=\$4; dealt[\$7]=1}((\$9>0)&&(dealt[\$7]!=1)){tailed[substr (\$7, 1, length (\$7)-\$9)]+=\$4;all[substr (\$7, 1, length (\$7)-\$9)]=1 ;coordinates[substr (\$7, 1, length (\$7)-\$9)]=\$1\":\"\$2\"-\"\$3\"\\t\"\$6; alluniq+=\$4; dealt[\$7]=1}END{print \"coordinates\\tstrand\\tsequenceASis\\ttotalRPM\\ttailedRPM\"; for (r in all){print coordinates[r]\"\\t\"r\"\\t\"(tailed[r]+nontailed[r])*1000000/alluniq\"\\t\"tailed[r]*1000000/alluniq}}' $i |sort -k5 -g -r > $i.rpm"; done

#8.align rpm files onto dm6 with 0 mismatches
module load bowtie/1.3.1
for i in $(find . -name "*.bed2.rpm"); do /pi/phillip.zamore-umw/home/ildar.gainetdinov-umw/common/pipelines/ms "awk '(length(\$3)>=20){print \">\"\$3\";\"\$4;print \$3}' $i > $i.20.fa; bowtie -f -v 0 -a /pi/phillip.zamore-umw/home/ildar.gainetdinov-umw/common/dm6_set/dm6.pipipes $i.20.fa | awk '{print \$3\"\\t\"\$4\"\\t\"(\$4+length(\$5))\"\\t\"\$1\"\\t\"\$7+1\"\\t\"\$2}' > $i.20.bed2; awk 'BEGIN{FS=OFS=\"\\t\"}(FNR==NR){a[\$4]++}(FNR<NR){\$5=a[\$4];split(\$4,b,\";\");\$4=b[2];\$7=b[1]; print}' $i.20.bed2 $i.20.bed2 > $i.20.realbed2; rm $i.20.bed2 $i.20.fa"; done

#9.get stranded length distribution of flam, all reads apportioned
module load bedtools/2.30.0
for i in $(ls *.20.realbed2); do /pi/phillip.zamore-umw/home/ildar.gainetdinov-umw/common/pipelines/ms_short_8G "bedtools intersect -wo -a $i -b flam.bed | awk 'BEGIN{FS=OFS=\"\\t\"}{if(\$6==\$13){ aplus[(\$3-\$2)]+=\$4/\$5}else{aminus[(\$3-\$2)]+=\$4/\$5}}END{printf (\"plus\");for (i=20;i<=30;i++){if(aplus[i]==\"\"){printf (\"\\t0\")}else{printf (\"\\t\"aplus[i])}};printf(\"\\n\"); printf (\"minus\");for (i=20;i<=30;i++){if(aminus[i]==\"\"){printf (\"\\t0\")}else{printf (\"\\t-\"aminus[i])}};printf(\"\\n\")}' >  $i.flam_lendist"; done

#10.get stranded length distribution of SuSte and Ste, both uniquely and multiply mapping
module load bedtools/2.30.0
for i in $(ls *.20.realbed2); do /pi/phillip.zamore-umw/home/ildar.gainetdinov-umw/common/pipelines/ms_short_8G "bedtools intersect -wo -a $i -b nonredste.bed > $i.nonredste.tmp; bedtools intersect -wo -a $i -b 24suste_w_protop_subtracted.bed > $i.24suste_w_protop_subtracted.tmp;
awk 'BEGIN{FS=OFS=\"\\t\"}(FNR==NR){found[\$7]=1}(FNR<NR)&&(found[\$7]==\"\"){if(\$6==\$13){ aplus[(\$3-\$2)]+=\$4/\$5}else{aminus[(\$3-\$2)]+=\$4/\$5}}END{printf (\"plus\");for (i=20;i<=30;i++){if(aplus[i]==\"\"){printf (\"\\t0\")}else{printf (\"\\t\"aplus[i])}};printf(\"\\n\"); printf (\"minus\");for (i=20;i<=30;i++){if(aminus[i]==\"\"){printf (\"\\t0\")}else{printf (\"\\t-\"aminus[i])}};printf(\"\\n\")}' $i.24suste_w_protop_subtracted.tmp $i.nonredste.tmp >  $i.nonredste_uniq_lendist;
awk 'BEGIN{FS=OFS=\"\\t\"}(FNR==NR){found[\$7]=1}(FNR<NR)&&(found[\$7]==1){if(\$6==\$13){ aplus[(\$3-\$2)]+=\$4/\$5}else{aminus[(\$3-\$2)]+=\$4/\$5}}END{printf (\"plus\");for (i=20;i<=30;i++){if(aplus[i]==\"\"){printf (\"\\t0\")}else{printf (\"\\t\"aplus[i])}};printf(\"\\n\"); printf (\"minus\");for (i=20;i<=30;i++){if(aminus[i]==\"\"){printf (\"\\t0\")}else{printf (\"\\t-\"aminus[i])}};printf(\"\\n\")}' $i.24suste_w_protop_subtracted.tmp $i.nonredste.tmp >  $i.nonredste_multi_lendist;
awk 'BEGIN{FS=OFS=\"\\t\"}(FNR==NR){found[\$7]=1}(FNR<NR)&&(found[\$7]==\"\"){if(\$6==\$13){ aplus[(\$3-\$2)]+=\$4/\$5}else{aminus[(\$3-\$2)]+=\$4/\$5}}END{printf (\"plus\");for (i=20;i<=30;i++){if(aplus[i]==\"\"){printf (\"\\t0\")}else{printf (\"\\t\"aplus[i])}};printf(\"\\n\"); printf (\"minus\");for (i=20;i<=30;i++){if(aminus[i]==\"\"){printf (\"\\t0\")}else{printf (\"\\t-\"aminus[i])}};printf(\"\\n\")}' $i.nonredste.tmp $i.24suste_w_protop_subtracted.tmp >  $i.24suste_w_protop_subtracted_uniq_lendist;
awk 'BEGIN{FS=OFS=\"\\t\"}(FNR==NR){found[\$7]=1}(FNR<NR)&&(found[\$7]==1){if(\$6==\$13){ aplus[(\$3-\$2)]+=\$4/\$5}else{aminus[(\$3-\$2)]+=\$4/\$5}}END{printf (\"plus\");for (i=20;i<=30;i++){if(aplus[i]==\"\"){printf (\"\\t0\")}else{printf (\"\\t\"aplus[i])}};printf(\"\\n\"); printf (\"minus\");for (i=20;i<=30;i++){if(aminus[i]==\"\"){printf (\"\\t0\")}else{printf (\"\\t-\"aminus[i])}};printf(\"\\n\")}' $i.nonredste.tmp $i.24suste_w_protop_subtracted.tmp >  $i.24suste_w_protop_subtracted_multi_lendist"; done

#11a.get all suste+1360 piRNAs for downstream 5P-seq data analyses
module load bedtools/2.30.0
for i in $(ls *.20.realbed2); do /pi/phillip.zamore-umw/home/ildar.gainetdinov-umw/common/pipelines/ms_short_8G "bedtools intersect -wa -a $i -b 24suste_protop.bed | sort | uniq > $i.suste_protop"; done

#11b.get all piRNAs from 42AB, petrel, nos, bam, and bgcn for downstream 5P-seq data analyses
module load bedtools/2.30.0
for i in $(ls *.20.realbed2); do /pi/phillip.zamore-umw/home/ildar.gainetdinov-umw/common/pipelines/ms_short_8G "bedtools intersect -wa -a $i -b top_testis_clusters_plus_genes.bed | sort | uniq > $i.other"; done

