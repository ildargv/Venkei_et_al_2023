#1.find and remove 5' UMI part
for i in `find . -maxdepth 1 -name "*.R1.fastq"`; do /pi/phillip.zamore-umw/home/ildar.gainetdinov-umw/common/pipelines/ms_short_64G "awk 'BEGIN{FS=OFS=\"\\t\"}(FNR%4==1)&&(FNR==NR){name=\$1;getline;if((index(\$1,\"N\")==0)&&(((substr(\$1,4,3)==\"ATC\")&&(substr(\$1,10,3)==\"AGT\"))||((substr(\$1,4,3)==\"CGA\")&&(substr(\$1,10,3)==\"TAC\")))){a[FNR-1]=1; split(name,b,\" \"); c[FNR-1]=b[1]\"_\"substr(\$1,1,15);
print b[1]\"_\"substr(\$1,1,15)\" \"b[2]> FILENAME\".rfmtd.fq\"; 
print substr(\$1,16)> FILENAME\".rfmtd.fq\"; getline;
print > FILENAME\".rfmtd.fq\"; getline
print substr(\$1,16)> FILENAME\".rfmtd.fq\"}}(FNR<NR)&&(a[FNR]==1)&&(FNR%4==1){split(\$1,b,\" \"); print c[FNR]\" \"b[2] > FILENAME\".rfmtd.fq\" ; getline;print > FILENAME\".rfmtd.fq\";getline;print > FILENAME\".rfmtd.fq\";getline;print > FILENAME\".rfmtd.fq\"}' $i ${i/.R1./.R2.}"; done

#2.remove rRNA reads and align onto mm10 with STAR in piPipes
for i in `find . -name "*.R1.fastq.rfmtd.fq"`; do /pi/phillip.zamore-umw/home/ildar.gainetdinov-umw/common/pipelines/ms "/pi/phillip.zamore-umw/home/ildar.gainetdinov-umw/common/pipelines/piPipes/piPipes deg -l $i -r ${i/R1/R2} -g dm6 -o $i.pipipes.deg -c 8"; done

#3.collapse 5ends for all apportioned mapped reads
for i in $(find . -maxdepth 1  -name "*.all.bed12"); do /pi/phillip.zamore-umw/home/ildar.gainetdinov-umw/common/pipelines/ms_short_64G "awk '{if (\$6==\"+\"){aplus[\$1,\$2]+=\$5}else{aminus[\$1,\$3]+=\$5}}END{for (ij in aplus) {split(ij,indices,SUBSEP);
 i=indices[1];
 j=indices[2]; print i\"\\t\"j\"\\t\"j+1\"\\t\"aplus[i,j]\"\\tna\\t+\"};for (ij in aminus) {split(ij,indices,SUBSEP);
 i=indices[1];
 j=indices[2]; print i\"\\t\"j-1\"\\t\"j\"\\t\"aminus[i,j]\"\\tna\\t-\"}}' $i > $i.1"; done

#4.normalize to seq.depth: make rpm files
for i in $(find . -name "*.1"); do /pi/phillip.zamore-umw/home/ildar.gainetdinov-umw/common/pipelines/ms_short_64G "awk '{a[FNR]=\$0; total+=\$4; num=FNR}END{for (i=1;i<=num;i++){split (a[i],b,\"\\t\");print b[1]\"\\t\"b[2]\"\\t\"b[3]\"\\t\"1000000*b[4]/total\"\\t\"b[4]\"\\t\"b[6]}}' $i > $i.rpm"; done

#5a.select reads from 1360 repeats
module load bedtools/2.30.0
for i in $(ls *.x_rRNA.dm6.sorted.f0x40.noS.all.bed12.1.rpm); do /pi/phillip.zamore-umw/home/ildar.gainetdinov-umw/common/pipelines/ms_short_8G "bedtools intersect -s -wo -a $i -b 24protop100.bed | awk 'BEGIN{FS=OFS=\"\\t\"}(\$12==\"+\"){print \$0\"\\t\"\$9-\$2}(\$12==\"-\"){print \$0\"\\t\"\$3-\$8}' > $i.24protop100"; done
bed
#5b.get ping-pong frequencies between 5P-seq reads from 1360 repeats and smallRNA-seq reads
for srs in $(ls *.20.realbed2.suste_protop); do for deg in $(ls *.24protop100); do echo $deg;
awk '(FNR==NR){if($6=="+"){aplus[$1,$2]+=1}else{aminus[$1,$3]+=1}}(FNR<NR){
if($6=="+"){for (k=0;k<=20;k++){res[k]+=aminus[$1,$2+k]}};if($6=="-"){for (k=0;k<=20;k++){res[k]+=aplus[$1,$3-k]}}}END{
 for (i=0;i<=20;i++){print i"\t"res[i]}}' $srs $deg > ${deg/.x_rRNA.dm6.sorted.f0x40.noS.all.bed12.1.rpm.24protop100/.24}.${srs/.fastq.trimmed.q20p100.deUMI.dedup.fq.x_filter_1mm.x_spikein.p20.bed2.rpm.20.realbed2.suste_protop/.20}.speciesunap.pp0to20;
 awk '(FNR==NR){if($6=="+"){aplus[$1,$2]+=1/$5}else{aminus[$1,$3]+=1/$5}}(FNR<NR){
if($6=="+"){for (k=0;k<=20;k++){res[k]+=aminus[$1,$2+k]}};if($6=="-"){for (k=0;k<=20;k++){res[k]+=aplus[$1,$3-k]}}}END{
 for (i=0;i<=20;i++){print i"\t"res[i]}}' $srs $deg > ${deg/.x_rRNA.dm6.sorted.f0x40.noS.all.bed12.1.rpm.24protop100/.24}.${srs/.fastq.trimmed.q20p100.deUMI.dedup.fq.x_filter_1mm.x_spikein.p20.bed2.rpm.20.realbed2.suste_protop/.20}.species.pp0to20;
awk '(FNR==NR){if($6=="+"){aplus[$1,$2]+=$4/$5}else{aminus[$1,$3]+=$4/$5}}(FNR<NR){
if($6=="+"){for (k=0;k<=20;k++){res[k]+=aminus[$1,$2+k]}};if($6=="-"){for (k=0;k<=20;k++){res[k]+=aplus[$1,$3-k]}}}END{
 for (i=0;i<=20;i++){print i"\t"res[i]}}' $srs $deg > ${deg/.x_rRNA.dm6.sorted.f0x40.noS.all.bed12.1.rpm.24protop100/.24}.${srs/.fastq.trimmed.q20p100.deUMI.dedup.fq.x_filter_1mm.x_spikein.p20.bed2.rpm.20.realbed2.suste_protop/.20}.reads.pp0to20; done; done

#5c.select 5P-seq reads explained by smallRNA-seq reads with 10nt overlap on opposite genomic strand
for srs in $(ls *.20.realbed2.suste_protop); do for deg in $(ls *.24protop100); do echo $deg; awk '(FNR==NR){f1=FILENAME;if($6=="+"){ aplus[$1,$2]=1}else{aminus[$1,$3]=1}}(FNR<NR){f2=FILENAME;total+=$4;if(($6=="+")&&(aminus[$1,($2+10)]==1)){print;pp+=$4}; if(($6=="-")&&(aplus[$1,($3-10)]==1)){print;pp+=$4}}END{print "24\t20\t"f1"\t"f2"\t"pp/total"\t"pp"\t"total >> "pp.fraction.txt"}' $srs $deg > ${deg/.x_rRNA.dm6.sorted.f0x40.noS.all.bed12.1.rpm.24protop100/.24}.${srs/.fastq.trimmed.q20p100.deUMI.dedup.fq.x_filter_1mm.x_spikein.p20.bed2.rpm.20.realbed2.suste_protop/.20}.pp; done; done

#5d.get metaplots for distribution of small RNA reads around 5P-seq reads on the same genomic strand explained by smallRNA-seq reads with 10nt overlap on opposite genomic strand
for srs in $(ls *.20.realbed2.suste_protop); do for deg in $(ls *.20.pp); do echo $deg;awk '(FNR==NR){gsub(".fastq.trimmed.q20p100.deUMI.dedup.fq.x_filter_1mm.x_spikein.p20.bed2.rpm.20.realbed2.suste_protop","",FILENAME);f1=FILENAME;if($6=="+"){ aplus[$1,$2]++}else{aminus[$1,$3]++}}(FNR<NR){f2=FILENAME;if($6=="+"){for (k=-100; k<=500; k++){res[k]+=aplus[$1,$2+k]}};if($6=="-"){for (k=-100; k<=500; k++){res[k]+=aminus[$1,$3-k]}}}END{printf ("24\t20\t"f1"\t"f2); for (i=-100;i<=500;i++){printf "\t"res[i]};printf("\n")}' $srs $deg >> srs_around_deg.txt; done; done;

#6a.select reads from 42AB, petrel, nos, bam, and bgcn
module load bedtools/2.30.0
for i in $(ls *.x_rRNA.dm6.sorted.f0x40.noS.all.bed12.1.rpm); do /pi/phillip.zamore-umw/home/ildar.gainetdinov-umw/common/pipelines/ms_short_8G "bedtools intersect -wo -a $i -b top_testis_clusters_plus_genes.bed | awk 'BEGIN{FS=OFS=\"\\t\"}(\$12==\"+\"){print \$0\"\\t\"\$9-\$2}(\$12==\"-\"){print \$0\"\\t\"\$3-\$8}' > $i.others"; done

#6b.calculate total abundance of 5P-seq reads from 42AB, petrel, nos, bam, and bgcn
for deg in $(ls *.others); do echo $deg;
awk '($11=="ss")&&($6==$12){total[$10]+=$4}($11=="na"){total[$10]+=$4}END{
 for (gene in total){print gene"\t"total[gene]}}' $deg > ${deg/.x_rRNA.dm6.sorted.f0x40.noS.all.bed12.1.rpm/}.sum; done

#6c.get ping-pong frequencies between 5P-seq reads from 42AB, petrel, nos, bam, and bgcn and smallRNA-seq reads
for srs in $(ls SRS_UMI_testes_XXmom_c*.20.realbed2.other); do for deg in $(ls DEG_UMI_testes_XXmom_*.others ID?_IG_deg_dpp_*.others); do echo $deg;
awk '(FNR==NR){if($6=="+"){aplus[$1,$2]+=$4/$5}else{aminus[$1,$3]+=$4/$5}}((FNR<NR)&&($11=="ss")&&($6==$12))||((FNR<NR)&&($11=="na")){
genes[$10]=1;if($6=="+"){for (k=0;k<=20;k++){res[$10,k]+=aminus[$1,$2+k]}};if($6=="-"){for (k=0;k<=20;k++){res[$10,k]+=aplus[$1,$3-k]}}}END{
 for(gene in genes){printf(gene"\treads");for (i=0;i<=20;i++){printf ("\t"res[gene,i])};printf ("\n")}}' $srs $deg > ${deg/.x_rRNA.dm6.sorted.f0x40.noS.all.bed12.1.rpm.others/}.${srs/.fastq.trimmed.q20p100.deUMI.dedup.fq.x_filter_1mm.x_spikein.p20.bed2.rpm.20.realbed2.other/}.ctrls; awk '(FNR==NR){if($6=="+"){aplus[$1,$2]=1}else{aminus[$1,$3]=1}}((FNR<NR)&&($11=="ss")&&($6==$12))||((FNR<NR)&&($11=="na")){
genes[$10]=1;if($6=="+"){for (k=0;k<=20;k++){res[$10,k]+=aminus[$1,$2+k]}};if($6=="-"){for (k=0;k<=20;k++){res[$10,k]+=aplus[$1,$3-k]}}}END{
 for(gene in genes){printf(gene"\tspecies");for (i=0;i<=20;i++){printf ("\t"res[gene,i])};printf ("\n")}}' $srs $deg >> ${deg/.x_rRNA.dm6.sorted.f0x40.noS.all.bed12.1.rpm.others/}.${srs/.fastq.trimmed.q20p100.deUMI.dedup.fq.x_filter_1mm.x_spikein.p20.bed2.rpm.20.realbed2.other/}.ctrls; done; done
