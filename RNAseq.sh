#1.reformat umi sequences
module load python2/2.7.9
for i in `find . -maxdepth 1 -name "*R1.fastq"`; do bsub -q short -W 240 -n 4 -o $i.out -e $i.err "python /pi/phillip.zamore-umw/home/ildar.gainetdinov-umw/common/pipelines/RNAseq_reformat_umi_fastq -l $i -r ${i/R1/R2} -L $i.reformated  -R ${i/R1/R2}.reformated"; done

#2.run piPpipes to remove rRNA and align with STAR
for i in `find . -maxdepth 1 -name "*R1.fastq.reformated"`; do /pi/phillip.zamore-umw/home/ildar.gainetdinov-umw/common/pipelines/ms "/pi/phillip.zamore-umw/home/ildar.gainetdinov-umw/common/pipelines/piPipes/piPipes rna -l $i -r ${i/R1/R2} -g dm6 -o $i.pipipes -c 8"; done

#3.move bam files from piPipes to current directory
for i in $(find . -maxdepth 3 -name "*.dm6.sorted.bam"); do mv $i .; done

#4a.in bam files from piPipes mark duplicates
module load python2/2.7.9
module load gcc/12.2.0
pip install --user pysam
for i in $(find . -maxdepth 1 -name "*.dm6.sorted.bam"); do /pi/phillip.zamore-umw/home/ildar.gainetdinov-umw/common/pipelines/ms  "python2 /pi/phillip.zamore-umw/home/ildar.gainetdinov-umw/common/pipelines/RNAseq_umi_mark_duplicates.py -f $i -p 8"; done

#4b.deduplicate bams
module load samtools/1.16.1
for i in $(find . -maxdepth 1 -name "*.deumi.sorted.bam"); do /pi/phillip.zamore-umw/home/ildar.gainetdinov-umw/common/pipelines/ms_short_64G "samtools view  -@ 8 -b -F 0x400 $i > $i.dedup"; done

#4c.sort deduplicated bams by chrom pos
for i in $(find . -maxdepth 1 -name "*.deumi.sorted.bam.dedup"); do /pi/phillip.zamore-umw/home/ildar.gainetdinov-umw/common/pipelines/ms_short_64G "samtools sort -@ 8 $i > $i.sorted.bam"; done

#5.calculate the number of all mapped reads from .bam
module load bedtools/2.26.0
for i in $(find . -maxdepth 1 -name "*.bam"); do /pi/phillip.zamore-umw/home/ildar.gainetdinov-umw/common/pipelines/ms_short_64G "bedtools bamtobed -bed12 -tag NH -i $i | awk -v f=$i '{a[\$4]++}END{for (j in a){total++};print f\"\\t\"total/2}' > $i.all.mapped.reads.txt"; done

#6.get bed12 files and intersect those with with Ste and find Su(Ste) only
#6a.get bed12
module load bedtools/2.30.0
for i in $(find . -name "*.deumi.sorted.bam.dedup.sorted.bam"); do /pi/phillip.zamore-umw/home/ildar.gainetdinov-umw/common/pipelines/ms_short_64G "bedtools bamtobed -bed12 -tag NH -i $i | awk 'BEGIN{FS=OFS=\"\\t\"}{if(substr(\$4,length(\$4))==1){\$6=(\$6==\"+\"?\"-\":\"+\")}; print \$0}' > $i.all.bed12; awk '{a[\$4]++}END{for (j in a){total++};print total/2}' $i.all.bed12 > $i.all.reads"; done

#6b.intersect with Ste and find Su(Ste) only
for i in $(find . -name "*.bed12"); do /pi/phillip.zamore-umw/home/ildar.gainetdinov-umw/common/pipelines/ms_short_64G "bedtools intersect -wo -a $i -b 24suste_w_protop_subtracted_n_ste.bed > $i.24"; done

#6d.metaplot: for ste - only keep gaps < 60-nt and for suste - no gaps 
for i in $(find . -name "*.bed12.??"); do /pi/phillip.zamore-umw/home/ildar.gainetdinov-umw/common/pipelines/ms_short_64G "awk 'BEGIN{FS=OFS=\"\\t\"}(FNR==NR)&&(\$16==\"na\"){n=split(\$12,fl,\",\");if(fl[n]<=60){ste[substr(\$4,1,length(\$4)-2)]=1}}(FNR<NR)&&(\$16!=\"na\")&&(ste[substr(\$4,1,length(\$4)-2)]!=1){if(\$10==1){if(\$6==\$18){same[substr(\$4,1,length(\$4)-2)]=\$2\";\"\$3\";\"\$5\";\"\$6\";\"\$11\";\"\$14\";\"\$15}else{opposite[substr(\$4,1,length(\$4)-2)]=\$2\";\"\$3\";\"\$5\";\"\$6\";\"\$11\";\"\$14\";\"\$15}}}END{for(read in same){if (opposite[read]==\"\"){split(same[read],fl,\";\");if(fl[4]==\"+\"){for(i=0;i<fl[5];i++){sametotal[(fl[1]-fl[6]+i)]+=1/fl[3]}};if(fl[4]==\"-\"){for(i=0;i<fl[5];i++){sametotal[(fl[7]-fl[2]+i)]+=1/fl[3]}}}};for(read in opposite){if (same[read]==\"\"){split(opposite[read],fl,\";\");if(fl[4]==\"-\"){for(i=0;i<fl[5];i++){oppositetotal[(fl[1]-fl[6]+i)]+=1/fl[3]}};if(fl[4]==\"+\"){for(i=0;i<fl[5];i++){oppositetotal[(fl[7]-fl[2]+i)]+=1/fl[3]}}}};for (i=-1000;i<=3000;i++){print i\"\\t\"(1*sametotal[i])\"\\t\" (-1*oppositetotal[i])}}' $i $i > $i.meta.txt"; done

paste *dpp_minus*.24.meta.txt > meta_stranded_RSQ-seq_suste_only_dpp_minus.24.txt